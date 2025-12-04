
"""
根据 GWH assemblyZhang API：
- 对每个 accession，优先使用 ftpPathRefseqGff 下载重注释 GFF（.gff.gz）；
- 若没有 ftpPathRefseqGff，则使用 ftpPathGff；
- 解压后，从 GFF 中提取“调控元件”相关注释行，写成新的 GFF3 文件。
"""

import os
import io
import gzip
import time
import argparse
from typing import Iterable, List, Optional

import requests

BASE_API = "https://ngdc.cncb.ac.cn/gwh/api/public/assemblyZhang"

# 常见的调控相关 feature type，可按需要增删
REGULATORY_TYPES = {
    "regulatory_region",
    "promoter",
    "enhancer",
    "silencer",
    "insulator",
    "TF_binding_site",
    "CTCF_binding_site",
    "DNaseI_hypersensitive_site",
    "promoter_flanking_region",
    "open_chromatin_region",
    "tss_region",
}


def fetch_metadata(accession: str) -> dict:
    """调用 assemblyZhang API 获取元数据 JSON。"""
    url = f"{BASE_API}/{accession}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.json()


def choose_gff_url(record: dict) -> Optional[str]:
    """
    根据 JSON 记录选择 GFF 下载 URL：
    1. 若 ftpPathRefseqGff 有值，优先使用；
    2. 否则若 ftpPathGff 有值，则使用；
    3. 否则返回 None。
    """
    refseq = record.get("ftpPathRefseqGff")
    if refseq:
        return refseq

    orig = record.get("ftpPathGff")
    if orig:
        return orig

    return None


def download_gff_text(gff_url: str) -> str:
    """
    下载 GFF 文件（.gff 或 .gff.gz）。
    如果是 .gz，自动解压；返回 UTF-8 文本。
    """
    resp = requests.get(gff_url, timeout=120)
    resp.raise_for_status()
    data = resp.content

    if gff_url.endswith(".gz"):
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
            text = gz.read().decode("utf-8")
    else:
        text = data.decode("utf-8")

    return text


def is_regulatory_type(feature_type: str) -> bool:
    """
    判断一个 GFF 第 3 列 feature type 是否属于“调控元件”。
    可根据实际数据调整。
    """
    if feature_type in REGULATORY_TYPES:
        return True

    # 宽松规则：包含 “regulatory” 的 type 也视作调控相关
    if "regulatory" in feature_type.lower():
        return True

    return False


def extract_regulatory_gff(gff_text: str) -> Optional[str]:
    """
    从原始 GFF 文本中提取调控元件注释：
    - 保留所有以 '#' 开头的头部行；
    - 数据行中，只保留 feature_type 属于调控类型的记录。
    若没有任何调控元件行，返回 None。
    """
    headers: List[str] = []
    reg_lines: List[str] = []

    for line in gff_text.splitlines():
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith("#"):
            headers.append(line)
            continue

        fields = line.split("\t")
        if len(fields) < 3:
            continue

        feature_type = fields[2]
        if is_regulatory_type(feature_type):
            reg_lines.append(line)

    if not reg_lines:
        return None

    # 如果没有 gff-version 头，补上一行
    if not any(h.startswith("##gff-version") for h in headers):
        headers.insert(0, "##gff-version 3")

    return "\n".join(headers + reg_lines) + "\n"


def iter_accessions_from_file(path: str) -> Iterable[str]:
    """从文本文件中逐行读取 accession，跳过空行和以 # 开头的注释行。"""
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            yield line


def process_accession(
    accession: str,
    out_dir: str,
    sleep_sec: float = 0.2,
) -> None:
    """
    处理单个 accession：
    - 获取 metadata；
    - 选择 GFF URL（重注释优先，其次原始 GFF）；
    - 下载 & 解压；
    - 提取调控元件注释；
    - 写出 <accession>_regulatory.gff3。
    """
    print(f"[{accession}] 正在获取元数据...")
    try:
        record = fetch_metadata(accession)
    except Exception as e:
        print(f"[{accession}] 错误：获取元数据失败：{e}")
        return

    gff_url = choose_gff_url(record)
    if not gff_url:
        print(f"[{accession}] 未找到 ftpPathRefseqGff 或 ftpPathGff，跳过。")
        return

    print(f"[{accession}] 正在下载 GFF 文件：{gff_url}")
    try:
        gff_text = download_gff_text(gff_url)
    except Exception as e:
        print(f"[{accession}] 错误：下载或解析 GFF 文件失败：{e}")
        return

    filtered = extract_regulatory_gff(gff_text)
    if filtered is None:
        print(f"[{accession}] 未在 GFF 中找到任何调控元件注释，跳过写出。")
        return

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{accession}_regulatory.gff3")
    with open(out_path, "w", encoding="utf-8") as out_f:
        out_f.write(filtered)

    print(f"[{accession}] 已将调控元件注释写入：{out_path}")

    if sleep_sec > 0:
        time.sleep(sleep_sec)


def main():
    parser = argparse.ArgumentParser(
        description="从 GWH GFF 文件中提取调控元件注释（优先使用重注释 GFF）。"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="包含 accession 列表的文本文件，每行一个。",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="输出目录，用于保存 *_regulatory.gff3 文件。",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.2,
        help="每个 accession 之间暂停的秒数（默认 0.2）。",
    )

    args = parser.parse_args()

    for acc in iter_accessions_from_file(args.input):
        process_accession(acc, args.outdir, sleep_sec=args.sleep)


if __name__ == "__main__":
    main()
