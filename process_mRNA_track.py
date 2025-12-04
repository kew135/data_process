# 根据API中mRnaAnnotation字段的值是否大于0，来判断当前这一条记录的非重注释的GFF文件中是否有mRNA注释
# 如果有，对该GFF文件提取出关于mRNA的注释
# 提取完成后还以GFF格式存储


import os
import io
import gzip
import time
import argparse
from typing import Iterable, List, Optional

import requests

BASE_API = "https://ngdc.cncb.ac.cn/gwh/api/public/assemblyZhang"


# 将 JSON 中的 mRnaAnnotation 字段转为 int
def parse_mrna_annotation(value) -> int:
    if value is None:
        return 0
    if isinstance(value, int):
        return value
    s = str(value).strip()
    if not s:
        return 0
    s = s.replace(",", "")  # 去掉千分位逗号
    try:
        return int(s)
    except ValueError:
        return 0

# 调用 API 获取元数据 JSON
def fetch_metadata(accession: str) -> dict:
    url = f"{BASE_API}/{accession}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.json()

# 从 JSON 记录中获取 GFF 文件的下载 URL(ftpPathGff 字段)
def get_gff_url(record: dict) -> Optional[str]:
    url = record.get("ftpPathGff")
    if not url:
        return None
    return url

# 下载 GFF 文件内容，支持 .gz 格式解压
def download_gff_text(gff_url: str) -> str:
    resp = requests.get(gff_url, timeout=120)
    resp.raise_for_status()
    data = resp.content

    if gff_url.endswith(".gz"):
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
            text = gz.read().decode("utf-8")
    else:
        text = data.decode("utf-8")

    return text


def extract_mrna_gff(gff_text: str) -> Optional[str]:
    """
    从原始 GFF 文本中提取 mRNA 注释：
    - 保留所有以 # 开头的头部行；
    - 数据行中，只保留第三列(type) == "mRNA" 的行。
    若没有任何 mRNA 行，返回 None。
    """
    headers: List[str] = []
    mrna_lines: List[str] = []

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
        if feature_type == "mRNA":
            mrna_lines.append(line)

    if not mrna_lines:
        return None

    # 确保有 gff-version 头，没有的话补一行
    if not any(h.startswith("##gff-version") for h in headers):
        headers.insert(0, "##gff-version 3")

    return "\n".join(headers + mrna_lines) + "\n"


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
    1. 取 metadata；
    2. 若 mRnaAnnotation > 0，则下载 ftpPathGff；
    3. 提取 mRNA 注释，写到 <accession>_mRNA.gff3。
    """
    print(f"[{accession}] Fetching metadata...")
    try:
        record = fetch_metadata(accession)
    except Exception as e:
        print(f"[{accession}] ERROR: failed to fetch metadata: {e}")
        return

    mrna_count = parse_mrna_annotation(record.get("mRnaAnnotation"))
    if mrna_count <= 0:
        print(f"[{accession}] mRnaAnnotation={mrna_count}, skip.")
        return

    gff_url = get_gff_url(record)
    if not gff_url:
        print(f"[{accession}] No ftpPathGff field, skip.")
        return

    print(f"[{accession}] mRnaAnnotation={mrna_count}, downloading GFF: {gff_url}")
    try:
        gff_text = download_gff_text(gff_url)
    except Exception as e:
        print(f"[{accession}] ERROR: failed to download/parse GFF: {e}")
        return

    filtered = extract_mrna_gff(gff_text)
    if filtered is None:
        print(f"[{accession}] No mRNA features found in GFF, skip writing.")
        return

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{accession}_mRNA.gff3")
    with open(out_path, "w", encoding="utf-8") as out_f:
        out_f.write(filtered)

    print(f"[{accession}] Wrote mRNA-only GFF to: {out_path}")

    if sleep_sec > 0:
        time.sleep(sleep_sec)


def main():
    parser = argparse.ArgumentParser(
        description="提取 mRNA 注释，以GFF文件存储"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="包含 accession 列表的文本文件，每行一个",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="输出目录，用于保存 *_mRNA.gff3 文件",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.2,
        help="每个 accession 之间暂停秒数（默认 0.2）",
    )

    args = parser.parse_args()

    for acc in iter_accessions_from_file(args.input):
        process_accession(acc, args.outdir, sleep_sec=args.sleep)


if __name__ == "__main__":
    main()
