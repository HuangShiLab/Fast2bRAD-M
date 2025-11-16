#!/usr/bin/env bash
set -euo pipefail

# 用法示例：
# scripts/run_pipeline.sh \
#   --samples /abs/path/samples.tsv \
#   --site BcgI \
#   --level species \
#   --outdir /abs/path/runs/run1 \
#   --prefix run1 \
#   --threads 8 \
#   --resume yes \
#   [--genome-list /abs/path/genomes.tsv] \
#   [--pre-digested-dir /abs/path/predig] \
#   [--database /abs/path/db_ready] \
#   [--mock mock1,mock2] \
#   [--control ctrl1,ctrl2]

BIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${BIN_DIR}/target/release/fast2bRAD-M"

if [[ ! -x "$BIN" ]]; then
  echo "未找到可执行文件：$BIN"
  echo "请先构建：cargo build --release"
  exit 1
fi

exec "$BIN" pipeline "$@"


