#!/usr/bin/env bash
# Build the search engine. Run from this directory (Linux, macOS or WSL).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p bin
gcc -O3 -march=native -o bin/sparse3 sparse3.c -lm
echo "built bin/sparse3"
