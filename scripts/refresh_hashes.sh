#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root"

if command -v shasum >/dev/null 2>&1; then
  sha256=(shasum -a 256)
elif command -v sha256sum >/dev/null 2>&1; then
  sha256=(sha256sum)
else
  echo "Neither shasum nor sha256sum is available." >&2
  exit 1
fi

find . -type f \
  ! -path './.git/*' \
  ! -path './capd_cpp/lppssi_cert_v6_mp' \
  ! -path './capd_cpp/capd_build_v6_mp_local_*.txt' \
  ! -path './capd_cpp/capd_certificate_v6_220bit_local_*.txt' \
  ! -path './capd_cpp/build/*' \
  ! -name '.DS_Store' \
  ! -name 'SHA256SUMS' \
  | sort \
  | sed 's#^\./##' \
  | xargs "${sha256[@]}" > SHA256SUMS
