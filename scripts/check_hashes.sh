#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root"

if command -v shasum >/dev/null 2>&1; then
  shasum -a 256 -c SHA256SUMS
elif command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c SHA256SUMS
else
  echo "Neither shasum nor sha256sum is available." >&2
  exit 1
fi
