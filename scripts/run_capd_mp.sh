#!/usr/bin/env bash
set -euo pipefail

REQUIRED_CAPD_VERSION="5.3.0"

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repo_root/capd_cpp"

usage() {
  cat <<'EOF'
Compile and run the complete 220-bit C++ CAPD numerical certificate.

Usage:
  bash scripts/run_capd_mp.sh

Run this command from the repository root.  If CAPD 5.3.0 was installed in a
nonstandard location, set MPCAPD_CONFIG to the absolute path of its executable
mpcapd-config file.
EOF
}

if [[ $# -eq 1 && ( "$1" == "--help" || "$1" == "-h" ) ]]; then
  usage
  exit 0
fi
if [[ $# -ne 0 ]]; then
  usage >&2
  exit 2
fi

sha256_file() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    echo "Neither shasum nor sha256sum is available." >&2
    exit 1
  fi
}

if grep -nE 'RF\(|std::stod|nextafter|I[[:space:]]*\([[:space:]]*R[[:space:]]*\("' src/lppssi_cert_v6.cpp; then
  echo "Unsafe scalar parsing of a decimal literal found in the certificate source." >&2
  exit 1
fi

if [[ -n "${LPPSSI_EXPORT_VORTEX_TUBE:-}" ||
      -n "${LPPSSI_EXPORT_VORTEX_TUBE_CELLS:-}" ]]; then
  echo "Unset the vortex-profile export variables before running the certificate." >&2
  exit 1
fi

if [[ -z "${MPCAPD_CONFIG:-}" ]]; then
  if command -v mpcapd-config >/dev/null 2>&1; then
    MPCAPD_CONFIG="$(command -v mpcapd-config)"
  elif [[ -x "$HOME/.local/share/aymh-capd/build-capd-5.3.0-mp/bin/mpcapd-config" ]]; then
    MPCAPD_CONFIG="$HOME/.local/share/aymh-capd/build-capd-5.3.0-mp/bin/mpcapd-config"
  else
    cat >&2 <<'EOF'
The CAPD multiprecision configuration program mpcapd-config was not found.

For a supported macOS, Debian, Ubuntu, or WSL2 installation, return to the
repository root and run

  bash scripts/install_capd_5_3_0.sh

If CAPD is already installed elsewhere, run this script as

  MPCAPD_CONFIG=/absolute/path/to/mpcapd-config bash scripts/run_capd_mp.sh
EOF
    exit 1
  fi
fi
export MPCAPD_CONFIG

if [[ ! -x "$MPCAPD_CONFIG" ]]; then
  echo "MPCAPD_CONFIG is not an executable file: $MPCAPD_CONFIG" >&2
  exit 1
fi
cflags="$("$MPCAPD_CONFIG" --cflags)"
link_flags="$("$MPCAPD_CONFIG" --libs)"
if [[ -z "$cflags" || -z "$link_flags" ]]; then
  echo "mpcapd-config did not return usable compiler and linker flags." >&2
  exit 1
fi
if [[ "$link_flags" != *mpfr* || "$link_flags" != *gmp* ]]; then
  echo "The selected CAPD configuration does not report the MPFR/GMP libraries." >&2
  echo "Returned flags: $link_flags" >&2
  exit 1
fi
capd_version="$("$MPCAPD_CONFIG" --modversion 2>/dev/null || printf 'unknown')"
if [[ "$capd_version" != "$REQUIRED_CAPD_VERSION" ]]; then
  echo "The selected CAPD version is $capd_version; CAPD $REQUIRED_CAPD_VERSION is required." >&2
  echo "Selected configuration program: $MPCAPD_CONFIG" >&2
  exit 1
fi

echo "Clean-building C++ CAPD certificate"
echo "MPCAPD_CONFIG=${MPCAPD_CONFIG}"
timestamp="$(date +%Y%m%d_%H%M%S)"
build_log="capd_build_v6_mp_local_${timestamp}.txt"
log="capd_certificate_v6_220bit_local_${timestamp}.txt"
if ! {
  make clean
  make mp
} > "$build_log" 2>&1; then
  echo "The certificate build failed.  The last 40 lines are:" >&2
  tail -n 40 "$build_log" >&2
  echo "Complete build log: capd_cpp/$build_log" >&2
  exit 1
fi
echo "Certificate build completed"
echo "Build log: capd_cpp/$build_log"

echo "Running C++ CAPD certificate"
source_sha="$(sha256_file src/lppssi_cert_v6.cpp)"
compiler="$(make -s print-cxx)"
compiler_version="$("$compiler" --version 2>/dev/null | sed -n '1p')"
{
  echo "certificate source SHA-256 = ${source_sha}"
  echo "CAPD version = ${capd_version}"
  echo "compiler = ${compiler_version}"
  ./lppssi_cert_v6_mp
} 2>&1 | tee "$log"

if ! grep -Fq "interval backend = MpInterval/MpFloat, requested precision bits = 220" "$log"; then
  echo "The run did not report the required 220-bit multiprecision backend." >&2
  exit 1
fi
if ! grep -Fq "Taylor order = 40" "$log"; then
  echo "The run did not report the required Taylor order 40." >&2
  exit 1
fi
if ! grep -Fq "OVERALL: PASS FOR THE NUMERICAL CERTIFICATE BLOCKS" "$log"; then
  echo "The numerical certificate did not finish with OVERALL: PASS." >&2
  exit 1
fi

echo
echo "Local log: capd_cpp/${log}"
printf 'Local log SHA-256: %s\n' "$(sha256_file "$log")"
