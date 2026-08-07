#!/usr/bin/env bash
set -euo pipefail

CAPD_VERSION="5.3.0"
CAPD_ARCHIVE="capd-${CAPD_VERSION}.tar.gz"
CAPD_URL="https://sourceforge.net/projects/capd/files/${CAPD_VERSION}/${CAPD_ARCHIVE}/download"
CAPD_SHA256="e4100959a5409d330f8907d050f101a0485489075b4ce0d5eb2e349a2f8bf228"

CAPD_HOME="${CAPD_HOME:-$HOME/.local/share/aymh-capd}"
CAPD_PREFIX="${CAPD_PREFIX:-$CAPD_HOME/install-prefix}"
DOWNLOAD_DIR="$CAPD_HOME/downloads"
SOURCE_DIR="$CAPD_HOME/capd-${CAPD_VERSION}"
BUILD_DIR="$CAPD_HOME/build-capd-${CAPD_VERSION}-mp"
LOG_DIR="$CAPD_HOME/logs"
ARCHIVE_PATH="$DOWNLOAD_DIR/$CAPD_ARCHIVE"
CONFIG_STAMP="$BUILD_DIR/.aymh_configure_complete"

usage() {
  cat <<EOF
Install CAPD ${CAPD_VERSION} with MPFR multiprecision interval arithmetic.

Usage:
  $0

Supported systems:
  macOS with Homebrew
  Debian, Ubuntu, or Ubuntu under WSL2

The CAPD source and build directories are placed under
  $CAPD_HOME

Advanced users may change the working directory by setting CAPD_HOME before
invoking this script.  Set AYMH_SKIP_SYSTEM_PACKAGES=1 to skip the
operating-system package installation.  Set CAPD_BUILD_JOBS to limit the
number of parallel compilation jobs.  Under Slurm, the script automatically
uses SLURM_CPUS_PER_TASK.
EOF
}

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--help" && "$1" != "-h" ) ]]; then
  usage >&2
  exit 2
fi
if [[ $# -eq 1 ]]; then
  usage
  exit 0
fi

say() {
  printf '\n%s\n' "$*"
}

have() {
  command -v "$1" >/dev/null 2>&1
}

sha256_file() {
  if have shasum; then
    shasum -a 256 "$1" | awk '{print $1}'
  elif have sha256sum; then
    sha256sum "$1" | awk '{print $1}'
  else
    echo "Neither shasum nor sha256sum is available." >&2
    exit 1
  fi
}

os="$(uname -s)"
jobs=1
cc_compiler=""
cxx_compiler=""
cppflags=""
ldflags=""
libs="-lmpfr -lgmpxx -lgmp"
configure_platform_args=()

case "$os" in
  Darwin)
    if ! xcode-select -p >/dev/null 2>&1; then
      cat >&2 <<'EOF'
Apple's command-line tools are required.  Run

  xcode-select --install

finish the installation window, and then run this script again.
EOF
      exit 1
    fi
    if ! have brew; then
      cat >&2 <<'EOF'
Homebrew is required for the macOS installation.  Install it using the command
shown at https://brew.sh/ and then run this script again.  At the time this
script was written, the Homebrew command was

  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
EOF
      exit 1
    fi
    if [[ "${AYMH_SKIP_SYSTEM_PACKAGES:-0}" != "1" ]]; then
      say "Installing the required Homebrew packages"
      brew install gcc gmp mpfr pkg-config
    fi
    gcc_prefix="$(brew --prefix gcc)"
    mpfr_prefix="$(brew --prefix mpfr)"
    gmp_prefix="$(brew --prefix gmp)"
    cc_compiler="$(find "$gcc_prefix/bin" -maxdepth 1 -type f -name 'gcc-[0-9]*' | sort | tail -1)"
    cxx_compiler="$(find "$gcc_prefix/bin" -maxdepth 1 -type f -name 'g++-[0-9]*' | sort | tail -1)"
    if [[ -z "$cc_compiler" || -z "$cxx_compiler" ]]; then
      echo "The Homebrew GCC executables could not be located." >&2
      exit 1
    fi
    cppflags="-I$mpfr_prefix/include -I$gmp_prefix/include"
    ldflags="-L$mpfr_prefix/lib -L$gmp_prefix/lib"
    machine="$(uname -m)"
    if [[ "$machine" == "arm64" ]]; then
      machine="aarch64"
    fi
    configure_platform_args=("--build=${machine}-apple-darwin" "--host=${machine}-apple-darwin")
    jobs="$(sysctl -n hw.ncpu 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')"
    ;;
  Linux)
    if [[ "${AYMH_SKIP_SYSTEM_PACKAGES:-0}" != "1" ]]; then
      if have apt-get; then
        say "Installing the required Debian/Ubuntu packages"
        sudo apt-get update
        sudo apt-get install -y build-essential curl git tar pkg-config libgmp-dev libmpfr-dev
      else
        cat >&2 <<'EOF'
This installer supports Debian and Ubuntu directly.  On another Linux
distribution, install a C++ compiler, make, curl, tar, pkg-config, GMP, and
MPFR, and rerun with AYMH_SKIP_SYSTEM_PACKAGES=1.
EOF
        exit 1
      fi
    fi
    cc_compiler="${CC:-gcc}"
    cxx_compiler="${CXX:-g++}"
    if have nproc; then
      jobs="$(nproc)"
    fi
    ;;
  *)
    cat >&2 <<EOF
Unsupported operating system: $os

On Windows, install Ubuntu under WSL2 and run this script from an Ubuntu
terminal.  Native Windows CAPD builds are not supported by this repository.
EOF
    exit 1
    ;;
esac

requested_jobs="${CAPD_BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-}}"
if [[ -n "$requested_jobs" ]]; then
  if [[ ! "$requested_jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "CAPD_BUILD_JOBS or SLURM_CPUS_PER_TASK must be a positive integer." >&2
    exit 1
  fi
  jobs="$requested_jobs"
fi

for program in curl tar make pkg-config "$cc_compiler" "$cxx_compiler"; do
  if ! have "$program"; then
    echo "Required program not found: $program" >&2
    exit 1
  fi
done

mkdir -p "$DOWNLOAD_DIR" "$LOG_DIR"
timestamp="$(date +%Y%m%d_%H%M%S)"

if [[ ! -f "$ARCHIVE_PATH" ]]; then
  say "Downloading the CAPD ${CAPD_VERSION} source archive"
  curl -fL --retry 3 -o "$ARCHIVE_PATH" "$CAPD_URL"
fi

actual_sha256="$(sha256_file "$ARCHIVE_PATH")"
if [[ "$actual_sha256" != "$CAPD_SHA256" ]]; then
  cat >&2 <<EOF
The CAPD archive has the wrong SHA-256 hash.

Expected: $CAPD_SHA256
Found:    $actual_sha256

Remove $ARCHIVE_PATH and run the script again.  Do not build an archive that
fails this check.
EOF
  exit 1
fi
say "CAPD source archive: SHA-256 verified"

if [[ ! -x "$SOURCE_DIR/configure" ]]; then
  if [[ -e "$SOURCE_DIR" ]]; then
    incomplete_source="${SOURCE_DIR}.incomplete.${timestamp}"
    say "An incomplete CAPD source directory was found"
    printf 'Moving it to %s\n' "$incomplete_source"
    mv "$SOURCE_DIR" "$incomplete_source"
  fi
  say "Unpacking CAPD ${CAPD_VERSION}"
  tar -xzf "$ARCHIVE_PATH" -C "$CAPD_HOME"
fi

build_log="$LOG_DIR/capd-${CAPD_VERSION}-build-${timestamp}.log"
{
  printf 'CAPD version: %s\n' "$CAPD_VERSION"
  printf 'source directory: %s\n' "$SOURCE_DIR"
  printf 'build directory: %s\n' "$BUILD_DIR"
  printf 'C compiler: %s\n' "$cc_compiler"
  printf 'C++ compiler: %s\n' "$cxx_compiler"
} > "$build_log"

# CAPD 5.3.0 can generate this helper without its executable bit on Linux.
# Normalize the permission both for an existing configuration and for a new
# one below.
if [[ -f "$BUILD_DIR/bin/mpcapd-config" ]]; then
  chmod u+x "$BUILD_DIR/bin/mpcapd-config"
fi

configuration_ready=0
if [[ -f "$CONFIG_STAMP" && -x "$BUILD_DIR/config.status" &&
      -x "$BUILD_DIR/bin/mpcapd-config" ]]; then
  configuration_ready=1
elif [[ ! -f "$CONFIG_STAMP" && -x "$BUILD_DIR/config.status" &&
        -x "$BUILD_DIR/bin/mpcapd-config" ]]; then
  # Accept a complete build directory made by an earlier version of this
  # installer, and record that its configure step is complete.
  touch "$CONFIG_STAMP"
  configuration_ready=1
fi

if [[ "$configuration_ready" -eq 0 ]]; then
  if [[ -d "$BUILD_DIR" &&
        -n "$(find "$BUILD_DIR" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
    incomplete_dir="${BUILD_DIR}.incomplete.${timestamp}"
    say "An incomplete CAPD build directory was found"
    printf 'Moving it to %s\n' "$incomplete_dir"
    mv "$BUILD_DIR" "$incomplete_dir"
  fi
  mkdir -p "$BUILD_DIR"
  say "Configuring CAPD with MPFR multiprecision support"
  if ! (
    cd "$BUILD_DIR"
    # CAPD 5.3.0's recursive configure files require this relative source
    # path.  An absolute path leads to malformed paths in the subdirectories.
    source_configure="../capd-${CAPD_VERSION}/configure"
    env \
      CC="$cc_compiler" \
      CXX="$cxx_compiler" \
      CPPFLAGS="$cppflags" \
      LDFLAGS="$ldflags" \
      LIBS="$libs" \
      "$source_configure" \
        --disable-option-checking \
        "${configure_platform_args[@]}" \
        --with-mpfr \
        --with-filib=no \
        --prefix="$CAPD_PREFIX"
  ) >> "$build_log" 2>&1; then
    echo "CAPD configuration failed.  The last 40 lines are:" >&2
    tail -n 40 "$build_log" >&2
    echo "Complete log: $build_log" >&2
    exit 1
  fi
  if [[ -f "$BUILD_DIR/bin/mpcapd-config" ]]; then
    chmod u+x "$BUILD_DIR/bin/mpcapd-config"
  fi
  if [[ ! -x "$BUILD_DIR/config.status" ||
        ! -x "$BUILD_DIR/bin/mpcapd-config" ]]; then
    echo "CAPD configuration did not create the expected files." >&2
    echo "Complete log: $build_log" >&2
    exit 1
  fi
  touch "$CONFIG_STAMP"
else
  say "Reusing the existing CAPD build directory"
  printf 'Reusing an existing configured build directory.\n' >> "$build_log"
fi

say "Building CAPD; this may take several minutes"
if ! make -C "$BUILD_DIR" -j"$jobs" lib >> "$build_log" 2>&1; then
  echo "The CAPD build failed.  The last 40 lines are:" >&2
  tail -n 40 "$build_log" >&2
  echo "Complete log: $build_log" >&2
  exit 1
fi

mpcapd_config="$BUILD_DIR/bin/mpcapd-config"
if [[ ! -x "$mpcapd_config" ]]; then
  echo "The CAPD build did not produce $mpcapd_config." >&2
  exit 1
fi

version="$($mpcapd_config --modversion)"
cflags="$($mpcapd_config --cflags)"
link_flags="$($mpcapd_config --libs)"
if [[ "$version" != "$CAPD_VERSION" ]]; then
  echo "The build reports CAPD $version; CAPD $CAPD_VERSION is required." >&2
  exit 1
fi
if [[ -z "$cflags" || -z "$link_flags" ]]; then
  echo "mpcapd-config did not return the required build flags." >&2
  exit 1
fi
if [[ "$link_flags" != *mpfr* || "$link_flags" != *gmp* ]]; then
  echo "mpcapd-config does not report the MPFR/GMP libraries." >&2
  echo "Returned flags: $link_flags" >&2
  exit 1
fi

say "CAPD multiprecision build completed"
printf 'CAPD version: %s\n' "$version"
printf 'mpcapd-config: %s\n' "$mpcapd_config"
printf 'Build log: %s\n' "$build_log"
cat <<EOF

Return to the certificate repository and run

  bash scripts/run_capd_mp.sh

The run script automatically recognizes the default build location.  If
CAPD_HOME was changed, use the following explicit command instead:

  MPCAPD_CONFIG="$mpcapd_config" bash scripts/run_capd_mp.sh
EOF
