# Release Manifest

This file lists the reviewed files in the certificate package.  Their hashes
are recorded in `SHA256SUMS`.  No precompiled executable is distributed; the
run script builds the program from the listed source and Makefile.

## C++ certificate

- `capd_cpp/src/lppssi_cert_v6.cpp`
- `capd_cpp/Makefile`
- `capd_cpp/transcripts/capd_certificate_v6_220bit_macos.txt`
- `capd_cpp/transcripts/capd_certificate_v6_220bit_linux.txt`

## Documentation and scripts

- `README.md`
- `docs/MANUAL.md`
- `CITATION.cff`
- `.gitignore`
- `LICENSE`
- `MANIFEST.md`
- `scripts/install_capd_5_3_0.sh`
- `scripts/run_capd_mp.sh`
- `scripts/check_hashes.sh`
- `scripts/refresh_hashes.sh`
- `SHA256SUMS`

## Paper reference

The paper itself is not included in this repository.  It is available as
[arXiv:YYMM.NNNNN](https://arxiv.org/abs/YYMM.NNNNN).

Timestamped build logs, local transcripts, and the compiled executable are
created during a local run.  They are intentionally excluded from the release
manifest and from `SHA256SUMS`.
