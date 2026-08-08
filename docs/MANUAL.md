# Manual for the C++ Numerical Certificate

This manual supplements the [README](../README.md).  Readers who only want to
reproduce the certificate can follow the README's quick start.  The sections
below provide the platform-specific details, explain the correspondence
between the paper and the program, and describe how to check a new transcript.

- [Scope](#scope)
- [Numerical environment](#numerical-environment)
- [Installation and first run](#installation-and-first-run)
- [Paper--code interface](#paper-code-interface)
- [Checking transcripts and file hashes](#checking-transcripts-and-file-hashes)
- [Repository scope and responsibility](#repository-scope-and-responsibility)

## Scope

This repository contains the C++ CAPD interval-arithmetic code supporting the
computer-assisted estimates in *Asymptotic stability of the degree-one vortex
in the abelian Yang-Mills-Higgs model: Spectral Theory and Numerics*; see
[arXiv:YYMM.NNNNN](https://arxiv.org/abs/YYMM.NNNNN).  The code certifies
finite numerical statements after the analytic reductions in the paper have
been made.  It does not prove the analytic shooting dichotomy, the vortex tail
lemmas, Setô's eigenvalue bound, or the distorted Fourier theory.  Those
mathematical inputs, and the finite statements checked by the program, are
identified in the [paper--code interface](#paper-code-interface).

The relevant source file is

```text
capd_cpp/src/lppssi_cert_v6.cpp
```

and the reference passing transcript is

```text
capd_cpp/transcripts/capd_certificate_v6_220bit_macos.txt
```

An independent Linux reproduction is recorded in

```text
capd_cpp/transcripts/capd_certificate_v6_220bit_linux.txt
```

The file label `v6` is only a source-version label.  The mathematical claims
are identified by the lemmas and intervals stated in the paper.

## Numerical Environment

The run uses CAPD's multiprecision interval types `MpInterval` and `MpFloat`,
with 220 requested bits.  Every decimal input is read by CAPD with downward
rounding at the left endpoint and upward rounding at the right endpoint.  The
program begins with a direct check of this convention for representative
non-dyadic inputs.  The transcript records the source SHA-256, CAPD version,
compiler version, interval-arithmetic types, and requested precision.  The
calculation uses outward interval arithmetic and validated Taylor propagation
rather than ordinary floating-point approximations.

Throughout this manual, an interval is said to contain or enclose a quantity
when directed rounding proves that the exact quantity belongs to that
interval.  For an ODE integration, CAPD produces interval boxes containing the
solution throughout each time or radial step.  A transcript is the complete
plain-text record written by one run of the certificate.

The build requires:

- a C++17 compiler;
- CAPD 5.3.0 built with MPFR/GMP multiprecision support;
- an `mpcapd-config` executable or script returning CAPD compile and link flags.

## Installation and first run

CAPD is a C++ library rather than an interactive application.  A user does not
have to write C++ or enter commands in a CAPD shell.  The repository provides
one script that installs the required CAPD version and another that compiles
and runs the certificate.

The reference transcript was produced with CAPD 5.3.0 and Homebrew GCC 15.2.0.
The installation script downloads the fixed
[CAPD 5.3.0 source archive](https://sourceforge.net/projects/capd/files/5.3.0/)
and checks that its SHA-256 hash is

```text
e4100959a5409d330f8907d050f101a0485489075b4ce0d5eb2e349a2f8bf228.
```

The current CAPD project and its general build instructions are available from
the [official CAPD repository](https://github.com/CAPDGroup/CAPD).  This
certificate deliberately uses the stable 5.3.0 archive.  That archive uses its
provided `configure` script; it should not be built with instructions intended
only for a newer development checkout.

An internet connection is needed to obtain prerequisite packages and the CAPD
source archive.  Before installing anything, verify the files in the downloaded
repository from its root directory:

```bash
bash scripts/check_hashes.sh
```

The tested CAPD source and build tree occupies about 200 MB.  The compiler and
the prerequisite packages require additional space, so at least 2 GB of free
space is prudent.  Building CAPD normally takes several minutes.  The complete
configuration and compiler output is saved under
`$HOME/.local/share/aymh-capd/logs`; the terminal displays the location of this
log when the build finishes or fails.

### macOS

Open the Terminal application.  Install Apple's command-line tools once by
running

```bash
xcode-select --install
```

and complete the installation window that appears.  The automated installer
also requires Homebrew.  If `brew --version` says that the command is not
found, install Homebrew using the command currently given at
[brew.sh](https://brew.sh/), then close and reopen Terminal.

Change to the downloaded or cloned certificate repository.  For example,

```bash
cd "$HOME/Downloads/AYMH-PartI-Cpp-Certificate"
```

with the path changed if the repository is elsewhere.  Then run

```bash
bash scripts/install_capd_5_3_0.sh
```

The script uses Homebrew to install GCC, GMP, MPFR, and `pkg-config`.  It then
downloads, verifies, and builds CAPD beneath `$HOME/.local/share`.  The CAPD
build itself does not require administrative privileges.

### Debian or Ubuntu Linux

Open a terminal, change to the certificate repository, and run

```bash
bash scripts/install_capd_5_3_0.sh
```

The script uses `sudo apt-get` to install the compiler, `make`, `curl`, `tar`,
`pkg-config`, GMP, and MPFR.  It then builds CAPD beneath `$HOME/.local/share`.
The system may ask for the user's password while installing the prerequisite
packages; no password is needed for the subsequent CAPD build.

### Managed Linux systems without administrative access

On a managed Linux system without `sudo`, provide a C++17 compiler, GMP,
MPFR, `make`, `curl`, `tar`, and `pkg-config` through the system's software
environment.  Then build CAPD in the user's home directory without requesting
system-package installation:

```bash
AYMH_SKIP_SYSTEM_PACKAGES=1 \
CAPD_HOME="$HOME/.local/share/aymh-capd" \
bash scripts/install_capd_5_3_0.sh
```

If the system assigns a fixed number of build cores, set
`CAPD_BUILD_JOBS` to that number before running the installer.  After CAPD has
been built, run the certificate with

```bash
bash scripts/run_capd_mp.sh
```

The complete certificate has been independently reproduced on x86-64 Linux
with CAPD 5.3.0 and GCC 12.2.0.  Apart from the compiler-identification line,
the Linux transcript included in this repository is identical to the macOS
reference transcript.

### Windows

The supported Windows route is Ubuntu under Windows Subsystem for Linux
(WSL2).  Native Windows compilation has not been used to produce or check the
reference certificate.  The current prerequisites are given in Microsoft's
[official WSL installation instructions](https://learn.microsoft.com/windows/wsl/install).
In an Administrator PowerShell window, run

```powershell
wsl --install -d Ubuntu
```

Restart Windows if requested, open the newly installed Ubuntu application, and
complete its request for a Linux username and password.  Place the repository
inside the Linux home directory rather than under `/mnt/c`; this avoids slow
compilation on the Windows-mounted filesystem.  One way to do this is

```bash
sudo apt update
sudo apt install -y git
cd "$HOME"
git clone https://github.com/willsg69-lgtm/AYMH-PartI-Cpp-Certificate.git
cd AYMH-PartI-Cpp-Certificate
```

Then follow the Linux route:

```bash
bash scripts/check_hashes.sh
bash scripts/install_capd_5_3_0.sh
bash scripts/run_capd_mp.sh
```

### Compiling and running the certificate

After the CAPD installer finishes, remain in or return to the repository root
and run

```bash
bash scripts/run_capd_mp.sh
```

The run script automatically finds the default CAPD build made by
`install_capd_5_3_0.sh`.  It removes any previous certificate executable,
compiles a fresh one, and runs it.  No editing of the Makefile or C++ source is
needed.  The compiler output is saved in an ignored, timestamped file named
`capd_cpp/capd_build_v6_mp_local_*.txt`.  The numerical output remains visible
in the terminal and is also saved in
`capd_cpp/capd_certificate_v6_220bit_local_*.txt`.

Internally, the run script invokes `make mp`, which defines
`LPPSSI_USE_MP`.  This selects CAPD's `MpFloat`, `MpInterval`, `MpIVector`, and
the corresponding validated ODE solver.  At the beginning of `main`, before
the interval data are constructed, the executable calls

```text
MpFloat::setDefaultPrecision(220).
```

Thus the script compiles the multiprecision version, while the executable sets
the precision to 220 bits at runtime.  The value 220 is not a compiler option.

Before compiling, the script requires the exact CAPD version 5.3.0 and checks
that the reported link flags contain MPFR and GMP.  It then checks automatically
that the transcript contains

```text
CAPD version = 5.3.0
interval backend = MpInterval/MpFloat, requested precision bits = 220
Taylor order = 40
```

and ends with

```text
OVERALL: PASS FOR THE NUMERICAL CERTIFICATE BLOCKS
```

If any required line is absent, the script exits with an error.  In particular,
output mentioning `DInterval/double` or 53 bits is not a valid proof run.

An advanced user who already has a suitable multiprecision CAPD installation
in a different location may bypass the installer and run

```bash
MPCAPD_CONFIG=/absolute/path/to/mpcapd-config bash scripts/run_capd_mp.sh
```

The `mpcapd-config` program must return nonempty results for `--cflags` and
`--libs`, must report version 5.3.0, and must link the MPFR/GMP version of CAPD.
The two transcripts in `capd_cpp/transcripts` are the reviewed passing runs.
A new local transcript should not replace either of them without review.

### If a command fails

The installation script can be run again safely after an interruption.  It
reuses a completed configuration; if it finds an incomplete source or build
directory, it moves that directory aside with an `.incomplete` suffix and
starts the affected step again.  It does not silently delete the earlier
directory.

If the shell reports `Permission denied`, invoke the scripts through Bash:

```bash
bash scripts/install_capd_5_3_0.sh
bash scripts/run_capd_mp.sh
```

If configuration or compilation fails, the script prints the last 40 lines of
the relevant log and gives the full path to that log.  Those lines normally
identify the missing package or compiler error.  Ordinary compiler warnings
are retained in the log rather than filling the terminal.  An error saying
that the CAPD version is not 5.3.0 usually means that a different
`mpcapd-config` was selected; use the explicit `MPCAPD_CONFIG=...` command above
to select the build made by the installer.

The time required depends on the processor and compiler.  Leave the terminal
open until the final `OVERALL: PASS` line and the local transcript hash have
appeared.  Closing the terminal or interrupting the program before that point
does not produce a certificate.

## Paper-Code Interface

The paper reduces the required numerical input to a finite list of interval
statements.  The C++ code certifies those statements in the following order.

### Vortex shooting bracket

The broad interval

```text
J_0 = [0.6032878545810, 0.6032878545819]
```

is checked by endpoint shooting.  At the lower endpoint the code proves
`a(27.75)>2`; at the upper endpoint it proves `U(31.6)>2`.  Together with the
shooting dichotomy in Proposition `prop:exuniqUa`, this gives the enclosure in
Proposition `prop:J0slope`.

### Barrier check at r = 4

Using the broad interval `J_0`, the code checks the four finite-radius
inequalities at `r=4` required by Lemma `lem_1-U2`.  Once these inequalities
are certified, that lemma and Corollary `cor:vortex_tail_bridge` provide the
global bounds for `r>=4`.

### Newton refinement of the shooting parameter

The code performs the interval Newton step in Lemma `lem:Newton` for

```text
F_20(c) = U(20;c)-1.
```

The transcript records the Newton image, the derivative lower bound for
`partial_c U(20;c)` on `J_0`, the analytic bound for `1-U(20;c_*)`, and the
resulting interval containing the true asymptotic shooting parameter.  The
working interval used later is

```text
I_cert = [0.6032878545816699, 0.6032878545816856].
```

### Frobenius initial data at r = 0.1

The code certifies the coefficient and series-tail estimates used to initialize:

- the vortex profile;
- the shooting derivatives `partial_c U` and `partial_c a`;
- the internal mode `psi`;
- the solutions `Phi_1` and `Phi_2` determined by their behavior at the origin;
- the threshold regular solution `Phi_{2,0}^{(0)}`.

For `psi` the calculation retains coefficients through degree `M=100`; for
each `Phi`-series it retains coefficients through degree `M=60`.  It computes
and checks the next `K=20` coefficients with `rho=0.95`.  A separate all-order
induction estimate then proves that the full coefficient recurrence preserves
the proposed geometric bound for every remaining index.  Only after this
all-order estimate has passed are the resulting value and derivative tail
bounds used to initialize an ODE calculation.  The terminal coefficient ratios
printed in the transcript are finite consistency checks; they are not
extrapolated to justify an infinite geometric series.

### Positivity of H_1

On the compact interval needed by the scalar operator `H_1`, the code evaluates
the potential and proves

```text
V_1(r)-1 > 0.0548.
```

This is the finite interval input used in Lemma `lem_appV1g1` to exclude bound
states and a threshold resonance for `H_1`.

### Setô eigenvalue count

The code evaluates the logarithmic-kernel integrals entering the estimate in
Proposition `theoLT2` and certifies

```text
Lambda < 3.182292 < 3.2.
```

Together with the counting theorem quoted in the paper, this gives at most one
bound state for `L_2`.

### Threshold nonresonance

The code starts the regular threshold solution at `r=0.1`.  A Volterra estimate
provides an interval containing the solution normalized at infinity and its
derivative at the finite starting radius.  CAPD then carries both solutions to
the matching radius, where the program computes their Wronskian.  The
transcript gives an interval contained in

```text
[-0.8003267, -0.8001801].
```

Since this interval excludes zero, Lemma `Lemspec0` rules out a threshold
resonance.

### Internal eigenvalue and Lambda_FGR

The broad Wronskian sign check in Lemma `lem:fgr_eigenvalue_box` gives an
eigenvalue in

```text
[0.77747, 0.77753].
```

Using Proposition `theoLT2`, based on Setô's bound, this eigenvalue is unique.
The code then evaluates the Wronskian at directed intervals containing the
exact endpoints of the sharper interval printed in the paper:

```text
Lambda_FGR = [0.777471875, 0.77747375].
```

The endpoint Wronskians have opposite strict signs.  An outward-rounded
interval containing these exact decimal endpoints is then used for all later
calculations.  Consequently, all subsequent spectral estimates hold uniformly
for every `mu` in `Lambda_FGR` and hence at `mu=lambda^2`.

### K_0 comparison and derivative comparison for the internal mode

For `delta_0=0.0005` and `kappa=sqrt(1-delta_0-mu)`, the code checks the
constants needed for

```text
psi(r)       <= 1.4 K_0(kappa r),
|psi'(r)|   <= 1.4 K_0(kappa r),
```

in the ranges stated in Lemmas `lem:psi_K0` and
`lem:capd_psi_prime_tail`.

### Outgoing data at R = 16

The code constructs outgoing scalar data from the truncated asymptotic
expansion at `R=16`, computes the exact error left by that truncation, and
verifies the constants in a Volterra remainder estimate.  The analytic vortex
tail controls the difference between the exact scalar equation and the limiting
outgoing equation.

The transcript records the matrix Weyl upper bound

```text
W_infty <= 1.253.
```

This is the bound proved in Lemma `lem:capd_outgoing_start` and used in the
FGR estimates.

### FGR source bound and final lower bounds

The code evaluates the FGR integrals on `[0.1,16]` by interval quadrature.  The
interval `(0,0.1)` is treated in Lemma `lem:capd_fgr_origin`, and
`[16,infinity)` is treated in Lemma `lem:capd_fgr_tail_upper_bound` using the
Bessel comparison, derivative comparison, and outgoing-solution bound.
For the compact calculation the program stores the complex conjugates of the
two-vector integrals denoted by `A_ij` in the paper; their Euclidean norms are
identical. With the corrected distorted-Fourier normalization, the common
scalar factor outside each radial vector integral on the shell is
`1/(2 lambda sqrt(2 pi))`. The direct estimate in the paper gives a constant
smaller than `8.02`, and hence the convenient bound `30` remains valid. The
exact prefactor of the tail contribution to `hat D_alpha^(0)` is
`pi k_lambda/(8 lambda^2)`. The certificate estimates the larger expression
with prefactor `pi k_lambda/(4 lambda^2)`, so its printed tail bound remains
rigorous and is slightly conservative.

The code uses the `k`-shell factor described in Remark
`rem:fgr_code_normalization` of the paper:

```text
k_mu = sqrt(4 mu - 1),
factor = k_mu/(64 mu^2).
```

After applying the triangle inequality to the compact, origin, and
large-radius vector integrals, and only then squaring, the transcript gives

```text
hat D_1^(0)  > 0.017297889
hat D_2^(0)  > 0.047248801
hat D_12^(0) > 0.022940050
```

These are the quantities formed with `psi_0(0)=1` and the factor
`k_mu/(64 mu^2)`.  In the notation of Proposition `prop:FGR1`,

```text
hat D_alpha^(0) = -(k_lambda/(4 lambda)) D_alpha^(0),
D_alpha = -(4 lambda/(k_lambda nu_psi^4)) hat D_alpha^(0).
```

The conversion factor is positive, so these three bounds prove the negativity
of the physical coefficients `D_1,D_2,D_12`.

## Checking transcripts and file hashes

### Comparing a new run with a reference transcript

The run script prints the path of the new transcript when it finishes.  From
the repository root, one may also select the most recent local transcript with

```bash
local_transcript="$(ls -t capd_cpp/capd_certificate_v6_220bit_local_*.txt | head -n 1)"
```

On macOS set

```bash
reference_transcript="capd_cpp/transcripts/capd_certificate_v6_220bit_macos.txt"
```

and on Linux or WSL2 set

```bash
reference_transcript="capd_cpp/transcripts/capd_certificate_v6_220bit_linux.txt"
```

First inspect the three identifying lines:

```bash
head -n 3 "$local_transcript"
```

They record the source hash, CAPD version, and compiler.  The compiler line may
differ when another compatible compiler release is used.  The remaining lines
can be compared with

```bash
diff -u \
  <(tail -n +4 "$reference_transcript") \
  <(tail -n +4 "$local_transcript")
```

The two versioned transcripts have identical mathematical contents in this
comparison.  A new run with a different compiler may produce slightly
different outward-rounded endpoints.  Such a difference should be examined,
not hidden: every strict inequality must still hold, every required block must
say `PASS`, and the run must end with the overall `PASS` statement.

### Checking the repository files

The file `SHA256SUMS` records SHA-256 hashes for the repository contents. Run

```bash
bash scripts/check_hashes.sh
```

from the repository root.  Each listed file should be reported as `OK`.  A
failure means that the local file differs from the reviewed package; do not
silently regenerate the hashes unless the change was intentional and has been
reviewed.

### For maintainers and contributors

If the C++ source or a numerical input is changed, first make
and review a complete clean run and replace the reference transcript.
Documentation-only changes do not require a numerical rerun.  After any
intentional change to a versioned file, regenerate `SHA256SUMS` with

```bash
bash scripts/refresh_hashes.sh
bash scripts/check_hashes.sh
```

and inspect the resulting diff before committing.

## Repository scope and responsibility

Only the C++ CAPD certificate is part of this repository.  Earlier exploratory
computations and independent experimental checks are not included and are not
used as proof inputs.  The paper source is also not included; the paper is
available as [arXiv:YYMM.NNNNN](https://arxiv.org/abs/YYMM.NNNNN).  The
repository is limited to the source, transcripts, hashes, and documentation
needed to reproduce the C++ interval certificate.

Codex assisted the authors in preparing and checking the code package and its
documentation.  The mathematical arguments, numerical claims, and final
verification decisions remain the responsibility of the authors.
