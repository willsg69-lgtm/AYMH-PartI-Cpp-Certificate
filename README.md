# AYMH Part I C++ Numerical Certificate

This repository contains the reproducibility package for the computer-assisted
part of

**Asymptotic stability of the degree-one vortex in the abelian
Yang-Mills-Higgs model: Spectral Theory and Numerics**

by Jonas Lührmann, José M. Palacios, Fabio Pusateri, Wilhelm Schlag, and
Sohrab Shahshahani; see
[arXiv:2608.10608](https://arxiv.org/abs/2608.10608).

The paper studies the spectral theory of the operator obtained by linearizing
the self-dual abelian Yang--Mills--Higgs equations about the degree-one vortex.
Most of the proof is analytic.  After those arguments have reduced the
remaining questions to explicit inequalities for solutions of radial ordinary
differential equations and associated integrals, the C++ program in this
repository verifies those inequalities by interval arithmetic.  Thus the code
is part of the proof, but it is not a substitute for the analytic lemmas in the
paper.

No knowledge of C++ or of CAPD is required to reproduce the calculation.  The
provided scripts install the required version of CAPD, compile the program,
and run the complete certificate.

## What the certificate verifies

The program follows the order of the mathematical argument.  It stops if any
required earlier check fails, so a later conclusion is never printed as
certified when one of its numerical hypotheses has not been proved.

The principal checks are the following.

1. **Vortex profile.**  Endpoint shooting, an interval Newton calculation,
   and the analytic tail estimates from the paper locate the unique vortex
   shooting parameter.  In particular,
   ```text
   |c_* - 0.6032878545816779| < 7e-17.
   ```
   The program also verifies at ρ₀ = 4 the finite inequalities from which the
   paper obtains the global vortex barriers.

2. **Frobenius data at the origin.**  The regular series for the vortex, its
   shooting derivatives, the internal mode, the distorted Fourier basis, and
   the regular threshold solution are bounded at r₀ = 0.1.  The program checks
   finitely many coefficients and then verifies an all-order induction for the
   remaining coefficients.  It does not infer an infinite geometric tail from
   a finite list of observed coefficient ratios.

3. **The first scalar channel.**  A validated enclosure of the potential gives
   $V_1(r)-1>0.0548$ on the compact interval required by the paper.  Together
   with the analytic argument there, this excludes bound states and a
   threshold resonance in this channel.

4. **The second scalar channel.**  The logarithmic-kernel integral in the
   Setô eigenvalue bound satisfies
   ```text
   Lambda < 3.182292 < 3.2.
   ```
   The paper therefore obtains at most one bound state in this channel.  A
   separate Wronskian calculation gives
   ```text
   W in [-0.8003267, -0.8001801],
   ```
   which excludes a threshold resonance.

5. **Internal mode.**  A Wronskian sign change, combined with the preceding
   eigenvalue count, proves existence and uniqueness of the internal
   eigenvalue.  Every later calculation is performed uniformly on the
   certified interval
   ```text
   lambda^2 in [0.777471875, 0.77747375].
   ```

6. **Large-radius estimates.**  The program verifies the finite inequalities
   used in the comparison with $K_0$, bounds the derivative of the internal
   mode, and constructs the outgoing solutions at $R=16$.  The error in the
   outgoing data is controlled by the Volterra estimates proved in the paper.

7. **Fermi Golden Rule.**  The radial integrals on
   [r₀,R] = [0.1,16] are evaluated by
   validated quadrature.  The omitted origin and large-radius pieces are
   bounded separately and combined with the compact integral before
   squaring.  The final certified lower bounds are
   ```text
   hat D_1^(0)  > 0.017297889,
   hat D_2^(0)  > 0.047248801,
   hat D_12^(0) > 0.022940050.
   ```
   With the sign convention and normalization of the paper, these inequalities
   prove that the physical coefficients $D_1,D_2,D_{12}$ are strictly negative.

The [manual](docs/MANUAL.md#paper-code-interface) identifies the corresponding
paper results and explains precisely which part is proved analytically and
which finite statement is checked by the program.

## Why the numerical statements are rigorous

The program uses CAPD 5.3.0 with its multiprecision interval types
`MpInterval` and `MpFloat`.  It requests 220 bits of precision and uses a
validated Taylor method of order 40.  The relevant features are:

- every decimal endpoint is converted with directed rounding;
- arithmetic operations produce intervals containing the exact result;
- CAPD encloses each ODE solution throughout every integration step, not only
  at the mesh points;
- the regular singular point at the origin is handled by proved Frobenius
  remainder bounds;
- the unbounded intervals are handled by the analytic comparison and Volterra
  arguments in the paper;
- every numerical target is tested by the program, and the complete run ends
  in failure unless all required checks pass.

The two versioned transcripts were produced independently on macOS and
x86-64 Linux.  Apart from the line identifying the compiler, their numerical
contents agree.

## Quick start

The supported systems are:

- macOS with Apple's command-line tools and Homebrew;
- Debian or Ubuntu Linux;
- Windows through Ubuntu under WSL2;
- managed Linux systems on which the required compiler, GMP, and MPFR are
  already available.

An internet connection and at least 2 GB of free disk space are recommended.
On macOS, install Apple's command-line tools and Homebrew first; on Windows,
install WSL2 with Ubuntu first.  The manual gives the exact preparation
commands.  CAPD is a C++ library; there is no separate CAPD application or
interactive CAPD shell to operate.

Open a terminal and obtain the repository:

```bash
git clone https://github.com/willsg69-lgtm/AYMH-PartI-Cpp-Certificate.git
cd AYMH-PartI-Cpp-Certificate
```

If `git` is not installed, use GitHub's **Code → Download ZIP** command,
unpack the archive, open a terminal, and change to the unpacked directory.
The commands below invoke each script through Bash and therefore do not depend
on executable permissions being preserved by the ZIP archive.

Before installing or compiling anything, check the downloaded files:

```bash
bash scripts/check_hashes.sh
```

Every listed file should be reported as `OK`.  Next install CAPD 5.3.0 with
MPFR/GMP multiprecision support:

```bash
bash scripts/install_capd_5_3_0.sh
```

On macOS the installer uses Homebrew.  On Debian, Ubuntu, and WSL2 Ubuntu it
uses `apt` and may request the user's password while installing system
packages.  CAPD itself is built under
`$HOME/.local/share/aymh-capd`; no administrative privilege is needed for that
part of the installation.  The downloaded CAPD archive is accepted only if
its SHA-256 hash agrees with the fixed value recorded in the installer.

Finally, compile and run the complete certificate:

```bash
bash scripts/run_capd_mp.sh
```

The script makes a clean build and automatically locates the CAPD installation
made by the preceding command.  At the beginning of a correct run one should
see

```text
interval backend = MpInterval/MpFloat, requested precision bits = 220
Taylor order = 40
```

The calculation has succeeded only if it finishes with

```text
OVERALL: PASS FOR THE NUMERICAL CERTIFICATE BLOCKS
```

The script saves the compiler output and the complete numerical output in
timestamped files under `capd_cpp/` and prints their names at the end.  A failed
or interrupted run proves nothing; consult the
[troubleshooting section of the manual](docs/MANUAL.md#if-a-command-fails)
before trying again.

Detailed preparation instructions for macOS, Linux, managed systems, and WSL2
are given in the
[installation section of the manual](docs/MANUAL.md#installation-and-first-run).
If CAPD is already installed elsewhere, use

```bash
MPCAPD_CONFIG=/absolute/path/to/mpcapd-config bash scripts/run_capd_mp.sh
```

## Reference transcripts and file verification

The repository contains two complete passing transcripts:

- `capd_cpp/transcripts/capd_certificate_v6_220bit_macos.txt`;
- `capd_cpp/transcripts/capd_certificate_v6_220bit_linux.txt`.

The local transcript produced by `run_capd_mp.sh` can be compared with the
transcript for the corresponding platform.  Compiler-identification lines may
differ when a different compiler release is used.  The two included
transcripts have identical mathematical contents.  Another compatible
compiler may produce slightly different outward-rounded endpoints, but every
required strict inequality and every `PASS` statement must remain valid.  The
manual gives a precise comparison procedure.

`SHA256SUMS` covers every versioned file in the package.  It can be checked at
any time with

```bash
bash scripts/check_hashes.sh
```

Changing a file necessarily changes its hash.  Maintainers who intentionally
modify the repository regenerate the list with

```bash
bash scripts/refresh_hashes.sh
bash scripts/check_hashes.sh
```

A documentation-only change does not require a new numerical run.  A change
to the C++ source or a numerical input requires a complete clean run and a new
reviewed reference transcript.

## Repository contents

- [`capd_cpp/src/lppssi_cert_v6.cpp`](capd_cpp/src/lppssi_cert_v6.cpp) is the
  complete certificate program.
- [`capd_cpp/Makefile`](capd_cpp/Makefile) performs the reviewed
  multiprecision build.
- [`capd_cpp/transcripts/`](capd_cpp/transcripts/) contains the macOS and Linux
  reference transcripts.
- [`scripts/install_capd_5_3_0.sh`](scripts/install_capd_5_3_0.sh) installs the
  fixed CAPD release.
- [`scripts/run_capd_mp.sh`](scripts/run_capd_mp.sh) clean-builds and runs the
  certificate.
- [`scripts/check_hashes.sh`](scripts/check_hashes.sh) verifies the package.
- [`docs/MANUAL.md`](docs/MANUAL.md) explains installation, troubleshooting,
  and the paper--code correspondence in detail.
- [`CITATION.cff`](CITATION.cff) gives the citation metadata.
- [`SHA256SUMS`](SHA256SUMS) records the versioned-file hashes.

## FGR normalization

For completeness, the code uses the internal eigenfunction normalized by
$\psi_0(0)=1$ and certifies positive auxiliary quantities
$\widehat D_\alpha^{(0)}$.  In the notation of the paper,

$$
\begin{aligned}
\widehat D_\alpha^{(0)}
  &=-\frac{k_\lambda}{4\lambda}D_\alpha^{(0)},\\
D_\alpha
  &=-\frac{4\lambda}{k_\lambda\nu_\psi^4}\widehat D_\alpha^{(0)},\\
k_\lambda&=\sqrt{4\lambda^2-1}.
\end{aligned}
$$

All factors displayed apart from the minus signs are positive.  This is why
the positive lower bounds printed by the program imply the negativity of the
physical FGR coefficients.  The derivation of these formulas and the
distorted-Fourier normalization are explained in the paper and in the
[paper--code interface](docs/MANUAL.md#fgr-source-bound-and-final-lower-bounds).

## Citation and license

If this package is used in other work, please cite the accompanying paper and
the repository using the metadata in [`CITATION.cff`](CITATION.cff).

The original code and documentation in this repository are released under the
[MIT License](LICENSE).  CAPD is an external dependency, is not distributed in
this repository, and is governed by its own license.
