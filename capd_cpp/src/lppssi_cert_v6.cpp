// LPPSSI C++ CAPD numerical certificate, v6
//
// This file is intentionally self-contained as a numerical certificate driver.
// CAPD supplies rigorous finite-time ODE enclosures; every global input that is
// not a finite-time ODE step is either produced by a named code block below or
// explicitly identified as an analytic lemma imported from the paper.  Thus this
// program certifies the finite numerical hypotheses of the paper's analytic
// lemmas; it is not a proof of those analytic lemmas themselves.
//
// Connection with the Part I paper:
//   * vortex shooting / J0 corresponds to the paper's vortex-profile interval;
//   * H1 positivity corresponds to the compact numerical part of the V1-1
//     positivity lemma;
//   * the Seto logarithmic-kernel block supplies the certified one-bound-state
//     input used for uniqueness of the internal mode;
//   * the threshold block supplies the non-resonance Wronskian certificate;
//   * the K0 comparison and FGR blocks certify the numerical hypotheses in the
//     internal-mode tail and Fermi Golden Rule statements.
//
// Version v5 made one structural correction over v4: the eigenvalue interval
// used by K0/FGR is no longer a fixed decimal.  The code first constructs it by
// validated Wronskian bisection inside the already-certified broad
// eigenvalue interval, installs that certified interval into Settings::mu_box,
// and only then runs the spectral-origin, K0, and FGR computations.
//
// Version v6 addresses two structural checks:
//   * the outgoing FGR Weyl coefficients are no longer hard-coded.  They are
//     generated from the limiting equations solved by vfOutgoing and checked by
//     the exact recurrence residual;
//   * the shooting-derivative tail for a_c no longer includes the unused
//     q[nMax+1] storage slot, so the tail certificate cannot accidentally use a
//     zero placeholder as its final coefficient.
//
// The present corrected v6 source also constructs every decimal input through
// CAPD's directed string parser.  Thus the mathematical decimal endpoints used
// in the paper are contained in the corresponding machine intervals.  It also
// checks the exact printed I_cert and Lambda_FGR endpoints explicitly.
//
// Build:
//   CAPD_CONFIG=/path/to/CAPD/build/bin/capd-config make
//
// Notes:
//  * The code supports the ordinary CAPD interval typedefs supplied by
//    capd-config, with guarded multiprecision aliases for installations that
//    provide mpcapd.
//  * The code is fail-closed and fail-fast: a dependent block is not run unless
//    all prerequisite blocks in the same run have passed.
//  * This source is written as literate proof code.  Search for "CODED LEMMA"
//    to see the exact auxiliary estimates being certified.

#include <capd/capdlib.h>
#ifdef LPPSSI_USE_MP
#include <capd/mpcapdlib.h>
#endif

#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace capd;

namespace lppssi {

// =============================================================================
// What is CAPD doing in this file?
// =============================================================================
//
// CAPD is the external rigorous-numerics library used here.  The mathematical
// inequalities, tail estimates, Frobenius recurrences, and FGR reductions below
// are our "coded lemmas"; CAPD supplies the validated interval arithmetic and
// finite-time ODE propagation used inside those lemmas.
//
// The CAPD objects used in this source are:
//
//   DInterval / MpInterval
//     A real interval scalar.  Operations such as +,-,*,/, sqrt, exp, log, sin,
//     cos are outward rounded by CAPD, so the result encloses every exact value
//     obtained from inputs in the given intervals.  In an ordinary installation
//     we use DInterval, based on directed-rounded doubles.  If the program is
//     compiled with LPPSSI_USE_MP and CAPD's mpcapd headers are available, the
//     aliases switch to MpInterval/MpFloat for multiprecision intervals.
//
//   IVector / MpIVector
//     A vector whose entries are intervals.  Initial data such as
//       r in [0.1,0.1], U in U0, a in a0, c in c_box
//     are stored as one such interval vector.  This vector represents a whole
//     box of initial conditions, not a single floating-point trajectory.
//
//   IMap / MpIMap
//     CAPD's parsed vector-field class.  We give it strings of the form
//       "var:r,u,a,c; fun:1,((1-a)/r)*u,r*(1-u^2)/2,0;"
//     and CAPD builds the interval Taylor model of that right-hand side,
//     including derivatives needed by the validated Taylor solver.
//
//   IOdeSolver / MpIOdeSolver
//     CAPD's interval Taylor ODE solver.  Given an IMap and a Taylor order, it
//     advances interval sets with a rigorous local truncation/remainder
//     enclosure.  This is the finite-time validation engine: if x0 is an
//     interval box, the output encloses all exact solutions starting from x0.
//
//   ITimeMap / MpITimeMap
//     A convenience layer around IOdeSolver which iterates validated time steps
//     until a requested time T is reached.  The call tm(T,set) returns a rigorous
//     endpoint enclosure.  The overload tm(T,set,curve) also stores a validated
//     solution curve; evaluating curve([0,T]) gives a tube enclosing the solution
//     during the whole time interval.  We use that for rectangle/quadrature
//     estimates where pointwise endpoint data are not enough.
//
//   C0TripletonSet / MpC0TripletonSet
//     CAPD's representation of the evolving interval set.  It is more refined
//     than a naive interval vector: the "tripleton" form tracks the main affine
//     part and remainders separately, reducing wrapping while preserving a
//     rigorous enclosure.  We construct it from an interval vector and let
//     ITimeMap update it.
//
// The short wrappers integrate(...) and integrateTube(...) below are the only
// places in this file where CAPD performs finite-time ODE propagation.  All
// other routines either manufacture rigorous initial boxes/tail bounds for CAPD
// or check algebraic interval inequalities from CAPD's output.

#ifdef LPPSSI_USE_MP
using R = MpFloat;
using I = MpInterval;
using Vec = MpIVector;
using Map = MpIMap;
using Solver = MpIOdeSolver;
using TimeMap = MpITimeMap;
using Set = MpC0TripletonSet;
#else
using R = double;
using I = DInterval;
using Vec = IVector;
using Map = IMap;
using Solver = IOdeSolver;
using TimeMap = ITimeMap;
using Set = C0TripletonSet;
#endif

// -----------------------------------------------------------------------------
// Basic interval utilities
// -----------------------------------------------------------------------------

// CAPD's two-string constructor reads the left endpoint downward and the
// right endpoint upward.  In particular, IV("0.1") contains the exact decimal
// 0.1.  Converting the string first to MpFloat (or double) and then making a
// singleton interval would instead retain only a nearest-rounded value.
static I IV(const char* s) { return I(s, s); }
static I IV(const char* a, const char* b) { return I(a, b); }
static I zero() { return IV("0"); }
static I one() { return IV("1"); }
static I two() { return IV("2"); }
static I half() { return IV("0.5"); }
static I piI() { return IV("3.14159265358979323846264338327950288419716939937510", "3.14159265358979323846264338327950288419716939937511"); }

static I sqr(const I& x) { return x * x; }

static I ipow(I x, unsigned n) {
  I y = one();
  while (n) {
    if (n & 1U) y *= x;
    x *= x;
    n >>= 1U;
  }
  return y;
}

static bool positive(const I& x) { return x.leftBound() > R(0); }
static bool negative(const I& x) { return x.rightBound() < R(0); }
static bool excludesZero(const I& x) { return positive(x) || negative(x); }
static bool upperLess(const I& x, const char* b) {
  I B = IV(b);
  return x.rightBound() < B.leftBound();
}
static bool lowerGreater(const I& x, const char* b) {
  I B = IV(b);
  return x.leftBound() > B.rightBound();
}
static bool upperAtMost(const I& x, const char* b) {
  I B = IV(b);
  return x.rightBound() <= B.leftBound();
}
static bool lowerAtLeast(const I& x, const char* b) {
  I B = IV(b);
  return x.leftBound() >= B.rightBound();
}
static bool strictSubsetOfDecimalInterval(const I& x,
                                          const char* left,
                                          const char* right) {
  I L = IV(left);
  I Rb = IV(right);
  return x.leftBound() > L.rightBound() && x.rightBound() < Rb.leftBound();
}
static bool subsetOfDecimalInterval(const I& x,
                                    const char* left,
                                    const char* right) {
  return lowerAtLeast(x, left) && upperAtMost(x, right);
}

static I midpointInterval(const I& x) {
  return I((x.leftBound() + x.rightBound()) / R(2));
}

static I widen(const I& x, const I& eps) {
  // Widen by interval arithmetic rather than by raw endpoint arithmetic.
  // Earlier drafts formed
  //   I(x.leftBound()-eps.rightBound(), x.rightBound()+eps.rightBound())
  // directly in the scalar type.  In a double build that scalar subtraction can
  // round inward before the interval constructor sees the endpoint.  The form
  // below asks CAPD to do the addition as interval arithmetic, hence with the
  // directed rounding policy of the active interval type.
  return x + I(-eps.rightBound(), eps.rightBound());
}

static I absUpperAsInterval(const I& x) {
  R a = x.leftBound();
  R b = x.rightBound();
  if (a < R(0)) a = -a;
  if (b < R(0)) b = -b;
  return I(R(0), (a > b ? a : b));
}

static I nonnegativePartUpper(const I& x) {
  R b = x.rightBound();
  if (b < R(0)) return zero();
  return I(R(0), b);
}

static std::string pass(bool ok) { return ok ? "PASS" : "FAIL"; }

// Diagnostic interval width.
//
// A small width is not by itself a proof, and a proof does not require widths
// to be "small" in a floating-point sense.  Still, widths are essential diagnostic
// data: if an FGR lower bound is positive only because an enclosure was silently
// assumed narrow, the reader should see the enclosure size.  These routines are
// therefore used mainly for reporting.  The eigenvalue-bisection block also
// uses one width as an explicit guard, so the subtraction itself must be an
// interval operation.
static I intervalWidth(const I& x) {
  return I(x.rightBound()) - I(x.leftBound());
}

static void printWidth(const std::string& label, const I& x) {
  std::cout << label << " width = " << intervalWidth(x) << "\n";
}

// CAPD interval elementary wrappers.
//
// These call CAPD's interval versions of sqrt, exp, log, sin, cos.  For example,
// if x=[a,b], then isqrt(x) is an interval enclosing sqrt(y) for every y in
// [a,b].  Keeping the calls behind wrappers makes it explicit which special
// functions are supplied by CAPD's interval arithmetic and also gives a single
// adjustment point if a local CAPD installation exposes them under a different
// namespace.
static I isqrt(const I& x) { return sqrt(x); }
static I iexp(const I& x) { return exp(x); }
static I ilog(const I& x) { return log(x); }
static I isin(const I& x) { return sin(x); }
static I icos(const I& x) { return cos(x); }
static I icosh(const I& x) { return (iexp(x) + iexp(-x)) / two(); }

// -----------------------------------------------------------------------------
// Complex interval arithmetic
// -----------------------------------------------------------------------------

struct CI {
  I re;
  I im;
  CI() : re(zero()), im(zero()) {}
  CI(I r, I i = zero()) : re(r), im(i) {}
};

static CI operator+(const CI& a, const CI& b) { return CI(a.re + b.re, a.im + b.im); }
static CI operator-(const CI& a, const CI& b) { return CI(a.re - b.re, a.im - b.im); }
static CI operator-(const CI& a) { return CI(-a.re, -a.im); }
static CI operator*(const CI& a, const CI& b) {
  return CI(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static CI operator*(const CI& a, const I& b) { return CI(a.re * b, a.im * b); }
static CI operator*(const I& b, const CI& a) { return a * b; }
static CI operator/(const CI& a, const I& b) { return CI(a.re / b, a.im / b); }
static CI conj(const CI& a) { return CI(a.re, -a.im); }
static CI mulI(const CI& a) { return CI(-a.im, a.re); }
static I norm2(const CI& a) { return sqr(a.re) + sqr(a.im); }
static CI inv(const CI& a) {
  I d = norm2(a);
  if (!positive(d)) throw std::runtime_error("complex interval division by box containing 0");
  return conj(a) / d;
}
static CI div(const CI& a, const CI& b) { return a * inv(b); }

static std::ostream& operator<<(std::ostream& os, const CI& z) {
  os << "(" << z.re << ") + i*(" << z.im << ")";
  return os;
}

static CI cis(const I& x) { return CI(icos(x), isin(x)); }
static I absUpper(const I& x) { return absUpperAsInterval(x); }
static I absUpper(const CI& z) {
  // We need an upper bound for |z|, not a sharp interval evaluation of the
  // function sqrt(re^2+im^2).  For boxes very close to zero, CAPD's interval
  // multiplication can leave a tiny negative lower endpoint in re^2+im^2 due to
  // dependency/rounding representation.  Taking sqrt of that interval is
  // needlessly fragile.  Instead, explicitly form the nonnegative upper box
  //
  //   0 <= |z| <= sqrt( max|Re z|^2 + max|Im z|^2 ).
  //
  // This is the right object for all places where absUpper is used: tail
  // constants, residual diagnostics, and denominator lower-bound checks.
  I ar = absUpperAsInterval(z.re);
  I ai = absUpperAsInterval(z.im);
  I upper = sqr(I(ar.rightBound())) + sqr(I(ai.rightBound()));
  return isqrt(I(R(0), upper.rightBound()));
}

static I maxUpperInterval(I a, const I& b) {
  if (b.rightBound() > a.rightBound()) return I(R(0), b.rightBound());
  return I(R(0), a.rightBound());
}

// Complex interval width diagnostic.  Complex quantities are represented as
// rectangular boxes Re z in I_re, Im z in I_im.  The printed widths tell the
// reader how much interval over-enclosure remains in each rectangular direction.
static void printCIWidth(const std::string& label, const CI& z) {
  std::cout << label << " Re width = " << intervalWidth(z.re)
            << ", Im width = " << intervalWidth(z.im) << "\n";
}

// Lower bound for |z|^2 when z is a complex interval rectangle.
//
// Plain interval arithmetic applied to re^2+im^2 is always rigorous, but it can
// be needlessly pessimistic if one component interval crosses zero.  For the
// FGR positivity certificate we want the sharp elementary lower bound for the
// distance from the rectangle to the origin:
//   min_{x in re, y in im} (x^2+y^2).
// If a real interval contains zero, its contribution to this minimum is zero;
// otherwise CAPD's interval square gives a rigorous lower endpoint.
static bool containsZero(const I& x) {
  return !(positive(x) || negative(x));
}

static I squareLowerBound(const I& x) {
  if (containsZero(x)) return zero();
  I y = sqr(x);
  R lb = y.leftBound();
  if (lb < R(0)) lb = R(0);
  return I(lb);
}

static I norm2LowerBound(const CI& z) {
  return squareLowerBound(z.re) + squareLowerBound(z.im);
}

static void requireTrue(bool ok, const std::string& msg) {
  if (!ok) throw std::runtime_error("coded lemma failed: " + msg);
}

static bool containsInterval(const I& outer, const I& inner) {
  return outer.leftBound() <= inner.leftBound()
      && outer.rightBound() >= inner.rightBound();
}

static bool checkDecimalInputEnclosures() {
  std::cout << "\n[0] Directed decimal input check\n";
  struct DecimalRational {
    const char* decimal;
    const char* numerator;
    const char* denominator;
  };
  const DecimalRational tests[] = {
    {"0.1", "1", "10"},
    {"0.95", "19", "20"},
    {"0.6032878545810", "6032878545810", "10000000000000"},
    {"0.6032878545819", "6032878545819", "10000000000000"},
    {"0.6032878545816699", "6032878545816699", "10000000000000000"},
    {"0.6032878545816856", "6032878545816856", "10000000000000000"},
    {"0.77747", "77747", "100000"},
    {"0.77753", "77753", "100000"},
    {"0.777471875", "777471875", "1000000000"},
    {"0.77747375", "77747375", "100000000"},
    {"8.001", "8001", "1000"},
    {"0.0005", "1", "2000"},
    {"1.4", "7", "5"},
    {"0.0548", "137", "2500"},
    {"1e-15", "1", "1000000000000000"}
  };
  for (const auto& t : tests) {
    I parsed = IV(t.decimal);
    I rational = IV(t.numerator) / IV(t.denominator);
    requireTrue(containsInterval(parsed, rational),
                std::string("directed parsing must enclose exact decimal ") + t.decimal);
  }
  I nonDyadic = IV("0.1");
  I dyadic = IV("0.5");
  requireTrue(nonDyadic.leftBound() < nonDyadic.rightBound(),
              "non-dyadic decimal 0.1 must not become a singleton");
  requireTrue(dyadic.leftBound() == dyadic.rightBound(),
              "dyadic decimal 0.5 should be represented exactly");
  requireTrue(!upperLess(IV("0.0548"), "0.0548")
              && !lowerGreater(IV("0.0548"), "0.0548"),
              "strict decimal comparisons must not accept equality");
  I tenth = IV("1") / IV("10");
  I fifth = IV("1") / IV("5");
  I decimalBox = IV("0.1", "0.2");
  requireTrue(decimalBox.leftBound() <= tenth.leftBound()
              && decimalBox.rightBound() >= fifth.rightBound(),
              "two-endpoint decimal interval must be rounded outward");
  I below = IV("99") / IV("1000");
  I above = IV("101") / IV("1000");
  requireTrue(upperLess(below, "0.1") && lowerGreater(above, "0.1"),
              "strict decimal comparisons must accept separated values");
  std::cout << "Representative proof-relevant decimal literals are enclosed "
               "by directed parsing.\n";
  std::cout << "PASS\n";
  return true;
}

static I intervalFromUpper(const I& upper) {
  return I(R(0), upper.rightBound());
}

static I upperSingleton(const I& x) {
  return I(x.rightBound());
}

static I symmetricError(const I& eps) {
  I e = intervalFromUpper(eps);
  return I(-e.rightBound(), e.rightBound());
}

static std::vector<I> firstCoeffs(const std::vector<I>& coeffs, int lastIndex) {
  // Return coefficients 0,...,lastIndex and deliberately discard any storage
  // slots beyond lastIndex.  This matters for the vortex tangent q-coefficients:
  // the q/dq arrays are sometimes allocated with one spare entry for recurrence
  // convenience.  That spare entry is not a computed coefficient and must never
  // be used as the final term in a tail certificate.
  requireTrue(lastIndex >= 0, "coefficient prefix needs nonnegative last index");
  requireTrue(lastIndex < static_cast<int>(coeffs.size()), "coefficient prefix out of range");
  return std::vector<I>(coeffs.begin(), coeffs.begin() + lastIndex + 1);
}

// -----------------------------------------------------------------------------
// Settings and constants
// -----------------------------------------------------------------------------

struct Settings {
  static constexpr int default_precision_bits = 220;

  int precision_bits = default_precision_bits;
  int taylor_order = 40;
  int series_order = 80;

  I c_box = IV("0.6032878545810", "0.6032878545819");
  I c_lower = IV("0.6032878545810");
  I c_upper = IV("0.6032878545819");
  // Refined profile box used by all later profile-dependent checks.
  //
  // This is deliberately wider than the high-precision midpoint printed in
  // earlier exploratory diagnostics.  The CAPD shooting block below proves
  // that the finite-radius interval Newton image for U(20;c)-1 lies strictly
  // inside this box.  Keeping the literal Newton image visible in the report is
  // more useful than silently shrinking this interval to a decimal guess.
  I c_cert = IV("0.6032878545816699", "0.6032878545816856");
  I c_mid = IV("0.6032878545816734968619270444122872610288");
  // mu_box is the active spectral-parameter interval used by the K0 and FGR
  // blocks.  It starts as the historical high-precision decimal only so the
  // input report can show what value motivated the computation.  In v5 this
  // interval is not trusted: certifyAndInstallFgrEigenvalueBox(...) replaces it
  // by a validated Wronskian sign-change bracket before any block that depends on
  // the internal eigenvalue is executed.
  I mu_box = IV("0.77749759351635044", "0.77749759351635045");
  I mu_eigen_box = IV("0.77747", "0.77753");
  I mu_eigen_left = IV("0.77747");
  I mu_eigen_right = IV("0.77753");
  I mu_fgr_box = IV("0.777471875", "0.77747375");
  I mu_fgr_left = IV("0.777471875");
  I mu_fgr_right = IV("0.77747375");

  I r0 = IV("0.1");
  I r_shoot = IV("20");
  I r_refined_shoot = IV("30");
  I r_h1_a = isqrt(IV("6"));
  I r_h1_b = IV("4");
  I r_match = IV("10");
  I r_threshold_inf = IV("16");
  I r_psi_tail = IV("16");
  I r_k0_start = IV("8.001");
  I r_fgr = IV("16");
  I r_lt = IV("15");

  // Constants used by coded lemmas.  These are not accepted silently: every
  // place they enter has a routine that verifies the required inequality and
  // prints the resulting interval margin.
  I vortex_tail_A = IV("8");       // candidate in |1-U^2| <= A exp(-r)/sqrt(r)
  I fgr_source_C = IV("30");       // candidate source constant in FGR tail
  I k0_comparison_C = IV("1.4");   // candidate psi <= C K0(kappa r)
  I delta0 = IV("0.0005");

  // Meshes.  Tighten if any check fails due to wrapping or quadrature width.
  int h1_cells = 800;
  int lt_cells = 3000;
  int fgr_cells = 3000;
  int origin_cells = 80;
  int asym_order = 10;
};

// -----------------------------------------------------------------------------
// Power series at r=0
// -----------------------------------------------------------------------------

struct VortexCoeffs {
  std::vector<I> p; // U(r)=r sum p[n] r^(2n)
  std::vector<I> q; // a(r)=sum_{n>=1} q[n] r^(2n)
};

struct VortexCoeffDerivs {
  VortexCoeffs C;
  std::vector<I> dp;
  std::vector<I> dq;
};

static VortexCoeffs vortexCoeffs(const I& c, int nMax) {
  VortexCoeffs C;
  C.p.assign(nMax + 1, zero());
  C.q.assign(nMax + 2, zero());
  C.p[0] = c;
  C.q[0] = zero();
  C.q[1] = IV("0.25");
  for (int n = 1; n <= nMax; ++n) {
    if (n >= 2) {
      I s = zero();
      for (int j = 0; j <= n - 2; ++j) s += C.p[j] * C.p[n - 2 - j];
      C.q[n] = -s / IV(std::to_string(4 * n).c_str());
    }
    I s = zero();
    for (int j = 1; j <= n; ++j) s += C.q[j] * C.p[n - j];
    C.p[n] = -s / IV(std::to_string(2 * n).c_str());
  }
  return C;
}

static VortexCoeffDerivs vortexCoeffDerivs(const I& c, int nMax) {
  VortexCoeffDerivs D;
  D.C.p.assign(nMax + 1, zero());
  D.C.q.assign(nMax + 2, zero());
  D.dp.assign(nMax + 1, zero());
  D.dq.assign(nMax + 2, zero());

  D.C.p[0] = c;
  D.dp[0] = one();
  D.C.q[0] = zero();
  D.C.q[1] = IV("0.25");
  D.dq[0] = zero();
  D.dq[1] = zero();

  for (int n = 1; n <= nMax; ++n) {
    if (n >= 2) {
      I s = zero();
      I ds = zero();
      for (int j = 0; j <= n - 2; ++j) {
        int k = n - 2 - j;
        s += D.C.p[j] * D.C.p[k];
        ds += D.dp[j] * D.C.p[k] + D.C.p[j] * D.dp[k];
      }
      I den = IV(std::to_string(4 * n).c_str());
      D.C.q[n] = -s / den;
      D.dq[n] = -ds / den;
    }
    I s = zero();
    I ds = zero();
    for (int j = 1; j <= n; ++j) {
      int k = n - j;
      s += D.C.q[j] * D.C.p[k];
      ds += D.dq[j] * D.C.p[k] + D.C.q[j] * D.dp[k];
    }
    I den = IV(std::to_string(2 * n).c_str());
    D.C.p[n] = -s / den;
    D.dp[n] = -ds / den;
  }
  return D;
}

static I vortexUTailBound(const I& r, const I& cMax, int N) {
  // From the paper's coefficient majorants on 0 <= r <= 1/4:
  // |p_n| <= cMax*2^{-n}.  Tail after n=N.
  I rho = sqr(r) / two();
  return cMax * r * ipow(rho, static_cast<unsigned>(N + 1)) / (one() - rho);
}

static I vortexATailBound(const I& r, int N) {
  // |q_n| <= 2^{-n-1}; tail after n=N.
  I rho = sqr(r) / two();
  return half() * ipow(rho, static_cast<unsigned>(N + 1)) / (one() - rho);
}

static I vortexDerivativeFactorTailBound(const I& r, int N) {
  // Lemma Wbderror in the paper gives, for q=r^2/2,
  //   sum_{n>N} 4 n q^n
  //     = 4 q^(N+1) ((N+1)-N q)/(1-q)^2.
  // This bounds the differentiated Frobenius factors for W and b.
  I q = sqr(r) / two();
  I n = IV(std::to_string(N).c_str());
  I np1 = IV(std::to_string(N + 1).c_str());
  return IV("4") * ipow(q, static_cast<unsigned>(N + 1)) * (np1 - n * q) / sqr(one() - q);
}

static I vortexUTangentTailBound(const I& r, int N) {
  // U(r,c)=c r W(r,c).  The first term below is the W-tail contribution from
  // differentiating the leading c factor; the second is c*r*partial_c W with
  // c<=1/sqrt(2), as in Lemma Wbderror.
  I q = sqr(r) / two();
  I wTail = ipow(q, static_cast<unsigned>(N + 1)) / (one() - q);
  I factorTail = vortexDerivativeFactorTailBound(r, N);
  return r * (wTail + factorTail / isqrt(two()));
}

static I vortexATangentTailBound(const I& r, int N) {
  // a(r,c)=r^2 b(r,c)/4.
  return sqr(r) * vortexDerivativeFactorTailBound(r, N) / IV("4");
}

struct VortexState { I u, a, up; };
struct VortexVarState { I u, a, uc, ac; };

static VortexState vortexSeries(const I& r, const I& c, int nMax, bool addTail = true) {
  VortexCoeffs C = vortexCoeffs(c, nMax);
  I r2 = sqr(r);
  I sU = zero();
  I sA = zero();
  for (int n = 0; n <= nMax; ++n) sU += C.p[n] * ipow(r2, static_cast<unsigned>(n));
  for (int n = 1; n <= nMax; ++n) sA += C.q[n] * ipow(r2, static_cast<unsigned>(n));
  I u = r * sU;
  I a = sA;
  if (addTail) {
    I cAbs = absUpperAsInterval(c);
    u = widen(u, vortexUTailBound(r, cAbs, nMax));
    a = widen(a, vortexATailBound(r, nMax));
  }
  I up = ((one() - a) / r) * u;
  return {u, a, up};
}

static VortexVarState vortexSeriesVar(const I& r, const I& c, int nMax) {
  int Nwork = nMax + 20;
  VortexCoeffDerivs D = vortexCoeffDerivs(c, Nwork);
  I r2 = sqr(r);
  I sU = zero(), sA = zero(), sUc = zero(), sAc = zero();
  for (int n = 0; n <= Nwork; ++n) {
    I pow = ipow(r2, static_cast<unsigned>(n));
    sU += D.C.p[n] * pow;
    sUc += D.dp[n] * pow;
  }
  for (int n = 1; n <= Nwork; ++n) {
    I pow = ipow(r2, static_cast<unsigned>(n));
    sA += D.C.q[n] * pow;
    sAc += D.dq[n] * pow;
  }
  I u = r * sU;
  I a = sA;
  I uc = r * sUc;
  I ac = sAc;

  u = widen(u, vortexUTailBound(r, absUpperAsInterval(c), Nwork));
  a = widen(a, vortexATailBound(r, Nwork));
  uc = widen(uc, vortexUTangentTailBound(r, Nwork));
  ac = widen(ac, vortexATangentTailBound(r, Nwork));
  return {u, a, uc, ac};
}

// Polynomial helpers for series recurrences in x=r^2.
static std::vector<I> polyMul(const std::vector<I>& a, const std::vector<I>& b, int N) {
  std::vector<I> c(N + 1, zero());
  for (int i = 0; i <= N; ++i) {
    if (i >= static_cast<int>(a.size())) break;
    for (int j = 0; j + i <= N; ++j) {
      if (j >= static_cast<int>(b.size())) break;
      c[i + j] += a[i] * b[j];
    }
  }
  return c;
}

static std::vector<I> polyShiftX(const std::vector<I>& a, int N) {
  std::vector<I> b(N + 1, zero());
  for (int i = 0; i < N && i < static_cast<int>(a.size()); ++i) b[i + 1] = a[i];
  return b;
}

static I polyEval(const std::vector<I>& a, const I& x) {
  I y = zero();
  for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) y = y * x + a[i];
  return y;
}

static I polyDerivEval(const std::vector<I>& a, const I& x) {
  I y = zero();
  for (int i = static_cast<int>(a.size()) - 1; i >= 1; --i) y = y * x + IV(std::to_string(i).c_str()) * a[i];
  return y;
}

static std::vector<I> internalModeH(const I& c, const I& mu, int N) {
  VortexCoeffs VC = vortexCoeffs(c, N + 4);
  std::vector<I> f(N + 3, zero());
  for (int j = 0; j < static_cast<int>(f.size()) && j < static_cast<int>(VC.p.size()); ++j) f[j] = VC.p[j];
  std::vector<I> f2 = polyMul(f, f, N + 1);
  std::vector<I> h(N + 1, zero());
  h[0] = one();
  for (int n = 0; n <= N - 1; ++n) {
    // h_{n+1} = coeff_n(x f^2 h - mu h)/(4(n+1)^2)
    std::vector<I> fh = polyMul(f2, h, n); // enough for coeff n-1 before shift
    std::vector<I> xfh = polyShiftX(fh, n);
    I coeff = xfh[n] - mu * h[n];
    I den = IV(std::to_string(4 * (n + 1) * (n + 1)).c_str());
    h[n + 1] = coeff / den;
  }
  return h;
}

struct PsiState { I psi, psip; };

// Proof object for the Frobenius tail of the regular internal-mode solution.
// Besides the two radii, retain the complete parameter domain on which the
// all-order invariant-ball estimate was proved.  This prevents a tail bound
// from being reused silently outside its certified domain.
struct InternalFrobeniusTailCertificate {
  int M;
  I rBox;
  I cBox;
  I muBox;
  I valueTail;
  I xDerivativeTail;
};

static PsiState psiSeriesUsingCertifiedTail(
    const I& r, const I& c, const I& mu, int N,
    const InternalFrobeniusTailCertificate& cert) {
  const int M = N + 20;
  requireTrue(cert.M == M,
              "internal Frobenius certificate has the wrong truncation index");
  requireTrue(containsInterval(cert.rBox, r),
              "internal Frobenius certificate does not contain the starting radius");
  requireTrue(containsInterval(cert.cBox, c),
              "internal Frobenius certificate does not contain the vortex parameter");
  requireTrue(containsInterval(cert.muBox, mu),
              "internal Frobenius certificate does not contain the spectral parameter");

  std::vector<I> h = internalModeH(c, mu, M);
  I x = sqr(r);
  I psi = polyEval(h, x);
  I psip = two() * r * polyDerivEval(h, x);
  psi = widen(psi, cert.valueTail);
  // psi(r)=h(r^2), hence the radial derivative of the tail is bounded
  // by 2|r| times the certified x-derivative tail.
  psip = widen(psip,
               two() * absUpperAsInterval(r) * cert.xDerivativeTail);
  return {psi, psip};
}

static std::vector<I> phiEnergyH(const I& c, const I& energy, int which, int N) {
  VortexCoeffs VC = vortexCoeffs(c, N + 6);
  std::vector<I> f(N + 6, zero()), g(N + 6, zero());
  for (int j = 0; j < static_cast<int>(f.size()) && j < static_cast<int>(VC.p.size()); ++j) f[j] = VC.p[j];
  for (int j = 0; j < static_cast<int>(g.size()) && j < static_cast<int>(VC.q.size()); ++j) g[j] = VC.q[j];
  std::vector<I> f2 = polyMul(f, f, N + 2);

  I nu = which == 1 ? IV("2.5") : IV("0.5");
  I s0 = which == 1 ? IV("4") : zero();
  std::vector<I> sTerm(N + 1, zero());
  if (which == 1) {
    std::vector<I> twoMinusG(N + 1, zero());
    twoMinusG[0] = two();
    for (int i = 1; i <= N && i < static_cast<int>(g.size()); ++i) twoMinusG[i] = -g[i];
    sTerm = polyMul(twoMinusG, twoMinusG, N);
  }
  std::vector<I> v0(N + 1, zero());
  if (which == 1) {
    v0[0] = half();
    std::vector<I> xf2 = polyShiftX(f2, N);
    for (int i = 0; i <= N; ++i) v0[i] += half() * xf2[i];
  } else {
    v0 = polyShiftX(f2, N);
  }

  std::vector<I> h(N + 1, zero());
  h[0] = one();

  auto residualCoeff = [&](const std::vector<I>& hh, int n) -> I {
    I coeff = zero();
    // -2(2nu+1)x h'
    if (n >= 1) coeff += -two() * (two() * nu + one()) * IV(std::to_string(n).c_str()) * hh[n];
    // -4 x^2 h''
    if (n >= 2) coeff += -IV("4") * IV(std::to_string(n * (n - 1)).c_str()) * hh[n];
    // +(sTerm-s0)h
    std::vector<I> ss = sTerm;
    ss[0] -= s0;
    std::vector<I> prod1 = polyMul(ss, hh, n);
    coeff += prod1[n];
    // + x(v0-energy)h
    std::vector<I> vv = v0;
    vv[0] -= energy;
    std::vector<I> prod2 = polyShiftX(polyMul(vv, hh, n - 1 >= 0 ? n - 1 : 0), n);
    coeff += prod2[n];
    return coeff;
  };

  for (int n = 1; n <= N; ++n) {
    h[n] = zero();
    I b0 = residualCoeff(h, n);
    h[n] = one();
    I b1 = residualCoeff(h, n);
    I lin = b1 - b0;
    if (!excludesZero(lin)) throw std::runtime_error("singular phi-energy recurrence interval");
    h[n] = -b0 / lin;
  }
  return h;
}

struct PhiState { I phi, phip; };

struct PhiFrobeniusTailCertificate {
  int M;
  int which;
  I rBox;
  I cBox;
  I energyBox;
  I valueTail;
  I xDerivativeTail;
};

static PhiState phiEnergySeriesUsingCertifiedTail(
    const I& r, const I& c, const I& energy, int which, int N,
    const PhiFrobeniusTailCertificate& cert) {
  int HN = N / 2 + 20;
  requireTrue(cert.M == HN && cert.which == which,
              "Phi Frobenius certificate has the wrong series or truncation index");
  requireTrue(containsInterval(cert.rBox, r),
              "Phi Frobenius certificate does not contain the starting radius");
  requireTrue(containsInterval(cert.cBox, c),
              "Phi Frobenius certificate does not contain the vortex parameter");
  requireTrue(containsInterval(cert.energyBox, energy),
              "Phi Frobenius certificate does not contain the energy parameter");
  std::vector<I> h = phiEnergyH(c, energy, which, HN);
  I x = sqr(r);
  I hx = polyEval(h, x);
  I hxp = polyDerivEval(h, x);
  I nu = which == 1 ? IV("2.5") : IV("0.5");
  I rnu = exp(nu * log(r));
  I phi = rnu * hx;
  I phip = nu * exp((nu - one()) * log(r)) * hx + rnu * two() * r * hxp;
  // Phi(r)=r^nu h(r^2).  Both factors contribute to Phi':
  //   Phi'_tail <= nu r^(nu-1) valueTail + r^nu 2r (x-derivative tail).
  // This is the precise radial chain rule missing from the earlier draft.
  I valueTail = rnu * cert.valueTail;
  I derivativeTail =
      nu * exp((nu - one()) * log(r)) * cert.valueTail
      + rnu * two() * absUpperAsInterval(r) * cert.xDerivativeTail;
  phi = widen(phi, valueTail);
  phip = widen(phip, derivativeTail);
  return {phi, phip};
}

static void printVortexOriginTailCertificates(const Settings& S) {
  std::cout << "origin vortex/tangent tail certificates at r0=" << S.r0 << "\n";
  int N = S.series_order + 20;
  std::cout << "  U,a value tails use Lemma Wberror coefficient majorants.\n";
  std::cout << "  U_c,a_c tangent tails use Lemma Wbderror derivative majorants.\n";
  std::cout << "  differentiated factor tail E_N = "
            << vortexDerivativeFactorTailBound(S.r0, N) << " with N=" << N << "\n";
  std::cout << "  U_c value tail bound = " << vortexUTangentTailBound(S.r0, N) << "\n";
  std::cout << "  a_c value tail bound = " << vortexATangentTailBound(S.r0, N) << "\n";
}

static void printSpectralOriginTailRadii(
    const Settings& S,
    const InternalFrobeniusTailCertificate& psi,
    const PhiFrobeniusTailCertificate& phi1,
    const PhiFrobeniusTailCertificate& phi2) {
  std::cout << "origin spectral Frobenius tail radii at r0=" << S.r0 << "\n";
  std::cout << "  These are the radii carried by the all-order invariant-ball\n"
               "  certificates used by the subsequent ODE calculations.\n";
  std::cout << "  psi: M=" << psi.M << ", value tail=" << psi.valueTail
            << ", x-derivative tail=" << psi.xDerivativeTail << "\n";
  std::cout << "  Phi_1: M=" << phi1.M << ", value tail=" << phi1.valueTail
            << ", x-derivative tail=" << phi1.xDerivativeTail << "\n";
  std::cout << "  Phi_2: M=" << phi2.M << ", value tail=" << phi2.valueTail
            << ", x-derivative tail=" << phi2.xDerivativeTail << "\n";
}

struct TailBootstrapRecord {
  std::string name;
  int M;
  I lastTerm;
  I rho;
  I maxCheckedRatio;
  I valueTail;
  I xDerivativeTail;
};

struct InfiniteTailMapRecord {
  std::string name;
  int M;
  int checkedSlices;
  I rho;
  I finiteMaxRatio;
  I headLeakage;
  I tailOperatorNorm;
  I totalNorm;
  I valueTail;
  I xDerivativeTail;
};

static I weightedTerm(const std::vector<I>& coeffs, int n, const I& x) {
  requireTrue(n >= 0 && n < static_cast<int>(coeffs.size()), "weighted term index out of range");
  return upperSingleton(absUpperAsInterval(coeffs[n]) * ipow(x, static_cast<unsigned>(n)));
}

static I geomSum1(const I& a) {
  return a / (one() - a);
}

static I geomSumN(const I& a) {
  return a / sqr(one() - a);
}

static I geomSumNMinus1(const I& a) {
  return sqr(a) / sqr(one() - a);
}

static I geomSumKPlus1From0(const I& a) {
  return one() / sqr(one() - a);
}

static TailBootstrapRecord checkComputedTailBootstrapSlices(const std::string& name,
                                                  const std::vector<I>& coeffs,
                                                  const I& x,
                                                  int M,
                                                  int slices,
                                                  const char* rhoText) {
  requireTrue(M >= 8, name + ": geometric tail check needs M>=8");
  requireTrue(M + slices < static_cast<int>(coeffs.size()),
              name + ": not enough extra coefficients for geometric tail slice check");
  I rho = IV(rhoText);
  requireTrue(upperLess(rho, "1"), name + ": rho must be <1");
  I last = weightedTerm(coeffs, M, x);
  I maxRatio = zero();
  for (int m = 1; m <= slices; ++m) {
    I term = weightedTerm(coeffs, M + m, x);
    I allowed = last * ipow(rho, static_cast<unsigned>(m));
    // The proposed ball uses the exact decimal rho.  Comparing with the
    // downward-rounded endpoint of `allowed` prevents outward rounding of the
    // right endpoint from making this finite inclusion pass spuriously.
    requireTrue(term.rightBound() <= allowed.leftBound(),
                name + ": computed tail coefficient not in proposed ball");
    I prev = weightedTerm(coeffs, M + m - 1, x);
    if (prev.rightBound() > R(0)) {
      I ratio = term / prev;
      if (ratio.rightBound() > maxRatio.rightBound()) maxRatio = ratio;
    }
  }
  I oneMinusRho = one() - rho;
  I valueTail = last * rho / oneMinusRho;
  I weightedSum =
      IV(std::to_string(M).c_str()) * rho / oneMinusRho + rho / sqr(oneMinusRho);
  I xDerivativeTail = last * weightedSum / x;
  return {name, M, last, rho, maxRatio, valueTail, xDerivativeTail};
}

static void printTailBootstrapRecord(const TailBootstrapRecord& r) {
  std::cout << "  " << r.name << ": M=" << r.M
            << ", T_M=" << r.lastTerm
            << ", rho=" << r.rho
            << ", max checked ratio=" << r.maxCheckedRatio << "\n";
  std::cout << "    value tail from ball = " << r.valueTail
            << ", x-derivative tail from ball = " << r.xDerivativeTail << "\n";
}

static void printInfiniteTailMapRecord(const InfiniteTailMapRecord& r) {
  std::cout << "  " << r.name << " infinite tail self-map: M=" << r.M
            << ", finite slices=" << r.checkedSlices
            << ", rho=" << r.rho
            << ", max finite ratio=" << r.finiteMaxRatio << "\n";
  std::cout << "    head leakage = " << r.headLeakage
            << ", tail operator norm = " << r.tailOperatorNorm
            << ", total = " << r.totalNorm << "\n";
  std::cout << "    value tail from invariant ball = " << r.valueTail
            << ", x-derivative tail from invariant ball = " << r.xDerivativeTail << "\n";
}

static I weightedF2Majorant(int k, const I& cAbs, const I& sigma) {
  // If |p_n| x_*^n <= |c| sigma^n, then the weighted coefficients of f^2 obey
  // |(f^2)_k| x_*^k <= |c|^2 (k+1) sigma^k.
  requireTrue(k >= 0, "f^2 majorant index must be nonnegative");
  return sqr(cAbs) * IV(std::to_string(k + 1).c_str()) * ipow(sigma, static_cast<unsigned>(k));
}

static I phiSMajorant(int which, int i, const I& sigma) {
  requireTrue(i >= 1, "Phi S majorant index must be positive");
  if (which == 2) return zero();
  // For Phi_1, S=(2-a)^2 and S_0=4 is removed.  Since
  // |q_i|x_*^i <= (1/2)sigma^i, i>=1,
  // |S_i|x_*^i <= 2 sigma^i + (i-1) sigma^i/4.
  I ii = IV(std::to_string(i).c_str());
  return (two() + (ii - one()) / IV("4")) * ipow(sigma, static_cast<unsigned>(i));
}

static I phiWMajorant(int which, int i, const I& cAbs, const I& sigma) {
  requireTrue(i >= 1, "Phi W majorant index must be positive");
  I ii = IV(std::to_string(i).c_str());
  if (which == 1) {
    // W_i=(1/2)(f^2)_{i-1} after the extra x shift, hence
    // |W_i|x_*^i <= |c|^2 i sigma^i because sigma=x_*/2.
    return sqr(cAbs) * ii * ipow(sigma, static_cast<unsigned>(i));
  }
  // Phi_2 has W_i=(f^2)_{i-1} after the x shift.
  return two() * sqr(cAbs) * ii * ipow(sigma, static_cast<unsigned>(i));
}

static I phiDenominator(int which, int n) {
  requireTrue(n >= 1, "Phi denominator index must be positive");
  I nn = IV(std::to_string(n).c_str());
  if (which == 1) return IV("4") * sqr(nn) + IV("8") * nn;
  return IV("4") * sqr(nn);
}

static InfiniteTailMapRecord checkInternalInfiniteTailSelfMap(const std::string& name,
                                                              const std::vector<I>& coeffs,
                                                              const I& x,
                                                              const I& cAbs,
                                                              const I& muAbs,
                                                              int M,
                                                              int slices,
                                                              const char* rhoText) {
  TailBootstrapRecord finite = checkComputedTailBootstrapSlices(name, coeffs, x, M, slices, rhoText);
  I rho = finite.rho;
  I sigma = x / two();
  I ratio = sigma / rho;
  requireTrue(upperLess(ratio, "1"), name + ": kernel decay must be below tail rho");

  int m0 = slices + 1;
  int j0i = M + m0;
  I j0 = IV(std::to_string(j0i).c_str());
  I den = IV("4") * sqr(j0);
  I T = finite.lastTerm;
  I rhoM = ipow(rho, static_cast<unsigned>(m0));

  I head = zero();
  for (int l = 0; l <= M; ++l) {
    int k = M + m0 - 2 - l;
    if (k >= 0) head += weightedTerm(coeffs, l, x) * weightedF2Majorant(k, cAbs, sigma);
  }
  I headLeakage = sqr(x) * head / (den * T * rhoM);
  I muNorm = muAbs * x / (den * rho);
  I convolutionNorm =
      sqr(x) * sqr(cAbs) * geomSumKPlus1From0(ratio) / (den * sqr(rho));
  I tailNorm = muNorm + convolutionNorm;
  I total = headLeakage + tailNorm;
  requireTrue(upperLess(total, "1"), name + ": infinite Frobenius tail self-map norm < 1");

  return {name, M, slices, rho, finite.maxCheckedRatio, headLeakage, tailNorm,
          total, finite.valueTail, finite.xDerivativeTail};
}

static InfiniteTailMapRecord checkPhiInfiniteTailSelfMap(const std::string& name,
                                                         const std::vector<I>& coeffs,
                                                         const I& x,
                                                         const I& cAbs,
                                                         const I& energy,
                                                         int which,
                                                         int M,
                                                         int slices,
                                                         const char* rhoText) {
  TailBootstrapRecord finite = checkComputedTailBootstrapSlices(name, coeffs, x, M, slices, rhoText);
  I rho = finite.rho;
  I sigma = x / two();
  I ratio = sigma / rho;
  requireTrue(upperLess(ratio, "1"), name + ": kernel decay must be below tail rho");

  int m0 = slices + 1;
  int j0i = M + m0;
  I den = phiDenominator(which, j0i);
  I T = finite.lastTerm;
  I rhoM = ipow(rho, static_cast<unsigned>(m0));

  I headS = zero();
  I headW = zero();
  for (int l = 0; l <= M; ++l) {
    int iS = M + m0 - l;
    if (iS >= 1) headS += weightedTerm(coeffs, l, x) * phiSMajorant(which, iS, sigma);
    int iW = M + m0 - 1 - l;
    if (iW >= 1) headW += weightedTerm(coeffs, l, x) * phiWMajorant(which, iW, cAbs, sigma);
  }
  I headLeakage = (headS + x * headW) / (den * T * rhoM);

  I tailS = zero();
  if (which == 1) {
    tailS = (two() * geomSum1(ratio) + geomSumNMinus1(ratio) / IV("4")) / den;
  }
  I w0 = which == 1 ? absUpperAsInterval(half() - energy) : absUpperAsInterval(-energy);
  I wTailConstant = which == 1 ? sqr(cAbs) : two() * sqr(cAbs);
  I tailW0 = x * w0 / (den * rho);
  I tailW = x * wTailConstant * geomSumN(ratio) / (den * rho);
  I tailNorm = tailS + tailW0 + tailW;
  I total = headLeakage + tailNorm;
  requireTrue(upperLess(total, "1"), name + ": infinite Frobenius tail self-map norm < 1");

  return {name, M, slices, rho, finite.maxCheckedRatio, headLeakage, tailNorm,
          total, finite.valueTail, finite.xDerivativeTail};
}

static void checkVortexValueMajorantsAtOrigin(const Settings& S) {
  I x = sqr(S.r0);
  I cAbs = upperSingleton(absUpperAsInterval(S.c_cert));
  I sigma = x / two();
  requireTrue(upperLess(sigma, "0.01"), "vortex origin majorant sigma too large");
  // The paper's origin majorants are
  //   |p_n| <= |c| 2^{-n},   |q_n| <= 2^{-n-1}.
  // At x=r0^2 this means |p_n|x^n <= |c| sigma^n and
  // |q_n|x^n <= (1/2) sigma^n.  The induction closes because
  // |c|^2 < 1/2 and the p-recurrence gains a harmless factor <=1/4.
  requireTrue((sqr(cAbs)).rightBound() < half().leftBound(),
              "vortex p/q coefficient majorant requires |c|^2<1/2");
  VortexCoeffs C = vortexCoeffs(S.c_cert, S.series_order + 20);
  I quarter = IV("0.25");
  requireTrue(C.q[1].leftBound() == quarter.leftBound()
              && C.q[1].rightBound() == quarter.rightBound(),
              "the exact leading a coefficient must be q_1=1/4");
  for (int n = 0; n <= S.series_order + 20; ++n) {
    I pAllowed = cAbs * ipow(sigma, static_cast<unsigned>(n));
    requireTrue(weightedTerm(C.p, n, x).rightBound() <= pAllowed.leftBound(),
                "computed U Frobenius coefficient exceeds analytic majorant");
    // For n=1 the bound is the exact identity q_1 x=x/4, checked above at
    // the coefficient level.  Comparing two separately rounded evaluations
    // of x/4 would turn that equality into overlapping narrow intervals.
    if (n >= 2) {
      I qAllowed = half() * ipow(sigma, static_cast<unsigned>(n));
      requireTrue(weightedTerm(C.q, n, x).rightBound() <= qAllowed.leftBound(),
                  "computed a Frobenius coefficient exceeds analytic majorant");
    }
  }
}

static void checkVortexDerivativeMajorantsAtOrigin(const Settings& S) {
  I x = sqr(S.r0);
  VortexCoeffDerivs D = vortexCoeffDerivs(S.c_cert, S.series_order + 40);
  for (int n = 1; n <= S.series_order + 40; ++n) {
    I allowed = IV("4") * IV(std::to_string(n).c_str()) * ipow(x / two(), static_cast<unsigned>(n));
    requireTrue(weightedTerm(D.dp, n, x).rightBound() <= allowed.leftBound(),
                "computed partial_c U Frobenius coefficient exceeds derivative majorant");
    requireTrue(weightedTerm(D.dq, n, x).rightBound() <= allowed.leftBound(),
                "computed partial_c a Frobenius coefficient exceeds derivative majorant");
  }
}

static bool checkInternalModeTailBootstrapData(
    const Settings& S, const I& muBox, const std::string& label,
    InternalFrobeniusTailCertificate& certificate) {
  std::cout << "\n[Origin Frobenius series-tail data: " << label << "]\n";
  I x = sqr(S.r0);
  int slices = 20;
  I cAbs = upperSingleton(absUpperAsInterval(S.c_cert));
  I muAbs = upperSingleton(absUpperAsInterval(muBox));
  std::vector<I> psiH = internalModeH(S.c_cert, muBox, S.series_order + 20 + slices);
  InfiniteTailMapRecord psi = checkInternalInfiniteTailSelfMap(
      "psi " + label, psiH, x, cAbs, muAbs, S.series_order + 20, slices, "0.95");
  requireTrue(upperLess(psi.totalNorm, "0.358538353958974"),
              "broad psi Frobenius beta <= 0.358538353958974");
  requireTrue(upperLess(psi.valueTail, "2.74366499059244e-279"),
              "broad psi Frobenius value tail bound");
  requireTrue(upperLess(psi.xDerivativeTail, "3.29239798871093e-275"),
              "broad psi Frobenius x-derivative tail bound");
  requireTrue(upperLess(
                  two() * absUpperAsInterval(S.r0) * psi.xDerivativeTail,
                  "6.585e-276"),
              "broad psi Frobenius radial-derivative tail bound");
  certificate = {S.series_order + 20, S.r0, S.c_cert, muBox,
                 psi.valueTail, psi.xDerivativeTail};
  printInfiniteTailMapRecord(psi);
  std::cout << "PASS\n";
  return true;
}

static bool checkThresholdOriginTailBootstrapData(
    const Settings& S, PhiFrobeniusTailCertificate& certificate) {
  std::cout << "\n[4a] Threshold origin Frobenius series-tail data\n";
  I x = sqr(S.r0);
  int slices = 20;
  int phiM = S.series_order / 2 + 20;
  I cAbs = upperSingleton(absUpperAsInterval(S.c_cert));
  std::vector<I> phi2H = phiEnergyH(S.c_cert, one(), 2, phiM + slices);
  InfiniteTailMapRecord phi2 = checkPhiInfiniteTailSelfMap(
      "threshold Phi_2^(0)", phi2H, x, cAbs, one(), 2, phiM, slices, "0.95");
  requireTrue(upperLess(phi2.totalNorm, "4.027e-7"),
              "threshold Phi_2 Frobenius beta < 4.027e-7");
  requireTrue(upperLess(phi2.valueTail, "3.122e-167"),
              "threshold Phi_2 Frobenius value tail bound");
  requireTrue(upperLess(phi2.xDerivativeTail, "2.497e-163"),
              "threshold Phi_2 Frobenius x-derivative tail bound");
  certificate = {phiM, 2, S.r0, S.c_cert, one(),
                 phi2.valueTail, phi2.xDerivativeTail};
  printInfiniteTailMapRecord(phi2);
  std::cout << "PASS\n";
  return true;
}

static bool checkFrobeniusTailBootstrapData(
    const Settings& S,
    InternalFrobeniusTailCertificate& psiCertificate,
    PhiFrobeniusTailCertificate& phi1Certificate,
    PhiFrobeniusTailCertificate& phi2Certificate) {
  std::cout << "origin Frobenius series-tail bootstrap data\n";
  std::cout << "  This block verifies the finite recurrence data used by\n"
               "  Lemma capd_origin_frobenius.  The vortex U,a tails use\n"
               "  the analytic coefficient majorants; the shooting derivative\n"
               "  tails use Lemma Wbderror.  For psi, Phi_1, Phi_2 this\n"
               "  block checks finite slices and the closed-form infinite\n"
               "  weighted geometric tail bootstrap.\n";

  checkVortexValueMajorantsAtOrigin(S);
  std::cout << "  U,a analytic coefficient majorants at r0: PASS\n";
  checkVortexDerivativeMajorantsAtOrigin(S);
  std::cout << "  partial_c U,a analytic derivative majorants at r0: PASS\n";

  I x = sqr(S.r0);
  int slices = 20;

  VortexCoeffDerivs D = vortexCoeffDerivs(S.c_cert, S.series_order + 20 + slices);
  std::vector<I> dqComputed = firstCoeffs(D.dq, S.series_order + 20 + slices);
  TailBootstrapRecord uc = checkComputedTailBootstrapSlices("partial_c U", D.dp, x,
                                                  S.series_order + 20, slices, "0.95");
  TailBootstrapRecord ac = checkComputedTailBootstrapSlices("partial_c a", dqComputed, x,
                                                  S.series_order + 20, slices, "0.95");
  printTailBootstrapRecord(uc);
  printTailBootstrapRecord(ac);

  I cAbs = upperSingleton(absUpperAsInterval(S.c_cert));
  I muAbs = upperSingleton(absUpperAsInterval(S.mu_box));
  std::vector<I> psiH = internalModeH(S.c_cert, S.mu_box, S.series_order + 20 + slices);
  InfiniteTailMapRecord psi = checkInternalInfiniteTailSelfMap(
      "psi", psiH, x, cAbs, muAbs, S.series_order + 20, slices, "0.95");
  printInfiniteTailMapRecord(psi);

  I energy = IV("4") * S.mu_box;
  int phiM = S.series_order / 2 + 20;
  std::vector<I> phi1H = phiEnergyH(S.c_cert, energy, 1, phiM + slices);
  std::vector<I> phi2H = phiEnergyH(S.c_cert, energy, 2, phiM + slices);
  InfiniteTailMapRecord phi1 = checkPhiInfiniteTailSelfMap(
      "Phi_1", phi1H, x, cAbs, energy, 1, phiM, slices, "0.95");
  InfiniteTailMapRecord phi2 = checkPhiInfiniteTailSelfMap(
      "Phi_2", phi2H, x, cAbs, energy, 2, phiM, slices, "0.95");
  requireTrue(upperLess(psi.totalNorm, "0.358546")
              && upperLess(psi.valueTail, "2.744e-279")
              && upperLess(psi.xDerivativeTail, "3.293e-275"),
              "psi Frobenius table bounds");
  requireTrue(upperLess(phi1.totalNorm, "1.417e-6")
              && upperLess(phi1.valueTail, "2.948e-167")
              && upperLess(phi1.xDerivativeTail, "2.358e-163"),
              "Phi_1 Frobenius table bounds");
  requireTrue(upperLess(phi2.totalNorm, "1.249e-6")
              && upperLess(phi2.valueTail, "1.096e-166")
              && upperLess(phi2.xDerivativeTail, "8.762e-163"),
              "Phi_2 Frobenius table bounds");
  psiCertificate = {S.series_order + 20, S.r0, S.c_cert, S.mu_box,
                    psi.valueTail, psi.xDerivativeTail};
  phi1Certificate = {phiM, 1, S.r0, S.c_cert, energy,
                     phi1.valueTail, phi1.xDerivativeTail};
  phi2Certificate = {phiM, 2, S.r0, S.c_cert, energy,
                     phi2.valueTail, phi2.xDerivativeTail};
  printInfiniteTailMapRecord(phi1);
  printInfiniteTailMapRecord(phi2);
  std::cout << "PASS\n";
  return true;
}

// -----------------------------------------------------------------------------
// CAPD ODE wrappers
// -----------------------------------------------------------------------------

// Validated endpoint propagation.
//
// Inputs:
//   vf    CAPD vector-field string, parsed by IMap/MpIMap.
//   x0    interval box of initial data.
//   T     interval time.  In this program T is normally a singleton interval.
//   order Taylor order used by CAPD's interval Taylor solver.
//
// CAPD work performed here:
//   1. Map map(vf) parses the vector field and prepares interval Taylor
//      arithmetic for it.
//   2. Solver solver(map,order) constructs the validated Taylor ODE solver.
//   3. Set set(x0) embeds the initial box in CAPD's C0TripletonSet format.
//   4. TimeMap tm(solver) advances the whole set to time T.
//
// Return value:
//   An interval vector enclosing every exact solution at final time T, for all
//   initial data in x0 and all parameters included in x0.
static Vec integrate(const std::string& vf, const Vec& x0, const I& T, int order) {
  Map map(vf);
  Solver solver(map, order);
  TimeMap tm(solver);
  Set set(x0);
  return tm(T, set);
}

struct StepResult {
  Vec endpoint;
  Vec tube;
};

// Validated endpoint plus tube propagation.
//
// This is the same CAPD computation as integrate(...), but we also request a
// SolutionCurve from ITimeMap.  CAPD stores a rigorous Taylor enclosure of the
// flow along the integration interval.  Evaluating that curve on [0,T] produces
// an interval vector enclosing all intermediate values, not merely the endpoint.
//
// We use the tube for estimates that integrate or maximize quantities on a
// rectangle/cell.  This is essential: a bound obtained only at the cell endpoint
// would not be a proof for the whole cell.
static StepResult integrateTube(const std::string& vf, const Vec& x0, const I& T, int order) {
  Map map(vf);
  Solver solver(map, order);
  TimeMap tm(solver);
  TimeMap::SolutionCurve curve(I(R(0)));
  Set set(x0);
  Vec endpoint = tm(T, set, curve);
  I domain(R(0), T.rightBound());
  Vec tube = curve(domain);
  return {endpoint, tube};
}

static const std::string VF_VORTEX_FWD =
  "var:r,u,a,c;"
  "fun:1,((1-a)/r)*u,r*(1-u^2)/2,0;";

static const std::string VF_VORTEX_BWD =
  "var:r,u,a,c;"
  "fun:-1,-((1-a)/r)*u,-r*(1-u^2)/2,0;";

static const std::string VF_VORTEX_VAR_FWD =
  "var:r,u,a,uc,ac,c;"
  "fun:1,((1-a)/r)*u,r*(1-u^2)/2,"
  "((1-a)/r)*uc-(u/r)*ac,-r*u*uc,0;";

static const std::string VF_THRESHOLD_FWD =
  "var:r,u,a,f,fp,c;"
  "fun:1,((1-a)/r)*u,r*(1-u^2)/2,fp,(-1/(4*r^2)+u^2-1)*f,0;";

static const std::string VF_THRESHOLD_BWD =
  "var:r,u,a,f,fp,c;"
  "fun:-1,-((1-a)/r)*u,-r*(1-u^2)/2,-fp,-(-1/(4*r^2)+u^2-1)*f,0;";

static const std::string VF_INTERNAL_FWD =
  "var:r,u,a,psi,psip,c,mu;"
  "fun:1,((1-a)/r)*u,r*(1-u^2)/2,psip,-psip/r+(u^2-mu)*psi,0,0;";

static const std::string VF_INTERNAL_BWD =
  "var:r,u,a,psi,psip,c,mu;"
  "fun:-1,-((1-a)/r)*u,-r*(1-u^2)/2,-psip,psip/r-(u^2-mu)*psi,0,0;";

static const std::string VF_ZERO_FGR_FWD =
  "var:r,u,a,psi,psip,p1,p1p,p2,p2p,c,mu;"
  "fun:1,((1-a)/r)*u,r*(1-u^2)/2,"
  "psip,-psip/r+(u^2-mu)*psi,"
  "p1p,(-1/(4*r^2)+(1+u^2)/2+(2-a)^2/r^2-4*mu)*p1,"
  "p2p,(-1/(4*r^2)+u^2-4*mu)*p2,0,0;";

static std::string vfOutgoing(int which, bool backward) {
  std::ostringstream ss;
  ss << "var:r,u,a,sr,si,spr,spi,c,mu;fun:";
  if (!backward) {
    ss << "1,((1-a)/r)*u,r*(1-u^2)/2,";
    ss << "spr,spi,";
  } else {
    ss << "-1,-((1-a)/r)*u,-r*(1-u^2)/2,";
    ss << "-spr,-spi,";
  }
  std::string coef;
  if (which == 1)
    coef = "(-1/(4*r^2)+(1+u^2)/2+(2-a)^2/r^2-4*mu)";
  else
    coef = "(-1/(4*r^2)+u^2-4*mu)";
  if (!backward) {
    ss << coef << "*sr," << coef << "*si,0,0;";
  } else {
    ss << "-" << coef << "*sr,-" << coef << "*si,0,0;";
  }
  return ss.str();
}

static Vec initialVortex(const Settings& S, const I& c) {
  VortexState v = vortexSeries(S.r0, c, S.series_order);
  Vec x(4);
  x[0] = S.r0; x[1] = v.u; x[2] = v.a; x[3] = c;
  return x;
}

static void writeIntervalCsv(std::ostream& os, const I& x) {
  os << x.leftBound() << "," << x.rightBound();
}

static void emitVortexProfileTubeCsv(const Settings& S,
                                     const std::string& path,
                                     int cells) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("could not open vortex tube CSV: " + path);

  out << std::setprecision(40);
  out << "# backend="
#ifdef LPPSSI_USE_MP
      << "MpInterval/MpFloat"
#else
      << "DInterval/double"
#endif
      << "\n";
  out << "# precision_bits=" << S.precision_bits << "\n";
  out << "# taylor_order=" << S.taylor_order << "\n";
  out << "# c_cert,";
  writeIntervalCsv(out, S.c_cert);
  out << "\n";
  out << "# r_origin," << S.r0.leftBound() << "," << S.r0.rightBound() << "\n";
  out << "# r_tail," << S.r_fgr.leftBound() << "," << S.r_fgr.rightBound() << "\n";

  out << "cell,r_left,r_right,"
         "endpoint_U_lo,endpoint_U_hi,endpoint_a_lo,endpoint_a_hi,"
         "tube_U_lo,tube_U_hi,tube_a_lo,tube_a_hi\n";

  Vec x = initialVortex(S, S.c_cert);
  I h = (S.r_fgr - S.r0) / IV(std::to_string(cells).c_str());
  for (int j = 0; j < cells; ++j) {
    I rLeft = x[0];
    StepResult st = integrateTube(VF_VORTEX_FWD, x, h, S.taylor_order);
    I rRight = st.endpoint[0];

    out << j << ",";
    out << rLeft.leftBound() << "," << rRight.rightBound() << ",";
    writeIntervalCsv(out, st.endpoint[1]); out << ",";
    writeIntervalCsv(out, st.endpoint[2]); out << ",";
    writeIntervalCsv(out, st.tube[1]); out << ",";
    writeIntervalCsv(out, st.tube[2]); out << "\n";

    x = st.endpoint;
  }
}

static Vec initialVortexVar(const Settings& S, const I& c) {
  VortexVarState v = vortexSeriesVar(S.r0, c, S.series_order);
  Vec x(6);
  x[0] = S.r0; x[1] = v.u; x[2] = v.a; x[3] = v.uc; x[4] = v.ac; x[5] = c;
  return x;
}

// -----------------------------------------------------------------------------
// Bessel-K0/K1 rigorous bounds from the integral representation
//
// K_nu(x) = int_0^infty exp(-x cosh(t)) cosh(nu t) dt, nu=0,1.
// For x>1, both integrands are decreasing in t.  We use left/right Riemann
// sums on [0,T] plus a closed-form tail bound.  This avoids relying on a
// non-interval libm Bessel function.
// -----------------------------------------------------------------------------

static I besselK01_integral_bound(const I& x, int nu, int N = 12000, const I& T = IV("8")) {
  if (!(nu == 0 || nu == 1)) {
    throw std::runtime_error("besselK01_integral_bound supports only nu=0 or nu=1");
  }
  if (!(x.leftBound() > R(1))) {
    throw std::runtime_error("besselK01_integral_bound requires x.leftBound()>1");
  }
  I h = T / IV(std::to_string(N).c_str());
  I lower = zero();
  I upper = zero();
  I xl(x.leftBound());
  I xu(x.rightBound());
  for (int i = 0; i < N; ++i) {
    I ti = h * IV(std::to_string(i).c_str());
    I tip1 = h * IV(std::to_string(i + 1).c_str());
    I cu = icosh(ti);
    I cl = icosh(tip1);
    I gu = iexp(-xl * cu);
    I gl = iexp(-xu * cl);
    if (nu == 1) {
      gu *= cu;
      gl *= cl;
    }
    upper += h * gu;
    lower += h * gl;
  }
  // For t>=T, cosh(t) >= exp(t)/2 and cosh(t)<=exp(t).  For nu=0 the K0
  // tail is smaller than the K1-style bound below as cosh(t)>=1.
  I tail = two() / x * iexp(-x * iexp(T) / two());
  return I(lower.leftBound(), (upper + tail).rightBound());
}

static I K0(const I& x) { return besselK01_integral_bound(x, 0); }
static I K1(const I& x) { return besselK01_integral_bound(x, 1); }

// K0 upper bound used for tails: cosh(t) >= 1+t^2/2 gives
// K0(x) <= sqrt(pi/(2x)) exp(-x).
static I K0_gaussian_upper(const I& x) {
  return isqrt(piI() / (two() * x)) * iexp(-x);
}

// From Lemma lem_1-U2: for the true asymptotic vortex and R>=4,
//   0 <= 1-U(R) <= (4 R^{-1/2}-2 R^{-3/2}+3 R^{-5/2}) exp(-R).
// This is used to convert the finite-radius Newton root U(R;c_R)=1 into an
// enclosure of the actual asymptotic shooting parameter c_*.
static I vortexOneMinusUTailBound(const I& R0) {
  I sR = isqrt(R0);
  return (IV("4") / sR - two() / (R0 * sR) + IV("3") / (sqr(R0) * sR)) * iexp(-R0);
}

// -----------------------------------------------------------------------------
// CODED LEMMA A: vortex tail potential bound
//
// Statement.  For r >= R >= 4, |1-U(r)^2| is bounded by
//     8 exp(-r)/sqrt(r).
//
// Analytic input.  Lemma lem_1-U2 in the paper proves, for r>4,
//     U(r) > 1 - (4 r^{-1/2}-2 r^{-3/2}+3 r^{-5/2}) exp(-r)
// and U(r)<1.  Hence
//     0 <= 1-U(r)^2 <= 2(1-U(r))
//       <= (8 r^{-1/2}-4 r^{-3/2}+6 r^{-5/2}) exp(-r)
//       <= 8 exp(-r)/sqrt(r),
// since -4/r+6/r^2 <= 0 for r >= 3/2.  The CAPD endpoint check printed
// here is a consistency check against the propagated vortex enclosure at the
// chosen cutoff, but the global tail is supplied by the analytic lemma.
// -----------------------------------------------------------------------------

struct TailPotentialBound {
  I R;
  I A;
  I observedAtR;
  I barrierAtR;
  I entryMargin;
  I logKernelNorm;
  I firstMomentNorm;
};

static TailPotentialBound proveVortexTailPotentialBound(const Settings& S, const I& R0,
                                                        bool verbose = true) {
  I bridgeMargin = R0 - IV("4");
  I algebraMargin = R0 - IV("1.5");
  requireTrue(lowerAtLeast(bridgeMargin, "0"),
              "vortex tail bridge requires R>=4 from Lemma lem_1-U2");
  requireTrue(positive(algebraMargin), "vortex tail algebra requires R>3/2");

  Vec xv = initialVortex(S, S.c_cert);
  Vec xR = integrate(VF_VORTEX_FWD, xv, R0 - S.r0, S.taylor_order);
  I qR = absUpperAsInterval(one() - sqr(xR[1]));
  I barrier = S.vortex_tail_A * iexp(-R0) / isqrt(R0);
  I margin = barrier - qR;
  if (verbose) {
    std::cout << "vortex tail lemma at R=" << R0 << "\n";
    std::cout << "  imported analytic bridge: Lemma lem_1-U2 gives A=8 for all r>=R\n";
    std::cout << "  R-4 margin = " << bridgeMargin << ", R-3/2 margin = " << algebraMargin << "\n";
    std::cout << "  observed |1-U^2| <= " << qR << "\n";
    std::cout << "  proposed A exp(-R)/sqrt(R) = " << barrier << "\n";
    std::cout << "  endpoint consistency margin = " << margin << "\n";
  }
  requireTrue(positive(margin), "vortex tail barrier entry |1-U^2| <= A exp(-r)/sqrt(r)");

  // For R>=1, sqrt(s)<=s and log(s/R)<=s/R-1.  These closed forms bound
  // the Volterra kernels used by the threshold and outgoing coded lemmas:
  //   ∫_R^∞ s A e^-s/sqrt(s) log(s/R) ds <= A e^-R (R+2)/R,
  //   ∫_R^∞ s A e^-s/sqrt(s) ds          <= A e^-R (R+1).
  I logNorm = S.vortex_tail_A * iexp(-R0) * (R0 + two()) / R0;
  I firstMoment = S.vortex_tail_A * iexp(-R0) * (R0 + one());
  requireTrue(upperAtMost(logNorm, "1.012816572473333e-6"),
              "paper bound for the logarithmic vortex-tail norm at R=16");
  requireTrue(upperAtMost(firstMoment, "1.530478376181925e-5"),
              "paper bound for the first vortex-tail moment at R=16");
  requireTrue(upperAtMost(qR, "1.288693152516394e-7"),
              "paper consistency bound for |1-U(16)^2|");
  requireTrue(upperAtMost(barrier, "2.250703494385183e-7"),
              "paper value of the analytic vortex-tail barrier at R=16");
  return {R0, S.vortex_tail_A, qR, barrier, margin, logNorm, firstMoment};
}

static bool checkVortexTailBarrierAtFour(const Settings& S,
                                         const I& cBox,
                                         const std::string& boxLabel) {
  std::cout << "\n[Vortex tail-barrier entry inequalities at r=4]\n";

  // This block is the missing finite CAPD ingredient in Lemma lem_1-U2 of the
  // Part I paper.  The analytic comparison argument in the paper says:
  //
  //   If, at r0=4,
  //      A_lower(4) < a(4) < A_upper(4),
  //      U_lower(4) < U(4) < U_upper(4),
  //   then the explicit lower/upper barriers propagate for every r>4.
  //
  // The code below verifies exactly those four strict inequalities, using the
  // specified shooting box and CAPD propagation from the Frobenius start at
  // r=0.1 to r=4.  This is not a numerical plot check; the printed margins are
  // interval margins.  A positive left endpoint means the inequality is proved
  // for every trajectory in the supplied c-box.
  //
  // Dependency discipline:
  //   * In the shooting block this routine is called with the broad J0 box,
  //     after Proposition 4.1's endpoint dichotomy has placed the true c_* in
  //     J0 but before the finite-radius Newton bridge uses Lemma lem_1-U2.
  //     This breaks the circularity flagged in review.
  //   * Later profile-dependent checks use the refined c_cert only after the
  //     broad J0 entry check and the Newton bridge have both passed.
  I r = IV("4");
  std::cout << "c-box used for this entry check (" << boxLabel << ") = " << cBox << "\n";
  Vec x0 = initialVortex(S, cBox);
  Vec x4 = integrate(VF_VORTEX_FWD, x0, r - S.r0, S.taylor_order);
  I u = x4[1];
  I a = x4[2];

  I sqrtR = isqrt(r);
  I invSqrtR = one() / sqrtR;
  I expR = iexp(-r);
  I k0 = K0(r);
  I k1 = K1(r);

  I A_lower = one() - IV("4") * sqrtR * expR;
  I A_upper = one() - IV("1.4") * r * k1;
  I U_lower = one() - (IV("4") * invSqrtR - two() * invSqrtR / r) * expR;
  I U_upper = one() - IV("1.4") * k0 - sqr(k0);

  I marginALower = a - A_lower;
  I marginAUpper = A_upper - a;
  I marginULower = u - U_lower;
  I marginUUpper = U_upper - u;

  std::cout << "profile enclosure at r=4: U=" << u << ", a=" << a << "\n";
  std::cout << "  a(4)-A_lower(4) = " << marginALower << "\n";
  std::cout << "  A_upper(4)-a(4) = " << marginAUpper << "\n";
  std::cout << "  U(4)-U_lower(4) = " << marginULower << "\n";
  std::cout << "  U_upper(4)-U(4) = " << marginUUpper << "\n";

  bool ok = positive(marginALower) && positive(marginAUpper) &&
            positive(marginULower) && positive(marginUUpper);
  if (ok) {
    std::cout << "These four positive margins supply the CAPD entry check\n"
                 "for Lemma lem_1-U2 at r0=4.\n";
  }
  std::cout << pass(ok) << "\n";
  return ok;
}

// -----------------------------------------------------------------------------
// CODED LEMMA B: threshold infinity data
//
// Statement.  At R, the threshold solution normalized by f(r)/sqrt(r)->-1
// has an explicit interval enclosure for f(R),f'(R).
//
// Formula.  Write f=sqrt(r) w.  Then
//       (r w')' = r (U^2-1) w,
// and the Volterra map on r>=R has contraction norm bounded by the log-kernel
// tail norm in TailPotentialBound.  The code computes the contraction and the
// resulting bounds for w+1 and w'.
// -----------------------------------------------------------------------------

struct ThresholdTailBox {
  I f;
  I fp;
  I contraction;
  I wError;
  I wpError;
};

static ThresholdTailBox thresholdInfinityData(const Settings& S) {
  TailPotentialBound q = proveVortexTailPotentialBound(S, S.r_threshold_inf);
  I contraction = q.logKernelNorm;
  requireTrue(upperLess(contraction, "0.5"), "threshold Volterra contraction < 1/2");

  I denom = one() - contraction;
  I wErr = contraction / denom;
  I wpErr = (q.firstMomentNorm / q.R) / denom;
  I sr = isqrt(q.R);
  I f0 = -sr;
  I fp0 = -one() / (two() * sr);
  I fErr = sr * wErr;
  I fpErr = wErr / (two() * sr) + sr * wpErr;
  ThresholdTailBox box;
  box.f = f0 + symmetricError(fErr);
  box.fp = fp0 + symmetricError(fpErr);
  box.contraction = contraction;
  box.wError = wErr;
  box.wpError = wpErr;
  return box;
}

// -----------------------------------------------------------------------------
// CODED LEMMA C: internal decaying data at infinity
//
// Statement.  For mu in the supplied endpoint interval, the exact decaying
// internal-mode solution at R is enclosed by the K0/K1 model plus a Volterra
// error derived from the same verified potential tail.
// -----------------------------------------------------------------------------

struct InternalDecayBox {
  I psi;
  I psip;
  I alpha;
  I valueError;
  I derivError;
};

static InternalDecayBox internalDecayData(const Settings& S, const I& mu,
                                          bool verbose = true) {
  TailPotentialBound q = proveVortexTailPotentialBound(S, S.r_psi_tail, verbose);
  I alpha = isqrt(one() - mu);
  I x = alpha * q.R;
  I k0 = K0(x);
  I k1 = K1(x);

  // A conservative Volterra norm for the modified Bessel resolvent.  The
  // exponential gap is controlled by the same log-kernel norm, and the K0
  // Gaussian upper bound supplies the scale of the decaying branch.
  I contraction = q.logKernelNorm / alpha;
  requireTrue(upperLess(contraction, "0.25"), "internal-mode tail contraction < 1/4");
  I scale = K0_gaussian_upper(x);
  I valErr = scale * contraction / (one() - contraction);
  I derErr = alpha * valErr + scale * q.firstMomentNorm / (one() - contraction);

  InternalDecayBox box;
  box.psi = k0 + symmetricError(valErr);
  box.psip = -alpha * k1 + symmetricError(derErr);
  box.alpha = alpha;
  box.valueError = valErr;
  box.derivError = derErr;
  return box;
}

static I internalWronskianAtMu(
    const Settings& S, const I& mu,
    const InternalFrobeniusTailCertificate& broadPsiTail,
    bool verbose = true) {
  VortexState v0 = vortexSeries(S.r0, S.c_cert, S.series_order);
  PsiState p0 = psiSeriesUsingCertifiedTail(
      S.r0, S.c_cert, mu, S.series_order, broadPsiTail);
  Vec y0(7);
  y0[0] = S.r0; y0[1] = v0.u; y0[2] = v0.a; y0[3] = p0.psi; y0[4] = p0.psip; y0[5] = S.c_cert; y0[6] = mu;
  Vec yf = integrate(VF_INTERNAL_FWD, y0, S.r_match - S.r0, S.taylor_order);

  Vec xv = initialVortex(S, S.c_cert);
  Vec xvR = integrate(VF_VORTEX_FWD, xv, S.r_psi_tail - S.r0, S.taylor_order);
  InternalDecayBox dec = internalDecayData(S, mu, verbose);
  Vec yR(7);
  yR[0] = S.r_psi_tail; yR[1] = xvR[1]; yR[2] = xvR[2];
  yR[3] = dec.psi; yR[4] = dec.psip; yR[5] = S.c_cert; yR[6] = mu;
  Vec yb = integrate(VF_INTERNAL_BWD, yR, S.r_psi_tail - S.r_match, S.taylor_order);

  return yf[3] * yb[4] - yf[4] * yb[3];
}

static bool checkInternalEigenvalueBox(
    const Settings& S,
    const InternalFrobeniusTailCertificate& broadPsiTail) {
  std::cout << "\n[5a] Internal eigenvalue bracket coded lemma\n";
  I DL = internalWronskianAtMu(S, S.mu_eigen_left, broadPsiTail);
  I DR = internalWronskianAtMu(S, S.mu_eigen_right, broadPsiTail);
  I product = DL * DR;
  std::cout << "eigenvalue search interval = " << S.mu_eigen_box << "\n";
  std::cout << "Wronskian D(mu_left)  = " << DL << "\n";
  std::cout << "Wronskian D(mu_right) = " << DR << "\n";
  std::cout << "D(mu_left)*D(mu_right) = " << product << "\n";
  std::cout << "Uniqueness is supplied by the certified Seto eigenvalue-count block;\n"
               "the sign change brackets the unique internal eigenvalue.\n";
  bool printedWronskiansOk =
      subsetOfDecimalInterval(
          DL, "-4.074288757607958e-7", "-3.539855515662376e-7")
      && subsetOfDecimalInterval(
          DR, "7.716058532830248e-6", "7.768215008615830e-6");
  bool ok = negative(product) && printedWronskiansOk;
  std::cout << pass(ok) << "\n";
  return ok;
}

static bool oppositeStrictSigns(const I& a, const I& b) {
  return (negative(a) && positive(b)) || (positive(a) && negative(b));
}

static bool sameStrictSign(const I& a, const I& b) {
  return (negative(a) && negative(b)) || (positive(a) && positive(b));
}

static bool certifyAndInstallFgrEigenvalueBox(
    Settings& S,
    const InternalFrobeniusTailCertificate& broadPsiTail) {
  std::cout << "\n[5b] Certified internal eigenvalue box used by FGR\n";
  std::cout << "Input search box from the eigenvalue-counting block = "
            << S.mu_eigen_box << "\n";
    std::cout << "The initial nominal mu_box is not trusted here.  This routine\n"
               "constructs a new FGR interval by validated bisection of the\n"
               "Wronskian function, then overwrites S.mu_box with that\n"
               "certified interval before any K0 or FGR computation is run.\n";

  I nominalL(S.mu_box.leftBound());
  I nominalR(S.mu_box.rightBound());
  I DnominalL = internalWronskianAtMu(S, nominalL, broadPsiTail, false);
  I DnominalR = internalWronskianAtMu(S, nominalR, broadPsiTail, false);
  I nominalProduct = DnominalL * DnominalR;
  std::cout << "initial nominal mu box Wronskian diagnostic = " << S.mu_box << "\n";
  std::cout << "  Wronskian D(nominal left)  = " << DnominalL << "\n";
  std::cout << "  Wronskian D(nominal right) = " << DnominalR << "\n";
  std::cout << "  D(nominal left)*D(nominal right) = " << nominalProduct << "\n";
  if (oppositeStrictSigns(DnominalL, DnominalR)) {
    std::cout << "  diagnostic: nominal box has a strict Wronskian sign change.\n";
  } else if (sameStrictSign(DnominalL, DnominalR)) {
    std::cout << "  diagnostic: nominal box has no endpoint sign-change\n"
                 "  certificate for the current Wronskian convention.\n";
  } else {
    std::cout << "  diagnostic: nominal endpoint enclosure contains zero;\n"
                 "  refine around the nominal value before drawing conclusions.\n";
  }

  R lo = S.mu_eigen_left.leftBound();
  R hi = S.mu_eigen_right.rightBound();
  I Dlo = internalWronskianAtMu(S, S.mu_eigen_left, broadPsiTail, false);
  I Dhi = internalWronskianAtMu(S, S.mu_eigen_right, broadPsiTail, false);
  std::cout << "Wronskian D(search left)  = " << Dlo << "\n";
  std::cout << "Wronskian D(search right) = " << Dhi << "\n";
  requireTrue(oppositeStrictSigns(Dlo, Dhi),
              "FGR eigenvalue search box must have strict Wronskian sign change");

  const R targetWidth = IV("1e-15").leftBound();
  int iterations = 0;
  int ambiguousMidpoints = 0;

  // Validated bisection.  Each midpoint is evaluated as a singleton interval
  // through the same CAPD endpoint propagations used by the broad eigenvalue
  // check.  The sign-change bracket is never replaced by an uncertified
  // decimal guess: at every step the two stored endpoints have opposite
  // strict signs.  The Seto eigenvalue-count block supplies uniqueness, so this
  // existence bracket is also the interval for the unique internal eigenvalue.
  while ((I(hi) - I(lo)).rightBound() > targetWidth && iterations < 80) {
    R mid = (lo + hi) / R(2);
    I Dmid = internalWronskianAtMu(S, I(mid), broadPsiTail, false);
    if (sameStrictSign(Dlo, Dmid)) {
      lo = mid;
      Dlo = Dmid;
    } else if (oppositeStrictSigns(Dlo, Dmid)) {
      hi = mid;
      Dhi = Dmid;
    } else {
      // A midpoint enclosure containing zero is already useful information,
      // but it cannot replace the sign-change bracket.  Stop with the last
      // certified endpoints.  This is a deliberate diagnostic choice: the program
      // never invents a narrower decimal box after CAPD has lost sign
      // separation.  The downstream FGR block must either tolerate this
      // certified width or fail visibly.
      ++ambiguousMidpoints;
      break;
    }
    ++iterations;
  }

  I certified(lo, hi);
  I product = Dlo * Dhi;
  std::cout << "validated bisection steps = " << iterations << "\n";
  std::cout << "ambiguous midpoint stops = " << ambiguousMidpoints << "\n";
  std::cout << "certified FGR mu box = " << certified << "\n";
  printWidth("certified FGR mu box", certified);
  std::cout << "Wronskian D(certified left)  = " << Dlo << "\n";
  std::cout << "Wronskian D(certified right) = " << Dhi << "\n";
  std::cout << "D(left)*D(right) = " << product << "\n";

  // Tie the computation directly to the exact decimal interval printed in the
  // paper.  Each endpoint is passed as a directed interval containing the exact
  // decimal, and all later computations use the outward enclosure of that same
  // printed interval.
  I DpaperLeft = internalWronskianAtMu(
      S, S.mu_fgr_left, broadPsiTail, false);
  I DpaperRight = internalWronskianAtMu(
      S, S.mu_fgr_right, broadPsiTail, false);
  I paperProduct = DpaperLeft * DpaperRight;
  std::cout << "exact printed FGR mu box = " << S.mu_fgr_box << "\n";
  std::cout << "Wronskian D(exact printed left)  = " << DpaperLeft << "\n";
  std::cout << "Wronskian D(exact printed right) = " << DpaperRight << "\n";
  std::cout << "D(exact printed left)*D(exact printed right) = "
            << paperProduct << "\n";

  // The width threshold is only a guardrail.  The sign changes at both the
  // bisection endpoints and the exact printed endpoints are the existence
  // certificates; the Seto bound already proved uniqueness.
  bool ok = negative(product)
         && oppositeStrictSigns(DpaperLeft, DpaperRight)
         && iterations == 5
         && ambiguousMidpoints == 1
         && subsetOfDecimalInterval(
              DpaperLeft,
              "-1.536081253326721e-7", "-1.002050065350549e-7")
         && subsetOfDecimalInterval(
              DpaperRight,
              "1.002150911811321e-7", "1.535780040260618e-7")
         && upperLess(intervalWidth(S.mu_fgr_box), "2e-6");
  if (ok) {
    S.mu_box = S.mu_fgr_box;
    std::cout << "Installed the outward enclosure of the exact printed FGR mu box "
                 "into Settings.mu_box.\n";
  } else {
    std::cout << "The FGR eigenvalue box was not narrow enough for the\n"
                 "downstream FGR certificate.  No decimal fallback is used.\n";
  }
  std::cout << pass(ok) << "\n";
  return ok;
}

// -----------------------------------------------------------------------------
// CODED/PAPER LEMMA D: outgoing asymptotic start for FGR Weyl solutions
//
// Statement.  The exact outgoing solution at R is enclosed by the differentiated
// finite outgoing ansatz plus a rigorously bounded residual/Volterra remainder.
//
// The proof is now a direct a posteriori argument, not a last-term rule for an
// asymptotic series.  Write
//
//   y(r)=k^{-1/2} exp(i*k*r) z(r).
//
// Then z solves
//
//   z'' + 2 i k z' + (alpha/r^2 + V(r)) z = 0,
//
// where |V(r)| is bounded by the vortex tail.  The truncated polynomial z_N
// has an explicitly computable residual R_N.  The error e=z-z_N satisfies the
// free outgoing Volterra equation with kernel bounded by 1/k for values and by
// 1 for derivatives.  The code below bounds the residual forcing R_N, the
// true-coefficient forcing V z_N, and the contraction
//
//   (1/k) int_R^infty (|alpha|/s^2+|V(s)|) ds.
//
  // This is the r >= 16 analogue of the Frobenius series-tail check: an infinite
// tail is reduced to finite coefficient data plus closed-form integral bounds.
// -----------------------------------------------------------------------------

static std::vector<CI> outgoingCoeffs(int which, int N, bool verbose) {
  // Generate the outgoing coefficients from the limiting equation actually
  // solved by vfOutgoing.  This replaces a previous table, which was copied
  // from a notebook convention and did not satisfy the recurrence for the ODE
  // in this file.
  //
  // At infinity, with k^2=4*mu-1, the two channel equations in vfOutgoing are:
  //
  //   channel 1:  y'' + (k^2 - 3/(4r^2)) y = 0,
  //   channel 2:  y'' + (k^2 + 1/(4r^2)) y = 0.
  //
  // In the dimensionless variable z=k*r these have the common form
  //
  //   y''_z + (1 + alpha/z^2) y = 0,
  //
  // with alpha=-3/4 in channel 1 and alpha=1/4 in channel 2.  For the outgoing
  // ansatz
  //
  //   y(r) = k^{-1/2} exp(i*k*r) sum_{n>=0} c_n (k*r)^{-n},
  //
  // the formal recurrence is
  //
  //   -2*i*n*c_n + ((n-1)*n + alpha)c_{n-1} = 0,  n>=1.
  //
  // Thus
  //
  //   c_n = -i * (((n-1)*n + alpha)/(2*n)) * c_{n-1}.
  //
  // The first terms are consequently
  //   channel 1: 1, 3i/8, 15/128, -105i/1024, ...
  //   channel 2: 1, -i/8, -9/128, 75i/1024, ...
  //
  // These are exactly the coefficients of the limiting Bessel/Hankel model
  // associated with the ODE solved below.  The residual check after generation
  // is intentionally redundant: it makes a transcription error visible in the
  // certificate transcript.
  requireTrue(N >= 2, "outgoing coefficient generation needs at least two terms");
  I alpha = (which == 1 ? -IV("0.75") : IV("0.25"));
  std::vector<CI> c;
  c.reserve(static_cast<size_t>(N + 1));
  c.push_back(CI(one()));
  for (int n = 1; n <= N; ++n) {
    I A = IV(std::to_string(n * (n - 1)).c_str()) + alpha;
    I den = IV(std::to_string(2 * n).c_str());
    c.push_back(mulI(c.back()) * (-(A / den)));
  }

  I maxResidual = zero();
  for (int n = 1; n <= N; ++n) {
    I A = IV(std::to_string(n * (n - 1)).c_str()) + alpha;
    I den = IV(std::to_string(2 * n).c_str());
    CI residual = mulI(c[n]) * (-den) + c[n - 1] * A;
    requireTrue(!excludesZero(residual.re) && !excludesZero(residual.im),
                "outgoing coefficient recurrence residual encloses zero");
    I resAbs = absUpper(residual);
    if (resAbs.rightBound() > maxResidual.rightBound()) maxResidual = resAbs;
  }

  if (verbose) {
    std::cout << "outgoing channel " << which << " coefficients generated from limiting recurrence\n";
    std::cout << "  alpha = " << alpha << ", order = " << N << "\n";
    std::cout << "  c0 = " << c[0] << "\n";
    std::cout << "  c1 = " << c[1] << "\n";
    std::cout << "  c2 = " << c[2] << "\n";
    std::cout << "  recurrence residual max enclosure = " << maxResidual << "\n";
  }
  return c;
}

struct OutgoingStart {
  CI y;
  CI yp;
  I valueError;
  I derivError;
  I contraction;
  I residualIntegral;
  I potentialIntegral;
  I zError;
  I zPrimeError;
  I matrixTailBound;
};

static OutgoingStart outgoingDataAtR(const Settings& S, int which, const I& k) {
  I R0 = S.r_fgr;
  I kr = k * R0;
  const int N = S.asym_order;
  std::vector<CI> coeffs = outgoingCoeffs(which, S.asym_order, true);
  CI s, dsdr;
  for (int n = 0; n < static_cast<int>(coeffs.size()); ++n) {
    I power = ipow(kr, static_cast<unsigned>(n));
    s = s + coeffs[n] / power;
    if (n > 0) {
      I factor = -IV(std::to_string(n).c_str()) * k / ipow(kr, static_cast<unsigned>(n + 1));
      dsdr = dsdr + coeffs[n] * factor;
    }
  }
  CI phase = cis(k * R0) * (one() / isqrt(k));
  CI y = phase * s;
  CI yp = phase * (mulI(s) * k + dsdr);

  TailPotentialBound q = proveVortexTailPotentialBound(S, R0);

  I alpha = (which == 1 ? -IV("0.75") : IV("0.25"));
  I alphaAbs = upperSingleton(absUpperAsInterval(alpha));
  I zSup = zero();
  I zPrimeSup = zero();
  for (int n = 0; n <= N; ++n) {
    I term = upperSingleton(absUpper(coeffs[n]) / ipow(kr, static_cast<unsigned>(n)));
    zSup += term;
    if (n > 0) zPrimeSup += IV(std::to_string(n).c_str()) * term / R0;
  }

  // Exact residual of z_N for the limiting equation.  The recurrence cancels
  // all lower powers; only ((N)(N+1)+alpha)c_N (kr)^(-N-2) remains after
  // applying z''+2ikz'+alpha r^{-2}z.
  I residualCoeff =
      (IV(std::to_string(N * (N + 1)).c_str()) + alpha) * absUpper(coeffs.back());
  residualCoeff = upperSingleton(absUpperAsInterval(residualCoeff));
  I residualIntegral =
      residualCoeff * k / (IV(std::to_string(N + 1).c_str()) *
                           ipow(kr, static_cast<unsigned>(N + 1)));

  // For |V(s)| <= A exp(-s)/sqrt(s), s>=R,
  // int_R^infty |V(s)| ds <= A exp(-R)/sqrt(R).
  I potentialIntegral = q.A * iexp(-R0) / isqrt(R0);
  I wIntegral = alphaAbs / R0 + potentialIntegral;
  I contraction = wIntegral / k;
  requireTrue(upperLess(contraction, "0.25"), "outgoing residual Volterra contraction < 1/4");

  I forcingValue = (residualIntegral + potentialIntegral * zSup) / k;
  I zErr = forcingValue / (one() - contraction);
  I zPrimeErr = residualIntegral + potentialIntegral * zSup + wIntegral * zErr;

  I valErr = zErr / isqrt(k);
  I derErr = (k * zErr + zPrimeErr) / isqrt(k);

  // Uniform upper bound for the matrix Weyl functions used in the FGR source
  // bound.  Lemma 4.6 gives |1-a| <= (33/8) exp(-R) sqrt(R) for r>=R,
  // hence |b|=|1-a|/r <= (33/8) exp(-R)/sqrt(R).  Then
  // | -y' - y/(2r) - b y | and |U y| are bounded by the following quantity.
  I zTotal = zSup + zErr;
  I zPrimeTotal = zPrimeSup + zPrimeErr;
  I yTail = zTotal / isqrt(k);
  I ypTail = (k * zTotal + zPrimeTotal) / isqrt(k);
  I bTail = (IV("33") / IV("8")) * iexp(-R0) / isqrt(R0);
  I matrixTail = maxUpperInterval(ypTail + (one()/(two()*R0) + bTail) * yTail, yTail);
  requireTrue(upperLess(matrixTail, "3"), "outgoing matrix Weyl tail bound < 3");
  requireTrue(upperLess(matrixTail, "1.253"), "outgoing matrix Weyl tail bound < 1.253");

  std::cout << "outgoing channel " << which << " residual/Volterra enclosure for r >= 16\n";
  std::cout << "  z_N sup on [R,inf] <= " << zSup
            << ", z_N' sup <= " << zPrimeSup << "\n";
  std::cout << "  residual integral <= " << residualIntegral
            << ", potential integral <= " << potentialIntegral << "\n";
  std::cout << "  Volterra contraction = " << contraction << "\n";
  std::cout << "  z error <= " << zErr << ", z' error <= " << zPrimeErr << "\n";
  std::cout << "  outgoing matrix upper bound for r >= 16 <= " << matrixTail << "\n";

  OutgoingStart out;
  out.y = CI(y.re + symmetricError(valErr), y.im + symmetricError(valErr));
  out.yp = CI(yp.re + symmetricError(derErr), yp.im + symmetricError(derErr));
  out.valueError = valErr;
  out.derivError = derErr;
  out.contraction = contraction;
  out.residualIntegral = residualIntegral;
  out.potentialIntegral = potentialIntegral;
  out.zError = zErr;
  out.zPrimeError = zPrimeErr;
  out.matrixTailBound = matrixTail;
  return out;
}

// -----------------------------------------------------------------------------
// Check 1: vortex bracket coded lemma
// -----------------------------------------------------------------------------

struct VortexBadBehaviorCert {
  bool ok;
  I radius;
  I margin;
};

static VortexBadBehaviorCert certifyVortexBadBehaviorAtRadius(
    const Settings& S,
    const I& c,
    bool wantABad,
    const I& radius) {
  // Proposition 4.1 in the paper is a shooting dichotomy/uniqueness theorem:
  // the desired vortex parameter c_* is the unique separator between the
  // "a-bad" solutions and the "U-bad" solutions.  The computer should not
  // re-prove that topological ODE theorem here; what it can certify is the
  // finite-radius hypothesis needed to place the two endpoints of J0 on the
  // opposite sides of that dichotomy.
  //
  // For the lower endpoint we try to prove a-badness by finding a finite radius
  // where a(r;c_lower)>2.  For the upper endpoint we try to prove U-badness by
  // finding a finite radius where U(r;c_upper)>2.  The threshold 2 is chosen
  // deliberately away from the separating value 1, so a positive interval
  // margin is a robust, directly checkable certificate of bad behavior.
  //
  // The two radii are the exact decimals stated in the paper.  Passing them as
  // directed intervals ensures that the CAPD propagation contains those exact
  // values and does not depend on an accumulated floating-point step size.
  Vec x = initialVortex(S, c);
  x = integrate(VF_VORTEX_FWD, x, radius - S.r0, S.taylor_order);
  I margin = (wantABad ? x[2] : x[1]) - IV("2");
  return {positive(margin), x[0], margin};
}

static bool checkVortexBracket(const Settings& S) {
  std::cout << "\n[1] Vortex shooting interval coded lemma\n";
  std::cout << "J0 = " << S.c_box << "\n";
  bool midpointInside = S.c_mid.leftBound() >= S.c_box.leftBound() &&
                        S.c_mid.rightBound() <= S.c_box.rightBound();
  std::cout << "midpoint value inside J0: " << pass(midpointInside) << "\n";

  Vec xl0 = initialVortex(S, S.c_lower);
  Vec xr0 = initialVortex(S, S.c_upper);
  Vec xl = integrate(VF_VORTEX_FWD, xl0, S.r_shoot - S.r0, S.taylor_order);
  Vec xr = integrate(VF_VORTEX_FWD, xr0, S.r_shoot - S.r0, S.taylor_order);
  I lowerEvent = xl[1] - one();
  I upperEvent = xr[1] - one();
  std::cout << "shooting event U(R)-1 at lower endpoint = " << lowerEvent << "\n";
  std::cout << "shooting event U(R)-1 at upper endpoint = " << upperEvent << "\n";

  VortexBadBehaviorCert lowerBad = certifyVortexBadBehaviorAtRadius(
      S, S.c_lower, true, IV("27.75"));
  VortexBadBehaviorCert upperBad = certifyVortexBadBehaviorAtRadius(
      S, S.c_upper, false, IV("31.6"));
  std::cout << "Proposition 4.1 shooting-dichotomy endpoint check\n";
  std::cout << "  lower endpoint a-bad certificate: a(r;c_lower)-2 = "
            << lowerBad.margin << " at r=" << lowerBad.radius
            << " -> " << pass(lowerBad.ok) << "\n";
  std::cout << "  upper endpoint U-bad certificate: U(r;c_upper)-2 = "
            << upperBad.margin << " at r=" << upperBad.radius
            << " -> " << pass(upperBad.ok) << "\n";
  bool endpointRadiiOk = containsInterval(lowerBad.radius, IV("27.75"))
                      && containsInterval(upperBad.radius, IV("31.6"));
  bool printedMarginsOk = subsetOfDecimalInterval(
                              lowerBad.margin,
                              "0.01048373219531", "0.01048373219532")
                       && subsetOfDecimalInterval(
                              upperBad.margin,
                              "0.00213149536392", "0.00213149536394");
  std::cout << "  propagated radius intervals contain the exact stated radii: "
            << pass(endpointRadiiOk) << "\n";
  std::cout << "  margins lie in the intervals printed in Proposition prop:J0slope: "
            << pass(printedMarginsOk) << "\n";
  std::cout << "  With Proposition 4.1's analytic dichotomy and uniqueness, these\n"
               "  opposite finite-radius bad behaviors certify c_* in J0.\n";
  bool broadOk = midpointInside && lowerBad.ok && upperBad.ok
              && endpointRadiiOk && printedMarginsOk;

  // This is the key acyclic ordering.  We first know c_* lies in J0 by the
  // shooting dichotomy.  We then verify the r=4 barrier entry inequalities on
  // the whole J0 box.  Only after that do we use Lemma lem_1-U2 to convert the
  // finite-radius Newton root U(20;c_R)=1 into a box for the true asymptotic
  // shooting parameter c_*.
  bool broadTailEntryOk = false;
  if (broadOk) {
    broadTailEntryOk = checkVortexTailBarrierAtFour(S, S.c_box, "broad J0, before Newton bridge");
  } else {
    std::cout << "\n[Vortex tail-barrier entry inequalities at r=4]\n"
                 "Skipped: the broad shooting dichotomy did not yet certify c_* in J0.\n";
  }

  I cMid = midpointInterval(S.c_cert);
  Vec xMid0 = initialVortex(S, cMid);
  Vec xMid = integrate(VF_VORTEX_FWD, xMid0, S.r_shoot - S.r0, S.taylor_order);
  I Fmid = xMid[1] - one();
  Vec xVar0 = initialVortexVar(S, S.c_cert);
  Vec xVar = integrate(VF_VORTEX_VAR_FWD, xVar0, S.r_shoot - S.r0, S.taylor_order);
  I dF = xVar[3];
  I newton = cMid - Fmid / dF;
  bool printedNewtonDataOk =
      subsetOfDecimalInterval(
          Fmid,
          "-4.08767126613755487333512480664e-9",
          "-4.08767126613755487333512480663e-9")
      && subsetOfDecimalInterval(
          dF,
          "2.85328651498644319085546375648e7",
          "2.85328868245274473217494026266e7")
      && subsetOfDecimalInterval(
          newton,
          "0.603287854581677893261748847078",
          "0.603287854581677893261857674219");
  bool newtonOk = excludesZero(dF)
               && strictSubsetOfDecimalInterval(
                    newton, "0.6032878545816699", "0.6032878545816856")
               && printedNewtonDataOk;
  std::cout << "finite-R interval Newton radius = " << S.r_shoot << "\n";
  std::cout << "  F(c_mid)=U(R;c_mid)-1 = " << Fmid << "\n";
  std::cout << "  dF/dc on refined box = " << dF << "\n";
  std::cout << "  Newton image = " << newton << "\n";
  printWidth("  Newton image", newton);
  std::cout << "  Newton image lies strictly inside the exact printed refined box: "
            << pass(newtonOk) << "\n";
  std::cout << "  Newton data lie in the intervals printed in Lemma lem:Newton: "
            << pass(printedNewtonDataOk) << "\n";

  Vec xVarBroad0 = initialVortexVar(S, S.c_box);
  Vec xVarBroad = integrate(VF_VORTEX_VAR_FWD, xVarBroad0, S.r_shoot - S.r0, S.taylor_order);
  I dFBroad = xVarBroad[3];
  I trueVortexTailAtR = upperSingleton(vortexOneMinusUTailBound(S.r_shoot));
  I cBridgeRadius = upperSingleton(trueVortexTailAtR / dFBroad);
  I trueCBox = widen(newton, cBridgeRadius);
  bool printedBridgeDataOk =
      lowerAtLeast(
          dFBroad, "2.85319410965073633711311481846e7")
      && upperLess(
          trueVortexTailAtR,
          "1.80091970832747686279582580145e-9")
      && upperLess(cBridgeRadius, "6.31194247259935e-17")
      && subsetOfDecimalInterval(
          trueCBox,
          "0.603287854581677830142324121084",
          "0.603287854581677956381282400213");
  bool bridgeOk = broadOk && broadTailEntryOk && newtonOk &&
                  positive(dFBroad)
                  && strictSubsetOfDecimalInterval(
                       trueCBox, "0.6032878545816699", "0.6032878545816856")
                  && printedBridgeDataOk;
  std::cout << "finite-R/asymptotic shooting bridge\n";
  std::cout << "  dF/dc on broad J0 = " << dFBroad << "\n";
  std::cout << "  true-vortex bound 1-U(R;c_*) <= " << trueVortexTailAtR << "\n";
  std::cout << "  c_* distance from finite-R root <= " << cBridgeRadius << "\n";
  std::cout << "  resulting c_* enclosure = " << trueCBox << "\n";
  printWidth("  resulting c_* enclosure", trueCBox);
  std::cout << "  true asymptotic c_* lies strictly inside the exact printed refined box: "
            << pass(bridgeOk) << "\n";
  std::cout << "  bridge data satisfy the bounds printed in Lemma lem:Newton: "
            << pass(printedBridgeDataOk) << "\n";

  Vec xcl0 = initialVortex(S, I(S.c_cert.leftBound()));
  Vec xcr0 = initialVortex(S, I(S.c_cert.rightBound()));
  Vec xcl = integrate(VF_VORTEX_FWD, xcl0, S.r_refined_shoot - S.r0, S.taylor_order);
  Vec xcr = integrate(VF_VORTEX_FWD, xcr0, S.r_refined_shoot - S.r0, S.taylor_order);
  I certLowerEvent = xcl[1] - one();
  I certUpperEvent = xcr[1] - one();
  std::cout << "refined c_* box = " << S.c_cert << "\n";
  std::cout << "refined shooting radius = " << S.r_refined_shoot << "\n";
  std::cout << "refined lower event U(R)-1 = " << certLowerEvent << "\n";
  std::cout << "refined upper event U(R)-1 = " << certUpperEvent << "\n";
  bool refinedSignDiagnostic = negative(certLowerEvent) && positive(certUpperEvent);
  std::cout << "long-radius endpoint-sign diagnostic: "
            << pass(refinedSignDiagnostic) << "\n";
  if (newtonOk) {
    std::cout << "The active refined-box certificate is the finite-radius interval\n"
                 "Newton inclusion above.  The longer-radius endpoint signs are\n"
                 "printed only as a conditioning diagnostic; the long propagation\n"
                 "to R=30 is deliberately not used as the refined-box proof.\n";
  } else {
    std::cout << "The broad paper interval is bracketed, but the refined box is\n"
                 "not yet certified by the finite-radius interval Newton test.\n"
                 "Later checks using the refined box are therefore conditional\n"
                 "diagnostics until this coded lemma is strengthened.\n";
  }
  printVortexOriginTailCertificates(S);
  bool ok = broadOk && broadTailEntryOk && newtonOk && bridgeOk;
  std::cout << pass(ok) << "\n";
  return ok;
}

// -----------------------------------------------------------------------------
// Check 2: H1 positivity V1-1 > 0 on [sqrt(6),4]
// -----------------------------------------------------------------------------

static bool checkH1(const Settings& S) {
  std::cout << "\n[2] H1 positivity on [sqrt(6),4]\n";
  Vec x = initialVortex(S, S.c_cert);
  x = integrate(VF_VORTEX_FWD, x, S.r_h1_a - S.r0, S.taylor_order);
  I h = (S.r_h1_b - S.r_h1_a) / IV(std::to_string(S.h1_cells).c_str());
  I minG = IV("1e9");
  for (int j = 0; j < S.h1_cells; ++j) {
    StepResult st = integrateTube(VF_VORTEX_FWD, x, h, S.taylor_order);
    I r = st.tube[0], u = st.tube[1], a = st.tube[2];
    I G = sqr(two() - a) / sqr(r) - (one() - sqr(u)) / two();
    if (G.leftBound() < minG.leftBound()) minG = G;
    x = st.endpoint;
  }
  std::cout << "enclosed min lower bound for V1-1: " << minG << "\n";
  printWidth("V1-1 minimum enclosure", minG);
  std::cout << "paper-level target checked here: V1-1 > 0.0548 on this interval\n";
  bool ok = lowerGreater(minG, "0.0548");
  std::cout << pass(ok) << "\n";
  return ok;
}

// -----------------------------------------------------------------------------
// Check 3: Seto logarithmic-kernel eigenvalue-count bound
// -----------------------------------------------------------------------------

struct CellBound {
  I a, b;      // r interval endpoints as intervals/singletons
  I fUpper;    // pointwise upper bound for f(r)=(1-U^2)r on this cell
};

static I intLog(const I& a, const I& b) {
  // integral_a^b log(x) dx = b log b - b - a log a + a,
  // with the limiting value a log(a)-a -> 0 at a=0.
  I upper = b * ilog(b) - b;
  if (a.rightBound() == R(0)) return upper;
  return upper - a * ilog(a) + a;
}

static bool sameCell(const I& a, const I& b, const I& c, const I& d) {
  return a.leftBound() == c.leftBound() && a.rightBound() == c.rightBound() &&
         b.leftBound() == d.leftBound() && b.rightBound() == d.rightBound();
}

static I diagonalKernelUpper(const I& a, const I& b) {
  // Exact integral over one diagonal cell:
  //   ∫_a^b ∫_a^b |log(r/s)| dr ds
  // = 2 ∫_a^b [s-a-a log(s/a)] ds.
  // This removes the large artificial sup-norm loss near r=s.
  I width = b - a;
  I logRatioIntegral = intLog(a, b) - width * ilog(a);
  I exact = (sqr(b) - sqr(a)) - two() * a * width - two() * a * logRatioIntegral;
  return nonnegativePartUpper(exact);
}

static I kernelRectUpper(const I& a, const I& b, const I& c, const I& d) {
  I widthR = b - a;
  I widthS = d - c;
  if (sameCell(a, b, c, d)) {
    return diagonalKernelUpper(a, b);
  }
  if (b.rightBound() <= c.leftBound()) {
    // s >= r: |log(r/s)| = log s - log r
    return widthR * intLog(c, d) - widthS * intLog(a, b);
  }
  if (d.rightBound() <= a.leftBound()) {
    // r >= s
    return widthS * intLog(a, b) - widthR * intLog(c, d);
  }
  // Any other overlap should only occur from outward endpoint rounding.  Use a
  // valid sup bound there; true diagonal cells are handled exactly above.
  I lo(a.leftBound() < c.leftBound() ? a.leftBound() : c.leftBound());
  I hi(b.rightBound() > d.rightBound() ? b.rightBound() : d.rightBound());
  I sup = ilog(hi / lo);
  return widthR * widthS * sup;
}

static I intAbsLog(const I& a, const I& b) {
  // ∫_a^b |log(x)| dx, split at x=1 if needed.
  if (b.rightBound() <= R(1)) return -intLog(a, b);
  if (a.leftBound() >= R(1)) return intLog(a, b);
  return -intLog(a, one()) + intLog(one(), b);
}

static I ltTailComplementUpper(const Settings& S, const I& R0,
                               const I& compactMass,
                               const I& compactLogMoment) {
  // For r>=R0>=15, the vortex tail lower bound gives roughly
  // 1-U^2 <= 8 exp(-r)/sqrt(r), hence f(r)=(1-U^2)r <= 8 r exp(-r).
  // The cross term uses compact mass and log moment computed from the same
  // rectangle cells, not a global flux/log substitute.
  I R = R0;
  I A = S.vortex_tail_A;
  I mt0 = A * (R + one()) * iexp(-R);                              // ∫ A r e^-r
  I mtlog = A * (sqr(R) + two()*R + two()) * iexp(-R);              // log r <= r
  I cross = two() * (compactMass * mtlog + compactLogMoment * mt0);
  I tailtail = two() * mt0 * mtlog;
  return cross + tailtail;
}

static I ltMissingOriginUpper(const I& eps, const I& R) {
  // Bound the contribution of (0,eps) x (0,R) and its transpose using
  // f(r)<=r, total compact mass <=2, and a crude compact log moment.
  I m0eps = sqr(eps) / two();
  I mlogeps = sqr(eps) * (-ilog(eps)) / two() + sqr(eps) / IV("4");
  I mc0 = IV("2");
  I mclog = sqr(R) * ilog(R) / two() + IV("1");
  return two() * (m0eps * mclog + mlogeps * mc0);
}

static bool checkLT(const Settings& S) {
  std::cout << "\n[3] Seto logarithmic-kernel upper bound\n";
  std::vector<CellBound> cells;
  cells.reserve(S.lt_cells + S.origin_cells);

  // Origin cells [eps,0.1]: use f(r) <= r.  The tiny omitted
  // (0,eps) contribution is added explicitly below.
  I eps0 = IV("1e-30");
  I h0 = (S.r0 - eps0) / IV(std::to_string(S.origin_cells).c_str());
  for (int i = 0; i < S.origin_cells; ++i) {
    I a = eps0 + h0 * IV(std::to_string(i).c_str());
    I b = eps0 + h0 * IV(std::to_string(i + 1).c_str());
    cells.push_back({a, b, b});
  }

  // CAPD cells [0.1,15].
  Vec x = initialVortex(S, S.c_cert);
  I hc = (S.r_lt - S.r0) / IV(std::to_string(S.lt_cells).c_str());
  for (int j = 0; j < S.lt_cells; ++j) {
    I a = S.r0 + hc * IV(std::to_string(j).c_str());
    I b = S.r0 + hc * IV(std::to_string(j + 1).c_str());
    StepResult st = integrateTube(VF_VORTEX_FWD, x, hc, S.taylor_order);
    I r = st.tube[0], u = st.tube[1];
    I f = (one() - sqr(u)) * r;
    cells.push_back({a, b, nonnegativePartUpper(f)});
    x = st.endpoint;
  }

  I sum = zero();
  I compactMass = zero();
  I compactLogMoment = zero();
  const int n = static_cast<int>(cells.size());
  for (const CellBound& cell : cells) {
    I width = cell.b - cell.a;
    compactMass += cell.fUpper * width;
    compactLogMoment += cell.fUpper * intAbsLog(cell.a, cell.b);
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      I K = kernelRectUpper(cells[i].a, cells[i].b, cells[j].a, cells[j].b);
      sum += cells[i].fUpper * cells[j].fUpper * K;
    }
  }
  I complement = ltTailComplementUpper(S, S.r_lt, compactMass, compactLogMoment);
  I missingOrigin = ltMissingOriginUpper(IV("1e-30"), S.r_lt);
  I total = sum + complement + missingOrigin;
  std::cout << "compact mass upper = " << compactMass << "\n";
  std::cout << "compact |log r| moment upper = " << compactLogMoment << "\n";
  std::cout << "compact upper = " << sum << "\n";
  std::cout << "tail/complement upper = " << complement << "\n";
  std::cout << "missing-origin upper = " << missingOrigin << "\n";
  std::cout << "total upper = " << total << "\n";
  bool ok = upperLess(compactMass, "2.008204")
         && upperLess(compactLogMoment, "1.476031")
         && upperLess(sum, "3.17965")
         && upperLess(complement, "0.002642")
         && upperLess(missingOrigin, "4.449e-58")
         && upperLess(total, "3.182292");
  std::cout << pass(ok) << "\n";
  return ok;
}

// -----------------------------------------------------------------------------
// Check 4: H2 threshold Wronskian
// -----------------------------------------------------------------------------

static bool checkThreshold(
    const Settings& S,
    const PhiFrobeniusTailCertificate& thresholdPhi2Tail) {
  std::cout << "\n[4] H2 threshold Wronskian\n";
  VortexState v0 = vortexSeries(S.r0, S.c_cert, S.series_order);
  PhiState th0 = phiEnergySeriesUsingCertifiedTail(
      S.r0, S.c_cert, one(), 2, S.series_order, thresholdPhi2Tail);
  Vec y0(6);
  y0[0] = S.r0; y0[1] = v0.u; y0[2] = v0.a; y0[3] = th0.phi; y0[4] = th0.phip; y0[5] = S.c_cert;
  Vec yf = integrate(VF_THRESHOLD_FWD, y0, S.r_match - S.r0, S.taylor_order);

  // Use the CAPD-propagated vortex box at the coded threshold-tail radius.
  Vec xv = initialVortex(S, S.c_cert);
  Vec xvR = integrate(VF_VORTEX_FWD, xv, S.r_threshold_inf - S.r0, S.taylor_order);
  I Rr = S.r_threshold_inf;
  ThresholdTailBox tail = thresholdInfinityData(S);
  std::cout << "threshold tail contraction = " << tail.contraction << "\n";
  std::cout << "threshold w error = " << tail.wError << ", w' error = " << tail.wpError << "\n";
  requireTrue(upperLess(tail.contraction, "1.012817e-6"),
              "threshold Volterra norm < 1.012817e-6");
  requireTrue(upperLess(tail.wError, "1.012818e-6"),
              "threshold value error < 1.012818e-6");
  requireTrue(upperLess(tail.wpError, "9.56550e-7"),
              "threshold derivative error < 9.56550e-7");

  Vec yb0(6);
  yb0[0] = Rr; yb0[1] = xvR[1]; yb0[2] = xvR[2]; yb0[3] = tail.f; yb0[4] = tail.fp; yb0[5] = S.c_cert;
  Vec yb = integrate(VF_THRESHOLD_BWD, yb0, Rr - S.r_match, S.taylor_order);
  I W = yf[3] * yb[4] - yf[4] * yb[3];
  std::cout << "W interval = " << W << "\n";
  bool ok = strictSubsetOfDecimalInterval(
      W, "-0.8003267", "-0.8001801");
  std::cout << pass(ok) << "\n";
  return ok;
}

// -----------------------------------------------------------------------------
// Check 5: Internal mode K0 comparison
// -----------------------------------------------------------------------------

static bool checkInternalK0(
    const Settings& S,
    const InternalFrobeniusTailCertificate& psiTail) {
  std::cout << "\n[5] Internal mode: one-point K0 comparison\n";
  I r8 = S.r_k0_start;
  VortexState v0 = vortexSeries(S.r0, S.c_cert, S.series_order);
  PsiState p0 = psiSeriesUsingCertifiedTail(
      S.r0, S.c_cert, S.mu_box, S.series_order, psiTail);
  Vec y0(7);
  y0[0] = S.r0; y0[1] = v0.u; y0[2] = v0.a; y0[3] = p0.psi; y0[4] = p0.psip; y0[5] = S.c_cert; y0[6] = S.mu_box;
  Vec y8 = integrate(VF_INTERNAL_FWD, y0, r8 - S.r0, S.taylor_order);

  I kappa = isqrt(one() - S.delta0 - S.mu_box);
  I x = kappa * r8;
  I K0x = K0(x);
  I K1x = K1(x);
  I C = S.k0_comparison_C;
  I valGap = C * K0x - y8[3];
  I derivGap = y8[4] - C * (-kappa * K1x);
  I u2Gap = sqr(y8[1]) - (one() - S.delta0);
  I xTail = kappa * S.r_fgr;
  I k1OverK0Majorant = one() + one() / (two() * xTail);
  I derivativeK0Factor = (one() - S.mu_box) * k1OverK0Majorant / kappa;
  std::cout << "comparison start radius = " << r8 << "\n";
  std::cout << "kappa = " << kappa << "\n";
  std::cout << "value gap 1.4*K0(kappa*r)-psi(r) = " << valGap << "\n";
  std::cout << "derivative gap psi'(r)-(1.4*K0)' = " << derivGap << "\n";
  std::cout << "U(r)^2-(1-delta0) = " << u2Gap << "\n";
  std::cout << "large-r derivative factor ((1-mu)/kappa)*(1+1/(2*kappa*16)) = "
            << derivativeK0Factor << "\n";
  std::cout << "  paper target < 0.50359 (hence < 1): "
            << pass(upperLess(derivativeK0Factor, "0.50359")) << "\n";
  bool ok = lowerGreater(valGap, "5.5262e-4")
            && lowerGreater(derivGap, "2.8955e-4")
            && lowerGreater(u2Gap, "3.8992e-7")
            && upperLess(derivativeK0Factor, "0.50359");
  std::cout << pass(ok) << "\n";
  return ok;
}

// -----------------------------------------------------------------------------
// FGR machinery
// -----------------------------------------------------------------------------

struct FgrSums {
  CI a11,a12,a21,a22,a31,a32,a41,a42,a51,a52,a61,a62;
};

static void fgrAddCell(FgrSums& Sums, const Vec& box, const I& width, const CI& iOverA1, const CI& iOverA2) {
  I r = box[0], u = box[1], a = box[2], psi = box[3], psip = box[4];
  I p1 = box[5], p1p = box[6], p2 = box[7], p2p = box[8];
  I mu = box[10];
  I b = (one() - a) / r;
  I up = b * u;
  I q11c1 = -IV("1.5") * ipow(u,3) * sqr(psi) - two()*up*psi*psip - u*sqr(psip);
  I q11c2 = u*up*sqr(psi) + two()*sqr(u)*psi*psip;
  I q22c1 = half()*ipow(u,3)*sqr(psi) + two()*up*psi*psip + u*sqr(psip) - mu*u*sqr(psi);
  I q22c2 = u*up*sqr(psi);
  I q12c1 = -ipow(u,3)*sqr(psi) - two()*up*psi*psip - u*sqr(psip) + half()*mu*u*sqr(psi);
  I q12c2 = sqr(u)*psi*psip;

  CI rawE11 = iOverA1 * (-p1p - p1/(two()*r) - b*p1);
  CI rawE21 = iOverA1 * (u*p1);
  CI rawE12 = iOverA2 * (u*p2);
  CI rawE22 = iOverA2 * (-p2p + p2/(two()*r));
  I wt = width * isqrt(r);

  // The distorted Fourier transform uses conjugate-transpose contraction:
  //   \tilde q = int overline{E(r,k)}^t q(r) r dr.
  // Therefore, for q=(q1,q2)^t, the two scalar amplitudes are
  //   A_col1 = int sqrt(r) (q1 E11 + q2 E21) dr,
  //   A_col2 = int sqrt(r) (q1 E12 + q2 E22) dr,
  // up to the common normalization.  Earlier drafts accidentally accumulated
  // the row-wise combinations q1 E11+q2 E12 and q1 E21+q2 E22.  The variables
  // below are kept with their historical names, but the assignments now store
  // column-wise pieces: a11/a21 form column 1 and a12/a22 form column 2.
  Sums.a11 = Sums.a11 + wt * q11c1 * rawE11;
  Sums.a12 = Sums.a12 + wt * q11c1 * rawE12;
  Sums.a21 = Sums.a21 + wt * q11c2 * rawE21;
  Sums.a22 = Sums.a22 + wt * q11c2 * rawE22;
  Sums.a31 = Sums.a31 + wt * q22c1 * rawE11;
  Sums.a32 = Sums.a32 + wt * q22c1 * rawE12;
  Sums.a41 = Sums.a41 + wt * q22c2 * rawE21;
  Sums.a42 = Sums.a42 + wt * q22c2 * rawE22;
  Sums.a51 = Sums.a51 + wt * q12c1 * rawE11;
  Sums.a52 = Sums.a52 + wt * q12c1 * rawE12;
  Sums.a61 = Sums.a61 + wt * q12c2 * rawE21;
  Sums.a62 = Sums.a62 + wt * q12c2 * rawE22;
}

static Vec fgrSeriesBox(
    const I& r, const I& c, const I& mu, int N,
    const InternalFrobeniusTailCertificate& psiTail,
    const PhiFrobeniusTailCertificate& phi1Tail,
    const PhiFrobeniusTailCertificate& phi2Tail) {
  VortexState v = vortexSeries(r, c, N);
  PsiState ps = psiSeriesUsingCertifiedTail(r, c, mu, N, psiTail);
  I energy = IV("4") * mu;
  PhiState p1 = phiEnergySeriesUsingCertifiedTail(
      r, c, energy, 1, N, phi1Tail);
  PhiState p2 = phiEnergySeriesUsingCertifiedTail(
      r, c, energy, 2, N, phi2Tail);
  Vec y(11);
  y[0]=r; y[1]=v.u; y[2]=v.a; y[3]=ps.psi; y[4]=ps.psip;
  y[5]=p1.phi; y[6]=p1.phip; y[7]=p2.phi; y[8]=p2.phip; y[9]=c; y[10]=mu;
  return y;
}

static I lowerNorm2Sum(const CI& z1, const CI& z2) {
  // This is used only after the complex amplitudes have already been enclosed
  // by CAPD tube quadrature.  It returns an interval whose left endpoint is a
  // rigorous lower bound for |z1|^2+|z2|^2 and whose right endpoint is the usual
  // interval-arithmetic upper bound.  The lower endpoint is computed from the
  // distance of each complex rectangle to the origin, not by blindly trusting
  // the lower endpoint of a wide interval square.
  I lower = norm2LowerBound(z1) + norm2LowerBound(z2);
  I upper = norm2(z1) + norm2(z2);
  return I(lower.leftBound(), upper.rightBound());
}

static I fgrOriginRawComponentBound(const Settings& S, const CI& iOverA1, const CI& iOverA2) {
  // Coded origin lemma for [0,0.1].  It deliberately avoids using the singular
  // point in the algebraic dFT formulas.  On this short interval we use:
  //   U <= r, |U'| <= 1, |psi| <= (1-r^2/4)^(-1), |psi'| <= M r/2,
  //   |Phi_1 terms| <= harmless powers of r, and
  //   |-Phi_2' + Phi_2/(2r)| <= 3 r^(-1/2).
  // The last bound ignores the leading cancellation, so it is crude but safe.
  // It bounds each raw Aij component before the final paper normalization.
  I R0 = S.r0;
  I G1bound = one() / (one() - IV("5") * sqr(R0) / IV("12"));
  I G2bound = one() / (one() - sqr(R0));
  requireTrue(upperLess(G1bound, "1.005"),
              "origin Frobenius bound G_1(0.1) < 1.005");
  requireTrue(upperLess(G2bound, "1.011"),
              "origin Frobenius bound G_2(0.1) < 1.011");
  I rho = sqr(R0) / IV("4");
  I psiM = one() / (one() - rho);
  I dpsiM = psiM * R0 / two();
  I uM = R0;
  I upM = one();
  I muM = one();

  I q11c1 = IV("1.5")*ipow(uM,3)*sqr(psiM) + two()*upM*psiM*dpsiM + uM*sqr(dpsiM);
  I q11c2 = uM*upM*sqr(psiM) + two()*sqr(uM)*psiM*dpsiM;
  I q22c1 = half()*ipow(uM,3)*sqr(psiM) + two()*upM*psiM*dpsiM + uM*sqr(dpsiM) + muM*uM*sqr(psiM);
  I q22c2 = uM*upM*sqr(psiM);
  I q12c1 = ipow(uM,3)*sqr(psiM) + two()*upM*psiM*dpsiM + uM*sqr(dpsiM) + half()*muM*uM*sqr(psiM);
  I q12c2 = sqr(uM)*psiM*dpsiM;
  I qM = maxUpperInterval(q11c1, q11c2);
  qM = maxUpperInterval(qM, q22c1);
  qM = maxUpperInterval(qM, q22c2);
  qM = maxUpperInterval(qM, q12c1);
  qM = maxUpperInterval(qM, q12c2);

  I invM = maxUpperInterval(absUpper(iOverA1), absUpper(iOverA2));
  I rawE22Integral = IV("3") * invM * R0; // ∫_0^R r^{-1/2} sqrt(r) dr = R
  I componentBound = qM * rawE22Integral;
  requireTrue(upperLess(qM, "0.201757527"),
              "FGR origin source bound < 0.201757527");
  requireTrue(upperLess(invM, "2.506371585"),
              "FGR inverse-Wronskian bound < 2.506371585");
  requireTrue(upperLess(componentBound, "0.151703800"),
              "FGR origin integral component bound < 0.151703800");
  requireTrue(upperLess(componentBound, "0.151708"),
              "FGR origin integral component bound < B0");
  std::cout << "FGR origin integral component bound on [0,0.1] = "
            << componentBound << "\n";
  std::cout << "  q-origin bound = " << qM << ", inverse-Wronskian bound = " << invM << "\n";
  return componentBound;
}

static I fgrTailSquaredBound(const Settings& S) {
  I kappa = isqrt(one() - S.delta0 - S.mu_box);
  I R0 = S.r_fgr;
  I k = isqrt(IV("4") * S.mu_box - one());
  // ∫_R^∞ (C K0(kappa r))^2 sqrt(r) dr
  // K0(kappa r)^2 <= pi/(2 kappa r) exp(-2 kappa r)
  // so integrand <= C^2*pi/(2*kappa)*r^{-1/2} exp(-2kappa r)
  // <= same with R^{-1/2}; integrate exponential exactly.
  I tailInt = sqr(S.k0_comparison_C) * piI() / (two()*kappa*isqrt(R0)) * iexp(-two()*kappa*R0) / (two()*kappa);
  I bound = piI() / (IV("4") * S.mu_box) * k * sqr(S.fgr_source_C * tailInt);
  return bound;
}

static I fgrSourceConstantCheck(const Settings& S,
                                const I& weylTailEnvelope) {
  // With the distorted Fourier normalization used in the paper, the exact
  // connection formula gives the source bound
  //
  //   10 * WeylBound / (lambda * sqrt(pi)).
  //
  // The outgoing Volterra calculation supplies WeylBound uniformly
  // for mu in the certified FGR interval.
  I lambda = isqrt(S.mu_box);
  I paperBound =
      IV("10") * weylTailEnvelope / (lambda * isqrt(piI()));

  std::cout << "FGR source constant numerical check\n";
  std::cout << "  outgoing Weyl upper bound for r >= 16 <= "
            << weylTailEnvelope << "\n";
  std::cout << "  10*upper_bound/(lambda*sqrt(pi)) <= "
            << paperBound << "\n";
  std::cout << "  targets < 8.02 and < 30: "
            << pass(upperLess(paperBound, "8.02")
                    && upperLess(paperBound, "30")) << "\n";

  requireTrue(upperLess(paperBound, "8.02"),
              "paper FGR tail source constant < 8.02");
  requireTrue(upperLess(paperBound, "30"),
              "FGR tail source constant < 30");
  return paperBound;
}

static bool checkFGR(
    const Settings& S,
    const InternalFrobeniusTailCertificate& psiTail,
    const PhiFrobeniusTailCertificate& phi1Tail,
    const PhiFrobeniusTailCertificate& phi2Tail) {
  std::cout << "\n[7] FGR coefficients: compact integral, origin, and tail\n";
  std::cout << "The outgoing start is certified by a residual/Volterra estimate\n"
               "for r >= 16.  The source constant uses the resulting\n"
               "Weyl tail bound, the finite channel count, and the explicit\n"
               "source-polynomial coefficient bound.\n";
  I k = isqrt(IV("4") * S.mu_box - one());

  // Regular zero-energy/FGR basis from origin to r=10 and r=16.
  Vec y0 = fgrSeriesBox(
      S.r0, S.c_cert, S.mu_box, S.series_order,
      psiTail, phi1Tail, phi2Tail);
  Vec yMatch = integrate(VF_ZERO_FGR_FWD, y0, S.r_match - S.r0, S.taylor_order);
  Vec yR = integrate(VF_ZERO_FGR_FWD, y0, S.r_fgr - S.r0, S.taylor_order);

  // Outgoing solutions from r=16 backward to r=10.  The start boxes use the
  // interval-evaluated outgoing ansatz plus the residual/Volterra tail bound
  // computed in outgoingDataAtR.
  OutgoingStart out1 = outgoingDataAtR(S, 1, k);
  OutgoingStart out2 = outgoingDataAtR(S, 2, k);
  requireTrue(upperLess(out1.residualIntegral, "1.642e-12")
              && upperLess(out2.residualIntegral, "1.500e-12"),
              "outgoing residual-integral bounds");
  requireTrue(upperLess(out1.potentialIntegral, "2.251e-7")
              && upperLess(out2.potentialIntegral, "2.251e-7"),
              "outgoing potential-integral bounds");
  requireTrue(upperLess(out1.contraction, "3.228e-2")
              && upperLess(out2.contraction, "1.076e-2"),
              "outgoing Volterra contraction bounds");
  requireTrue(upperLess(out1.valueError, "1.350273e-7")
              && upperLess(out1.derivError, "3.922666e-7")
              && upperLess(out2.valueError, "1.306810e-7")
              && upperLess(out2.derivError, "3.796404e-7"),
              "outgoing initial-data error bounds");
  std::cout << "outgoing channel 1 value error = " << out1.valueError
            << ", derivative error = " << out1.derivError
            << ", contraction = " << out1.contraction << "\n";
  std::cout << "outgoing channel 2 value error = " << out2.valueError
            << ", derivative error = " << out2.derivError
            << ", contraction = " << out2.contraction << "\n";

  Vec b10(9), b20(9);
  b10[0]=S.r_fgr; b10[1]=yR[1]; b10[2]=yR[2]; b10[3]=out1.y.re; b10[4]=out1.y.im; b10[5]=out1.yp.re; b10[6]=out1.yp.im; b10[7]=S.c_cert; b10[8]=S.mu_box;
  b20[0]=S.r_fgr; b20[1]=yR[1]; b20[2]=yR[2]; b20[3]=out2.y.re; b20[4]=out2.y.im; b20[5]=out2.yp.re; b20[6]=out2.yp.im; b20[7]=S.c_cert; b20[8]=S.mu_box;
  Vec b1 = integrate(vfOutgoing(1,true), b10, S.r_fgr - S.r_match, S.taylor_order);
  Vec b2 = integrate(vfOutgoing(2,true), b20, S.r_fgr - S.r_match, S.taylor_order);

  CI Psi1(b1[3], b1[4]), Psi1p(b1[5], b1[6]);
  CI Psi2(b2[3], b2[4]), Psi2p(b2[5], b2[6]);
  CI W1 = CI(yMatch[5]) * conj(Psi1p) - CI(yMatch[6]) * conj(Psi1);
  CI W2 = CI(yMatch[7]) * conj(Psi2p) - CI(yMatch[8]) * conj(Psi2);
  CI a1 = mulI(W1) / two();
  CI a2 = mulI(W2) / two();
  std::cout << "Wronskian convention: the regular origin-normalized basis is real;\n"
               "  the outgoing Weyl basis is represented by real/imaginary CAPD\n"
               "  components, and the following a_j are complex interval boxes.\n";
  std::cout << "a1 = " << a1 << "\n";
  std::cout << "a2 = " << a2 << "\n";
  printCIWidth("a1", a1);
  printCIWidth("a2", a2);
  I weylTailEnvelope = maxUpperInterval(out1.matrixTailBound, out2.matrixTailBound);
  requireTrue(upperAtMost(weylTailEnvelope, "1.253"),
              "uniform outgoing Weyl-vector bound <= 1.253");
  fgrSourceConstantCheck(S, weylTailEnvelope);
  CI iOverA1 = mulI(inv(a1));
  CI iOverA2 = mulI(inv(a2));
  I originComponentBound = fgrOriginRawComponentBound(S, iOverA1, iOverA2);

  FgrSums sums;

  // Compact [0.1,16] by CAPD tubes.  The omitted origin contribution [0,0.1]
  // is handled by fgrOriginRawComponentBound and subtracted from the final
  // norm lower bounds.
  Vec x = y0;
  I hc = (S.r_fgr - S.r0) / IV(std::to_string(S.fgr_cells).c_str());
  for (int j = 0; j < S.fgr_cells; ++j) {
    StepResult st = integrateTube(VF_ZERO_FGR_FWD, x, hc, S.taylor_order);
    fgrAddCell(sums, st.tube, hc, iOverA1, iOverA2);
    x = st.endpoint;
  }

  I rawFactor = piI() / (IV("4") * S.mu_box);
  I norm2Factor = k / (IV("16") * piI() * S.mu_box);
  I Ffactor = rawFactor * norm2Factor;

  CI amp11 = sums.a11 + sums.a21;
  CI amp12 = sums.a12 + sums.a22;
  CI amp21 = sums.a31 + sums.a41;
  CI amp22 = sums.a32 + sums.a42;
  CI amp121 = sums.a51 + sums.a61;
  CI amp122 = sums.a52 + sums.a62;
  // These stored vectors are the complex conjugates of the A_ij vectors in
  // the paper.  Their component widths and Euclidean norms are identical.
  printCIWidth("conjugate(A_11([0.1,16])) component 1", amp11);
  printCIWidth("conjugate(A_11([0.1,16])) component 2", amp12);
  printCIWidth("conjugate(A_22([0.1,16])) component 1", amp21);
  printCIWidth("conjugate(A_22([0.1,16])) component 2", amp22);
  printCIWidth("conjugate(A_12([0.1,16])) component 1", amp121);
  printCIWidth("conjugate(A_12([0.1,16])) component 2", amp122);

  I S1 = lowerNorm2Sum(amp11, amp12);
  I S2 = lowerNorm2Sum(amp21, amp22);
  I S12 = lowerNorm2Sum(amp121, amp122);
  I F1 = Ffactor * S1;
  I F2 = Ffactor * S2;
  I F12 = Ffactor * S12;
  I tailB = fgrTailSquaredBound(S);
  I originRawNorm = isqrt(two()) * originComponentBound;
  I originB = Ffactor * sqr(originRawNorm);

  auto lowerAfterTail = [&](const I& F) -> I {
    R fl = F.leftBound();
    if (fl < R(0)) fl = R(0);
    I rootF = isqrt(I(fl));
    I rootB = isqrt(tailB) + isqrt(originB);
    I d = rootF - rootB;
    if (!positive(d)) return zero();
    return sqr(d);
  };

  I F1low = lowerAfterTail(F1);
  I F2low = lowerAfterTail(F2);
  I F12low = lowerAfterTail(F12);

  std::cout << "hat D_1,[0.1,16]^(0) interval = " << F1 << "\n";
  std::cout << "hat D_2,[0.1,16]^(0) interval = " << F2 << "\n";
  std::cout << "hat D_12,[0.1,16]^(0) interval = " << F12 << "\n";
  printWidth("hat D_1,[0.1,16]^(0) interval", F1);
  printWidth("hat D_2,[0.1,16]^(0) interval", F2);
  printWidth("hat D_12,[0.1,16]^(0) interval", F12);
  std::cout << "c(mu) = k/(64*mu^2) = " << Ffactor << "\n";
  std::cout << "B_origin = " << originB << "\n";
  std::cout << "B_tail = " << tailB << "\n";
  std::cout << "The vector integrals over (0,0.1), [0.1,16], and (16,infinity)\n"
               "  are combined by the triangle inequality before the\n"
               "  resulting lower bound is squared.\n";
  std::cout << "lower bound for hat D_1^(0) = " << F1low << "\n";
  std::cout << "lower bound for hat D_2^(0) = " << F2low << "\n";
  std::cout << "lower bound for hat D_12^(0) = " << F12low << "\n";

  bool compactOk = strictSubsetOfDecimalInterval(
                       F1, "0.0176262802754", "0.0332324717579")
                && strictSubsetOfDecimalInterval(
                       F2, "0.0477905323094", "0.0546253884722")
                && strictSubsetOfDecimalInterval(
                       F12, "0.0233179905311", "0.0287364386010");
  bool remainderOk = upperLess(originComponentBound, "0.151708")
                  && upperLess(originB, "1.728239729e-3")
                  && upperLess(originB, "1.72834e-3")
                  && upperLess(tailB, "3.174707028e-10")
                  && upperLess(tailB, "1e-9");
  bool finalOk = lowerGreater(F1low, "0.008312713")
              && lowerGreater(F2low, "0.031336270")
              && lowerGreater(F12low, "0.012345959");
  bool ok = compactOk && remainderOk && finalOk;
  std::cout << pass(ok) << "\n";
  return ok;
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

} // namespace lppssi

int main() {
  using namespace lppssi;
  try {
#ifdef LPPSSI_USE_MP
    R::setDefaultPrecision(Settings::default_precision_bits);
#endif
    Settings S;
    if (const char* tubePath = std::getenv("LPPSSI_EXPORT_VORTEX_TUBE")) {
      int cells = 1590; // h=0.01 on [0.1,16]
      if (const char* rawCells = std::getenv("LPPSSI_EXPORT_VORTEX_TUBE_CELLS")) {
        cells = std::atoi(rawCells);
      }
      if (cells <= 0) {
        throw std::runtime_error("LPPSSI_EXPORT_VORTEX_TUBE_CELLS must be positive");
      }
      emitVortexProfileTubeCsv(S, tubePath, cells);
      std::cout << "wrote vortex profile tube CSV: " << tubePath
                << " with cells=" << cells << "\n";
      return EXIT_SUCCESS;
    }

    std::cout.precision(30);
    std::cout << "LPPSSI C++ CAPD numerical certificate\n";
#ifdef LPPSSI_USE_MP
    std::cout << "interval backend = MpInterval/MpFloat, requested precision bits = "
              << S.precision_bits << "\n";
#else
    std::cout << "interval backend = DInterval/double, effective significand bits = 53\n";
    std::cout << "requested multiprecision setting (" << S.precision_bits
              << " bits) is inactive because LPPSSI_USE_MP is not defined\n";
#endif
    std::cout << "Taylor order = " << S.taylor_order << "\n";
    std::cout << "c box = " << S.c_box << "\n";
    std::cout << "initial nominal mu box = " << S.mu_box << "\n";
    std::cout << "  (not trusted for FGR until block [5b] certifies and installs\n"
                 "   a validated Wronskian sign-changing sub-box)\n";
    std::cout << "dependency policy: fail-fast.  A later block is run only if\n"
                 "all earlier prerequisite blocks in this same run have PASSed.\n";

    auto requireCheck = [&](const char* label, const auto& fn) {
      try {
        bool ok = fn();
        if (!ok) {
          throw std::runtime_error(std::string(label) + " returned FAIL");
        }
      } catch (const std::exception& e) {
        std::cerr << "\n[" << label << "] exception: " << e.what() << "\n";
        std::cerr << "Status for this block: FAIL\n";
        std::cerr << "Dependent later checks are suppressed.\n";
        throw;
      }
    };

    requireCheck("0 directed decimal inputs", [&]() { return checkDecimalInputEnclosures(); });
    requireCheck("1 vortex shooting", [&]() { return checkVortexBracket(S); });
    requireCheck("2 H1 positivity", [&]() { return checkH1(S); });
    requireCheck("3 Seto logarithmic-kernel bound", [&]() { return checkLT(S); });
    PhiFrobeniusTailCertificate thresholdPhi2Tail{
        0, 0, zero(), zero(), zero(), zero(), zero()};
    requireCheck("4a threshold-origin Frobenius series tail", [&]() {
      return checkThresholdOriginTailBootstrapData(S, thresholdPhi2Tail);
    });
    requireCheck("4 threshold Wronskian", [&]() {
      return checkThreshold(S, thresholdPhi2Tail);
    });
    InternalFrobeniusTailCertificate broadPsiTail{
        0, zero(), zero(), zero(), zero(), zero()};
    requireCheck("4b internal-origin Frobenius series tail", [&]() {
      return checkInternalModeTailBootstrapData(
          S, S.mu_eigen_box, "broad eigenvalue search", broadPsiTail);
    });
    requireCheck("5 internal eigenvalue", [&]() {
      return checkInternalEigenvalueBox(S, broadPsiTail);
    });
    requireCheck("5 FGR eigenvalue localization", [&]() {
      return certifyAndInstallFgrEigenvalueBox(S, broadPsiTail);
    });
    InternalFrobeniusTailCertificate fgrPsiTail{
        0, zero(), zero(), zero(), zero(), zero()};
    PhiFrobeniusTailCertificate fgrPhi1Tail{
        0, 0, zero(), zero(), zero(), zero(), zero()};
    PhiFrobeniusTailCertificate fgrPhi2Tail{
        0, 0, zero(), zero(), zero(), zero(), zero()};
    requireCheck("5c Frobenius series-tail data", [&]() {
      return checkFrobeniusTailBootstrapData(
          S, fgrPsiTail, fgrPhi1Tail, fgrPhi2Tail);
    });
    printSpectralOriginTailRadii(
        S, fgrPsiTail, fgrPhi1Tail, fgrPhi2Tail);
    requireCheck("6 K0 comparison", [&]() {
      return checkInternalK0(S, fgrPsiTail);
    });
    std::cout << "\nPrerequisites for FGR are certified in this run:\n"
                 "  refined vortex c_* box, H1 positivity, LT uniqueness,\n"
                 "  the broad-J0 r=4 entry check for Lemma lem_1-U2,\n"
                 "  threshold non-resonance, the broad internal eigenvalue\n"
                 "  bracket, the installed FGR eigenvalue sub-bracket, and the\n"
                 "  K0 comparison, including the psi' bound for r >= 16.\n"
                 "  The outgoing Weyl bound for r >= 16 and the remaining FGR\n"
                 "  source constant are checked inside the FGR block.\n";
    requireCheck("7 FGR", [&]() {
      return checkFGR(S, fgrPsiTail, fgrPhi1Tail, fgrPhi2Tail);
    });

    std::cout << "\nOVERALL: PASS FOR THE NUMERICAL CERTIFICATE BLOCKS\n";
    std::cout << "Interpretation: this verifies the CAPD/interval numerical\n"
                 "hypotheses, conditional on the analytic lemmas explicitly\n"
                 "identified in the paper and comments.\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "\nException: " << e.what() << "\n";
    std::cerr << "\nOVERALL: FAIL\n";
    return EXIT_FAILURE;
  }
}
