# A certified lower bound for Erdős' minimum overlap constant

This directory contains an explicit certificate for the bound

> **μ ≥ 0.380020**

for the constant μ = lim M(n)/n in Erdős' minimum overlap problem, together with a self-contained
script that checks it. The previous certified lower bound was 0.37912 (arXiv:2606.31182); the best
known upper bound is 0.380868. This closes about 51% of the remaining gap.

The accompanying paper is hosted separately; this repository contains only the certificate
data and the verifier, which are self-contained and need nothing from the paper to run.

## Verify it

```
python3 verify_certificate.py data/certificates.json
```

Requires Python 3.8+ and `numpy`. It prints one line per box and ends with the verified bound, or
fails loudly.

The defaults are exactly the settings used in the paper, so the
command above reproduces the paper's table. Almost all of that time is the rigorous interval
evaluation of the envelope integral on two million grid points.

```
verify_certificate.py <json> [N_final] [N_search] [digits]
                              2000000    = N_final     5
```

* `N_final` — grid for the rigorous check. Controls the Lipschitz discretisation term `L·delta`,
  and so how much margin is left over. Runtime is linear in it. Smaller values still *verify*;
  they just leave less margin. Passing `400000` is about five times faster.
* `N_search` — grid used to propose a level, which is then re-checked rigorously at `N_final`.
  It defaults to `N_final` and should normally be left alone: a coarser search grid charges the
  Lipschitz term at its own spacing and so proposes a needlessly low level, capping the reported
  bound for reasons unrelated to the certificate.
* `digits` — decimals to which the level is rounded down. Coarser rounding reports a slightly
  smaller bound but leaves proportionally more margin. At `digits=6` this certificate reaches
  `0.380026` with margin `2.5e-6`; the published `digits=5` reports `0.380020` with margin
  `1.8e-5`, and we prefer the headroom to the 6e-6.

A second certificate at the coarser spacing 1/16 is included as `data/certificates_1_16.json`.
It verifies the weaker bound `0.379815` by the same route and is useful as an independent
cross-check of the verifier:

```
python3 verify_certificate.py data/certificates_1_16.json 400000 200000 6    # -> mu >= 0.379815
```

## What is being verified

Write `f: [-1,1] → [0,1]` with `∫f = 1`, `g = (1-f)·1_[-1,1]`, `M(x) = ∫f(t)g(x+t)dt`,
`Ω = ‖M‖_∞`, `E = ∫ t f(t) dt`, and split `M = M̄ + O` into even and odd parts. Let
`S(ξ) = sin(2πξ)/(πξ)` and `Λ(x) = (2-|x|)₊`.

Every admissible `f` satisfies these **exact** facts:

| | |
|---|---|
| F1 | `∫ M̄ = 1` |
| F2 | `∫ x² M̄ = 2/3 + 2E²` |
| F3 | `0 ≤ M(x), M(-x) ≤ Ω` and `M(x) + M(-x) ≤ 2 - \|x\|` |
| F4 | `∫ M̄(x) cos(2πξx) dx ≤ S(ξ)²/4` for **every real** ξ |
| F5 | `∫ x O(x) dx = -2E` |
| F6 | `4S(ξ)²·M̄̂(ξ) + 8τ·Õ(ξ) - 4τ² ≤ S(ξ)⁴` for every real ξ and every real τ |

F4 is the identity `4M̄ = Λ - C_ψ` for `ψ = 2f-1` combined with `Ĉ_ψ = |ψ̂|² ≥ 0`. Because `C_ψ` is
compactly supported this holds at every real ξ, not only at the integer indices of the classical
period-4 expansion — that is the new ingredient. F6 is the tangent-line form of
`4S²M̄̂ + 4Õ² ≤ S⁴`, and holds for *any* τ, so the numerical τ in the certificate needs no
justification.

A certificate is a vector of multipliers

```
λ ∈ ℝ,   ν₊, ν₋, κ₊, κ₋ ≥ 0,   c₁..c_m ≥ 0,   γ₁..γ_m ≥ 0
```

together with frequencies `ξ_j = j/64`, `j = 1..768`, and tangent points `τ_j`. Setting

```
G_e(x) = λ + (ν₊-ν₋)x² - Σ c_j cos(2πξ_j x) - Σ 4γ_j S_j² cos(2πξ_j x)
G_o(x) = (κ₊-κ₋)x       - Σ 8γ_j τ_j sin(2πξ_j x)
```

and `α = (G_e+G_o)/2`, `β = (G_e-G_o)/2`, the facts give `V ≤ Σ` and F3 gives `Σ ≤ 2∫₀²h`, where
`Σ = ∫(M̄G_e + O G_o)`, `V` is the explicit constant assembled from the right-hand sides, and

```
h(x) = max { α(x)p + β(x)q : 0 ≤ p, q ≤ Ω,  p + q ≤ 2 - x }
```

is a maximum over the five vertices of a pentagon. **If `V > 2∫₀²h` then no admissible `f` with
`E` in the box and `‖M‖_∞ ≤ Ω` exists.** The certificate has **no feasibility constraints beyond
the signs above**, which is what makes exact verification easy.

`E` is confined by the bathtub principle: `Var(M) = 2/3 - 2E² ≥ 1/(12Ω²)`, so
`E² ≤ 1/3 - 1/(24Ω²)`; and `f(t) ↦ f(-t)` lets us take `E ≥ 0`. The verifier checks that the boxes
genuinely tile `[0, E_max]` — that they start at 0, leave no gap, and reach past `E_max`.

## What the verifier trusts

Only correctly rounded IEEE-754 binary64 basic arithmetic and `numpy.nextafter`; the verifier
asserts the expected binary64 format at runtime. In particular, it does not trust a library
implementation of `cos` or `sin`:

* **π** is derived inside the script from Machin's formula `π/4 = 4·arctan(1/5) - arctan(1/239)`
  with an explicit alternating-series truncation bound, giving an enclosure of width ~10⁻⁸⁶.
* **√2** is a hardcoded rational enclosure that the script checks by squaring.
* **`S(ξ_j)²`** is enclosed in exact rational arithmetic for any denominator divisible by 8, via
  `sin²t = (1-cos 2t)/2` (which needs no sign handling). The reduction to the first octant is exact
  integer arithmetic; the residual angle is enclosed by an alternating Taylor series with an
  explicit truncation bound, rounded outwards to a fixed denominator. As an independent check the
  script confirms at startup that this routine reproduces the closed-form values
  `0, (2-√2)/4, 1/2, (2+√2)/4, 1` that are available at spacing 1/16 — which is what the √2
  enclosure above is for.
* **`V`** is computed in exact rational arithmetic, using the *upper* enclosure of `S²` wherever it
  appears with a nonpositive coefficient, so the computed `V` is a rigorous lower bound.
* **`∫h`** is bounded by the midpoint rule plus `L·δ/2` with `L` an explicitly bounded Lipschitz
  constant for `h`. Exact integer range reduction puts every required angle in `[0,π/4]`; alternating
  Taylor partial sums enclose sine and cosine there. Every subsequent operation is rounded outwards,
  and the nonnegative grid sum uses upward-rounded pairwise addition.
* **Sign conditions** are checked exactly on the stored rationals.
* The ordinary floating-point search grid only *proposes* a level. The final interval pass
  independently checks it, so search-grid arithmetic cannot validate a false bound.

## Files

```
data/certificates.json        the published certificate (spacing 1/64, 768 frequencies): one
                              record per E-box giving the box, the frequency denominator, and
                              λ, ν±, κ±, c, γ, τ as exact rationals [numerator, denominator]
data/certificates_1_16.json   an earlier certificate at spacing 1/16, kept as a cross-check
verify_certificate.py         the verifier described above
requirements.txt              numpy
```

Each record's `c`, `gam`, `tau` arrays are indexed by `j = 1..n_freq` with `ξ_j = j/xi_den`.

## Reproducing the certificate

The certificate was produced by a linear programme (normalise `V = 1`, minimise `2∫h` on a grid),
whose solution was rounded to dyadic rationals with denominator 2³⁰ and sign-clamped, at
frequency spacing 1/64 up to ξ = 12. The
optimisation grid affects only the *quality* of the multipliers, never the validity of the
conclusion — any sign-respecting vector yields a true implication. Generation code is not required
to check the result.

## Verifier output

```
pi enclosed to width 7.03e-86 by Machin's formula (checked internally)
sqrt(2) enclosure checked by squaring
S^2 Taylor enclosures agree with closed-form sqrt(2) values (den=16 self-test)
search grid N = 2000000, final verification grid N = 2000000

  E in [0.00000,0.02690]  certified Omega = 0.380020  margin +1.83e-05  OK  (L=12.6)
  E in [0.02690,0.05379]  certified Omega = 0.382780  margin +5.60e-06  OK  (L=21.8)
  E in [0.05379,0.08069]  certified Omega = 0.392960  margin +2.02e-05  OK  (L=23.5)
  E in [0.08069,0.10758]  certified Omega = 0.399990  margin +9.75e-03  OK  (L=23.3)
  E in [0.10758,0.16137]  certified Omega = 0.399990  margin +2.49e-02  OK  (L=24.1)
  E in [0.16137,0.21516]  certified Omega = 0.399990  margin +7.35e-02  OK  (L=22.8)

  all boxes excluded: True
  E-range: any f with ||M||_inf <= 0.380020 has E^2 <= 0.04481338,
           i.e. |E| <= 0.211692; boxes tile [0,0.215164] with no gap: True   OK

  VERIFIED:  mu >= 0.380020
```

The margin on each box is `V - 2∫h`, in the units of the normalisation `V = 1`, after both the
outward-rounded interval grid evaluation and the Lipschitz discretisation term.
