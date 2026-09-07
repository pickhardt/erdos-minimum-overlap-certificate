#!/usr/bin/env python3
"""
Independent verification of the certified lower bound for Erdos' minimum overlap constant.

    python3 verify_certificate.py data/certificates.json

Everything the proof depends on is recomputed here from the stored rational data:
  * pi is obtained from Machin's formula with an explicit alternating-series truncation bound,
    and sqrt(2) from a rational enclosure that the script checks by squaring.  No constant is
    taken on trust.
  * the quantity V is evaluated in exact rational arithmetic, rounded in the safe direction.
  * the envelope integral is evaluated by outward-rounded IEEE-754 interval arithmetic.
    The only transcendental values needed are sin(theta) and cos(theta) for
    theta in [0,pi/4]; these are enclosed by rationally justified Taylor bounds.
    Every harmonic is evaluated after exact integer range reduction, and the grid
    sum is accumulated with upward-rounded pairwise addition.
The certificate itself carries NO feasibility constraints beyond the sign conditions
lam free, nu+ >= 0, nu- >= 0, kap+ >= 0, kap- >= 0, c_j >= 0, gam_k >= 0, which are checked
exactly on the rationals.  See README.md for the mathematics.
"""
import sys, json
from fractions import Fraction as F

# ---------------------------------------------------------------- exact constants
def pi_enclosure(nterms=60):
    """Machin: pi/4 = 4 arctan(1/5) - arctan(1/239).  Both series alternate with terms
    decreasing in absolute value, so truncating after n terms leaves an error bounded by
    the first omitted term.  Returns (lo, hi) rationals with lo < pi < hi."""
    def atan_inv(x, n):                      # arctan(1/x), alternating
        s = F(0); lo = hi = None
        for k in range(n):
            t = F((-1)**k, (2*k+1) * x**(2*k+1))
            s += t
        nxt = F(1, (2*n+1) * x**(2*n+1))      # magnitude of the first omitted term
        return s, nxt
    a5, e5 = atan_inv(5, nterms)
    a239, e239 = atan_inv(239, nterms)
    mid = 4*(4*a5 - a239)
    err = 4*(4*e5 + e239)
    return mid - err, mid + err

SQRT2_LO = F(1414213562373095048801688724209698078569, 10**39)
SQRT2_HI = SQRT2_LO + F(1, 10**39)
assert SQRT2_LO**2 < 2 < SQRT2_HI**2, "sqrt(2) enclosure fails its own check"

# Fixed-denominator outward rounding.  Every rational result is snapped to a multiple of
# 1/RD, away from the true value.  Without this the exact Taylor series below would build
# fractions with thousands of digits and be unusably slow.
RD = 1 << 192                                       # resolution ~1.6e-58

def _rfl(q):
    r = q*RD; return F(r.numerator // r.denominator, RD)
def _rce(q):
    r = q*RD; return F(-((-r.numerator) // r.denominator), RD)
def _radd(a, b): return (_rfl(a[0]+b[0]), _rce(a[1]+b[1]))
def _rsub(a, b): return (_rfl(a[0]-b[1]), _rce(a[1]-b[0]))
def _rmul(a, b):
    p = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return (_rfl(min(p)), _rce(max(p)))
def _rdivn(a, n): return (_rfl(F(a[0], n)), _rce(F(a[1], n)))

_NT = 12            # for |a| <= pi/4 the first omitted term a^24/24! is ~5e-27

def _trig_rational(a):
    """Enclose (cos a, sin a) for rational a in [0, pi/4].  Both series alternate with terms
    decreasing in absolute value, so the first omitted term bounds the truncation error."""
    A = (a, a); A2 = _rmul(A, A)
    term = (F(1), F(1)); s = (F(1), F(1))
    for k in range(1, _NT):
        term = _rdivn(_rmul(term, A2), (2*k-1)*(2*k))
        s = _rsub(s, term) if k & 1 else _radd(s, term)
    E = _rce(_rdivn(_rmul(term, A2), (2*_NT-1)*(2*_NT))[1])
    cos = (s[0]-E, s[1]+E)
    term = A; s = A
    for k in range(1, _NT):
        term = _rdivn(_rmul(term, A2), (2*k)*(2*k+1))
        s = _rsub(s, term) if k & 1 else _radd(s, term)
    E = _rce(_rdivn(_rmul(term, A2), (2*_NT)*(2*_NT+1))[1])
    sin = (s[0]-E, s[1]+E)
    return cos, sin

def _cos_2pi(s, den, pilo, pihi):
    """Enclose cos(2*pi*s/den) for integer s, requiring 8 | den.  Reduction to the first
    octant is exact integer arithmetic; only the final small angle involves pi."""
    assert den % 8 == 0, "frequency denominator must be divisible by 8"
    q, hq = den//4, den//8
    num = s % den
    k = ((num + hq)//q) % 4
    rem = num - k*q
    if rem >  den//2: rem -= den
    if rem < -den//2: rem += den
    a = abs(rem); sg = -1 if rem < 0 else 1
    ylo = _rfl(2*pilo*a/den); yhi = _rce(2*pihi*a/den)           # |angle| in [0, pi/4]
    (clo_h, _), (_, shi_h) = _trig_rational(yhi)                 # cos decreasing,
    (_, chi_l), (slo_l, _) = _trig_rational(ylo)                 # sin increasing
    cy = (clo_h, chi_l); sy = (slo_l, shi_h)
    if   k == 0: return cy
    elif k == 1: return (-sy[1], -sy[0]) if sg > 0 else sy
    elif k == 2: return (-cy[1], -cy[0])
    else:        return sy if sg > 0 else (-sy[1], -sy[0])

def S2_bounds(j, den, pilo, pihi):
    """Rational enclosure of S(xi)^2 = sin(2 pi xi)^2/(pi xi)^2 at xi = j/den, any 8 | den.
    Uses sin^2 t = (1 - cos 2t)/2, which needs no sign handling."""
    clo, chi = _cos_2pi(2*j, den, pilo, pihi)
    lo2 = max((1-chi)/2, F(0)); hi2 = (1-clo)/2
    k = F(den*den, j*j)
    return _rfl(lo2*k/(pihi*pihi)), _rce(hi2*k/(pilo*pilo))

def s2_selftest(pilo, pihi):
    """Cross-check the general Taylor path against closed-form algebraic values.  At den = 16,
    2 pi xi = pi j/8, so sin^2 is one of 0, (2-sqrt2)/4, 1/2, (2+sqrt2)/4, 1 -- built here from
    the independently squared sqrt(2) enclosure.  An independent check of the trig routine."""
    for j in range(1, 17):
        r = j % 16
        if r in (0, 8):            e = (F(0), F(0))
        elif r in (4, 12):         e = (F(1), F(1))
        elif r in (2, 6, 10, 14):  e = (F(1,2), F(1,2))
        elif r in (1, 7, 9, 15):   e = ((2-SQRT2_HI)/4, (2-SQRT2_LO)/4)
        else:                      e = ((2+SQRT2_LO)/4, (2+SQRT2_HI)/4)
        k = F(256, j*j)
        ref = (e[0]*k/(pihi*pihi), e[1]*k/(pilo*pilo))
        got = S2_bounds(j, 16, pilo, pihi)
        assert got[0] <= got[1], f"S2 enclosure inverted at j={j}"
        assert got[0] <= ref[1] and ref[0] <= got[1], f"S2 self-test disagrees at j={j}"
    return True

# ---------------------------------------------------------------- the two bounds
def V_lower(C, S2hi):
    """Lower bound on V = lam + nu+(2/3+2Elo^2) - nu-(2/3+2Ehi^2) - sum c_j S_j^2/4
       - 2 kap+ Ehi + 2 kap- Elo - sum gam_k (S_k^4 + 4 tau_k^2).
       c, gam >= 0 multiply subtracted terms, so use the UPPER bounds on S^2, S^4."""
    Elo, Ehi = C['Elo'], C['Ehi']
    V = C['lam'] + C['nup']*(F(2,3)+2*Elo*Elo) - C['num']*(F(2,3)+2*Ehi*Ehi)
    V += -2*C['kapp']*Ehi + 2*C['kapm']*Elo
    for cj, hi in zip(C['c'], S2hi):
        V -= cj*hi/4
    for gk, hi, t in zip(C['gam'], S2hi, C['tau']):
        V -= gk*(hi*hi + 4*t*t)
    return V

def sup_bounds(C, S2hi, pihi):
    """Rigorous suprema over x in [0,2] of |Ge|, |Go|, |Ge'|, |Go'| (used for the Lipschitz
    constant), from |cos|,|sin| <= 1 and |x| <= 2."""
    nu = C['nup'] - C['num']; kap = C['kapp'] - C['kapm']
    sc  = sum(C['c']); sg4 = sum(g*4*h for g, h in zip(C['gam'], S2hi))
    Ge  = abs(C['lam']) + 4*abs(nu) + sc + sg4
    Go  = 2*abs(kap) + sum(g*8*abs(t) for g, t in zip(C['gam'], C['tau']))
    Ged = 4*abs(nu) + sum(c*2*pihi*x for c, x in zip(C['c'], C['xis'])) \
                    + sum(g*8*pihi*h*x for g, h, x in zip(C['gam'], S2hi, C['xis']))
    God = abs(kap) + sum(g*16*pihi*abs(t)*x for g, t, x in zip(C['gam'], C['tau'], C['xis']))
    return Ge, Go, Ged, God

def envelope_sum(C, Om, N, S2mid, max_elems=8_000_000):
    """Fast non-rigorous grid sum, used only to propose a level for final verification.

    Chunked so the transient m x chunk trig arrays stay near max_elems entries.  With 768
    frequencies a fixed 100k chunk would allocate ~600 MB per array, which is what made this
    path unusable at the finer spacing."""
    import numpy as np
    xis = np.array([float(x) for x in C['xis']]); S2 = np.array([float(s) for s in S2mid])
    c = np.array([float(x) for x in C['c']]); gam = np.array([float(x) for x in C['gam']])
    tau = np.array([float(x) for x in C['tau']])
    lam = float(C['lam']); nu = float(C['nup']-C['num']); kap = float(C['kapp']-C['kapm'])
    Omf = float(Om); dx = 2.0/N; tot = 0.0
    ce = c + 4*gam*S2; co = 8*gam*tau                  # combined per-harmonic coefficients
    chunk = min(100000, N)
    ang = np.empty(chunk); tmp = np.empty(chunk)
    for s0 in range(0, N, chunk):
        n = min(s0+chunk, N) - s0
        x = (np.arange(s0, s0+n) + 0.5)*dx
        Ge = lam + nu*x**2; Go = kap*x
        a = ang[:n]; t = tmp[:n]
        for j in range(len(xis)):                      # O(chunk) memory, not O(m*chunk)
            if not (ce[j] or co[j]):                   # zero multiplier contributes nothing
                continue
            np.multiply(x, 2*np.pi*xis[j], out=a)
            if ce[j]: np.cos(a, out=t); Ge -= ce[j]*t
            if co[j]: np.sin(a, out=t); Go -= co[j]*t
        al = (Ge+Go)/2; be = (Ge-Go)/2
        U2 = 2-x; A = np.minimum(Omf, U2)
        B = np.minimum(Omf, np.maximum(U2-Omf, 0.0)); mk = U2 > Omf
        h = np.maximum.reduce([np.zeros_like(x), al*A, be*A,
                               np.where(mk, al*Omf+be*B, -np.inf),
                               np.where(mk, al*B+be*Omf, -np.inf)])
        tot += float(h.sum())
    return dx*tot

# -------------------------------------------------- outward-rounded interval grid evaluation
def _down(x):
    import numpy as np
    return np.nextafter(x, -np.inf)

def _up(x):
    import numpy as np
    return np.nextafter(x, np.inf)

def _qinterval(q):
    """Smallest easily computed binary64 interval containing the rational q."""
    import numpy as np
    v = float(q)
    fv = F.from_float(v)
    lo = np.nextafter(v, -np.inf) if fv > q else v
    hi = np.nextafter(v,  np.inf) if fv < q else v
    return lo, hi

def _iv_add(a, b):
    return _down(a[0] + b[0]), _up(a[1] + b[1])

def _iv_sub(a, b):
    return _down(a[0] - b[1]), _up(a[1] - b[0])

def _iv_mul(a, b):
    import numpy as np
    p = np.stack((a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]))
    return _down(np.min(p, axis=0)), _up(np.max(p, axis=0))

def _series_at_point(t):
    """Rigorous sin/cos intervals at nonnegative binary64 points t <= pi/4.

    Odd Taylor partial sums through degrees 19 and 20 bracket sin; the cosine
    partial sums through degrees 18 and 20 bracket cos.  Every elementary
    operation is rounded outwards.
    """
    import numpy as np
    t2 = (np.maximum(0.0, _down(t*t)), _up(t*t))

    # sin: S_9 (last term negative) is below sin; S_10 is above sin.
    sm = (t.copy(), t.copy())
    s = (t.copy(), t.copy())
    sin_lo = None
    for k in range(1, 11):
        d = float((2*k)*(2*k+1))
        sm = (_down(_down(sm[0]*t2[0])/d), _up(_up(sm[1]*t2[1])/d))
        s = _iv_sub(s, sm) if k & 1 else _iv_add(s, sm)
        if k == 9:
            sin_lo = s[0].copy()
    sin_hi = s[1]

    # cos: C_9 is below cos; C_10 is above cos.
    one = np.ones_like(t)
    cm = (one.copy(), one.copy())
    c = (one.copy(), one.copy())
    cos_lo = None
    for k in range(1, 11):
        d = float((2*k-1)*(2*k))
        cm = (_down(_down(cm[0]*t2[0])/d), _up(_up(cm[1]*t2[1])/d))
        c = _iv_sub(c, cm) if k & 1 else _iv_add(c, cm)
        if k == 9:
            cos_lo = c[0].copy()
    cos_hi = c[1]
    return (cos_lo, cos_hi), (sin_lo, sin_hi)

def _harmonic_sin_cos(odd, j, N, den, pilo_f, pihi_f):
    """Enclose cos(2*pi*(j/den)*(2i+1)/N) and sin of the same angle.

    In units of 2*pi/(den*N) the angle is the integer j*(2i+1), so the reduction to
    |angle| <= pi/4 is exact integer arithmetic.  Taylor evaluation is therefore stable
    and no platform trigonometric routine is used.
    """
    import numpy as np
    period = den*N
    quadrant = period//4
    half_quadrant = period//8
    num = (j*odd) % period
    k = ((num + half_quadrant)//quadrant) % 4
    rem = num - k*quadrant
    rem = np.where(rem > period//2, rem-period, rem)
    rem = np.where(rem < -period//2, rem+period, rem)
    sg = np.where(rem < 0, -1.0, 1.0)
    ar = np.abs(rem)
    assert np.all(ar <= half_quadrant)

    rv = ar.astype(np.float64) / float(period//2)
    rlo, rhi = np.maximum(0.0, _down(rv)), _up(rv)
    ylo = _down(rlo*pilo_f)
    yhi = _up(rhi*pihi_f)
    # sin is increasing and cos decreasing on [0,pi/4].
    cylo, sylo = _series_at_point(ylo)
    cyhi, syhi = _series_at_point(yhi)
    sy = (sylo[0], syhi[1])
    cy = (cyhi[0], cylo[1])

    # angle = k*pi/2 + sign(rem)*y.
    clo = np.empty_like(ylo); chi = np.empty_like(ylo)
    slo = np.empty_like(ylo); shi = np.empty_like(ylo)
    for q in range(4):
        m = k == q
        if q == 0:
            clo[m], chi[m] = cy[0][m], cy[1][m]
            slo[m] = np.where(sg[m] > 0, sy[0][m], -sy[1][m])
            shi[m] = np.where(sg[m] > 0, sy[1][m], -sy[0][m])
        elif q == 1:
            clo[m] = np.where(sg[m] > 0, -sy[1][m], sy[0][m])
            chi[m] = np.where(sg[m] > 0, -sy[0][m], sy[1][m])
            slo[m], shi[m] = cy[0][m], cy[1][m]
        elif q == 2:
            clo[m], chi[m] = -cy[1][m], -cy[0][m]
            slo[m] = np.where(sg[m] > 0, -sy[1][m], sy[0][m])
            shi[m] = np.where(sg[m] > 0, -sy[0][m], sy[1][m])
        else:
            clo[m] = np.where(sg[m] > 0, sy[0][m], -sy[1][m])
            chi[m] = np.where(sg[m] > 0, sy[1][m], -sy[0][m])
            slo[m], shi[m] = -cy[1][m], -cy[0][m]
    return (clo, chi), (slo, shi)

def _upward_pairwise_sum(a):
    """Upper bound for a sum of nonnegative binary64 values."""
    import numpy as np
    a = np.asarray(a, dtype=np.float64)
    while a.size > 1:
        n = a.size
        paired = _up(a[:n//2] + a[n//2:2*(n//2)])
        if n & 1:
            a = np.concatenate((paired, a[-1:]))
        else:
            a = paired
    return float(a[0]) if a.size else 0.0

def envelope_sum_interval_upper(C, Om, N, S2bounds, pilo, pihi, chunk=50000):
    """Rigorous upper bound for delta*sum h(x_i).

    Basic arithmetic is enclosed with np.nextafter.  No sin/cos implementation,
    BLAS reduction, or unverified blanket slack enters this calculation.
    """
    import numpy as np
    info = np.finfo(np.float64)
    assert np.dtype(np.float64).itemsize == 8 and info.bits == 64 and info.nmant == 52
    den = C['den']
    assert all((den*x).denominator == 1 for x in C['xis'])
    assert (den*N) % 8 == 0, "den*N must be divisible by 8 for exact octant reduction"

    lam = _qinterval(C['lam'])
    nu = _qinterval(C['nup'] - C['num'])
    kap = _qinterval(C['kapp'] - C['kapm'])
    Omiv = _qinterval(Om)

    # Ge has coefficients -(c_j + 4 gamma_j S_j^2); Go has -8 gamma_j tau_j.
    # A harmonic with c_j = gam_j = 0 contributes exactly zero to both G_e and G_o, so it can
    # be dropped with no effect on the result.  Certificates are sparse -- typically 3-30% of
    # harmonics carry a nonzero multiplier -- and the Taylor enclosure per harmonic per grid
    # point is the dominant cost, so this is the difference between minutes and hours.
    active = []
    for idx, (cj, gj, tj, (slo, shi)) in enumerate(
            zip(C['c'], C['gam'], C['tau'], S2bounds), start=1):
        if cj == 0 and gj == 0:
            continue
        ea = (np.float64(_qinterval(-(cj + 4*gj*shi))[0]),
              np.float64(_qinterval(-(cj + 4*gj*slo))[1]))
        oa = tuple(np.float64(v) for v in _qinterval(-8*gj*tj))
        active.append((idx, ea, oa))
    pilo_f = _qinterval(pilo)[0]
    pihi_f = _qinterval(pihi)[1]
    total_hi = 0.0
    om_num, om_den = Om.numerator, Om.denominator

    for s0 in range(0, N, chunk):
        ii = np.arange(s0, min(s0+chunk, N), dtype=np.int64)
        odd = 2*ii + 1

        # x_i=(2i+1)/N and theta=pi*x_i/8.  Enclose both without transcendental calls.
        xv = odd.astype(np.float64) / float(N)
        x = (_down(xv), _up(xv))
        trig_e = (np.zeros_like(xv), np.zeros_like(xv))
        trig_o = (np.zeros_like(xv), np.zeros_like(xv))
        for j, ea, oa in active:
            cost, sint = _harmonic_sin_cos(odd, j, N, den, pilo_f, pihi_f)
            trig_e = _iv_add(trig_e, _iv_mul(ea, cost))
            trig_o = _iv_add(trig_o, _iv_mul(oa, sint))
        Ge = _iv_add(lam, _iv_add(_iv_mul(nu, _iv_mul(x, x)), trig_e))
        Go = _iv_add(_iv_mul(kap, x), trig_o)
        al = _iv_mul(_iv_add(Ge, Go), (0.5, 0.5))
        be = _iv_mul(_iv_sub(Ge, Go), (0.5, 0.5))

        U2 = (_down(2.0-x[1]), _up(2.0-x[0]))
        A = (np.minimum(Omiv[0], U2[0]), np.minimum(Omiv[1], U2[1]))
        zero = np.zeros_like(xv)
        UmO = _iv_sub(U2, Omiv)
        pos = (np.maximum(zero, UmO[0]), np.maximum(zero, UmO[1]))
        B = (np.minimum(Omiv[0], pos[0]), np.minimum(Omiv[1], pos[1]))

        cand = [zero,
                _iv_mul(al, A)[1],
                _iv_mul(be, A)[1]]
        c4 = _iv_add(_iv_mul(al, Omiv), _iv_mul(be, B))[1]
        c5 = _iv_add(_iv_mul(al, B), _iv_mul(be, Omiv))[1]
        # Exact test of 2-x_i > Omega using int64 arithmetic.
        mk = (2*N-odd)*om_den > om_num*N
        cand.extend((np.where(mk, c4, -np.inf), np.where(mk, c5, -np.inf)))
        hhi = np.maximum.reduce(cand)
        chunk_hi = _upward_pairwise_sum(hhi)
        total_hi = float(_up(np.float64(total_hi) + np.float64(chunk_hi)))

    return F.from_float(total_hi) * F(2, N)

def verify_box(rec, pilo, pihi):
    # NB: rec['Omega'] is deliberately NOT read.  The level is re-derived below by bisecting the
    # rigorous criterion on this verifier's own grid, so nothing about the certified number is
    # taken from the stored file; the stored value is informational and may differ slightly.
    den = rec['xi_den']; m = rec['n_freq']
    C = dict(lam=F(*rec['lam']), nup=F(*rec['nup']), num=F(*rec['num']),
             kapp=F(*rec['kapp']), kapm=F(*rec['kapm']),
             c=[F(*z) for z in rec['c']], gam=[F(*z) for z in rec['gam']],
             tau=[F(*z) for z in rec['tau']],
             Elo=F(*rec['box'][0]), Ehi=F(*rec['box'][1]),
             xis=[F(j, den) for j in range(1, m+1)], den=den)
    for nm in ('nup','num','kapp','kapm'):
        assert C[nm] >= 0, f"sign condition violated: {nm} < 0"
    assert all(x >= 0 for x in C['c']),   "sign condition violated: some c_j < 0"
    assert all(x >= 0 for x in C['gam']), "sign condition violated: some gam_k < 0"
    b = [S2_bounds(j, den, pilo, pihi) for j in range(1, m+1)]
    S2hi = [t[1] for t in b]; S2mid = [(l+h)/2 for l, h in b]
    V = V_lower(C, S2hi)
    Ge, Go, Ged, God = sup_bounds(C, S2hi, pihi)
    return C, V, Ge, Go, Ged, God, b, S2mid

def excluded_fast(C, V, sups, S2mid, Om, N):
    """Non-rigorous search predicate; its proposed level is checked rigorously later."""
    Ge, Go, Ged, God = sups
    L = Om*(Ged+God) + (Ge+Go)
    delta = F(2, N)
    ub = 2*F(envelope_sum(C, Om, N, S2mid)) + L*delta
    return V > ub, ub, L

def excluded_interval(C, V, sups, S2bounds, Om, N, pilo, pihi):
    """Rigorous exclusion predicate using interval evaluation of every midpoint."""
    Ge, Go, Ged, God = sups
    L = Om*(Ged+God) + (Ge+Go)
    delta = F(2, N)
    grid_ub = envelope_sum_interval_upper(C, Om, N, S2bounds, pilo, pihi)
    ub = 2*grid_ub + L*delta
    return V > ub, ub, L

def certified_level(rec, N_search, N_final, pilo, pihi, verbose=True, digits=6):
    C, V, Ge, Go, Ged, God, S2bounds, S2mid = verify_box(rec, pilo, pihi)
    sups = (Ge, Go, Ged, God)
    lo, hi = F(30,100), F(40,100)
    for _ in range(24):                       # fast proposal only; final check is rigorous
        mid = (lo+hi)/2
        ok, _, _ = excluded_fast(C, V, sups, S2mid, mid, N_search)
        if ok: lo = mid
        else:  hi = mid
    # Round the level down to `digits` decimals.  Coarser rounding reports a slightly smaller
    # bound but leaves proportionally more margin in the rigorous check, which matters once the
    # search grid equals the final grid and the bisection pushes right to the edge.
    scale = 10**digits
    Om = F(int(lo*scale), scale)
    ok, ub, L = excluded_interval(C, V, sups, S2bounds, Om, N_final, pilo, pihi)
    if verbose:
        print(f"  E in [{float(C['Elo']):.5f},{float(C['Ehi']):.5f}]  certified Omega = {float(Om):.6f}"
              f"   V={float(V):.8f}  rigorous 2*int h <= {float(ub):.8f}"
              f"   margin {float(V-ub):+.2e}   {'OK' if ok else 'FAIL'}   (L={float(L):.1f})")
    return ok, Om, V, ub

def main(path, N=2000000, Ns=None, digits=5):
    # The bisection in certified_level charges the Lipschitz term L*delta at the SEARCH grid.
    # If that grid is coarser than the one which finally has to certify, the proposed level is
    # pessimistic by L*(2/Ns - 2/N) for reasons having nothing to do with the certificate, and
    # that pessimism caps the reported bound.  Default the search grid to the final grid.
    if Ns is None: Ns = N
    D = json.load(open(path))
    pilo, pihi = pi_enclosure()
    print(f"pi enclosed to width {float(pihi-pilo):.2e} by Machin's formula (checked internally)")
    print(f"sqrt(2) enclosure checked by squaring")
    s2_selftest(pilo, pihi)
    print(f"S^2 Taylor enclosures agree with closed-form sqrt(2) values (den=16 self-test)")
    print(f"search grid N = {Ns}, final verification grid N = {N}\n")
    ok_all = True; best = None; ivals = []
    for rec in D['boxes']:
        ok, Om, V, ub = certified_level(rec, Ns, N, pilo, pihi, digits=digits)
        ok_all &= ok
        best = Om if best is None else min(best, Om)
        ivals.append((F(*rec['box'][0]), F(*rec['box'][1])))
    # The boxes must actually TILE [0, Etop] -- start at 0 and leave no gap.  Checking only that
    # some box reaches past Emax would accept e.g. [0,0.01] u [0.2,0.2152], which excludes nothing
    # on (0.01, 0.2).  Overlaps are harmless; gaps are not.
    ivals.sort()
    tiles = bool(ivals) and ivals[0][0] == 0
    Etop = ivals[0][1] if ivals else F(0)
    for lo, hi in ivals[1:]:
        if lo > Etop:                     # gap between the covered prefix and this box
            tiles = False
            break
        Etop = max(Etop, hi)
    Emax_needed_sq = F(1,3) - F(1, 24)/(best*best)
    covers = tiles and Emax_needed_sq <= Etop*Etop
    print(f"\n  all boxes excluded: {ok_all}")
    print(f"  E-range: any f with ||M||_inf <= {float(best):.6f} has E^2 <= {float(Emax_needed_sq):.8f},")
    print(f"           i.e. |E| <= {float(Emax_needed_sq)**0.5:.6f}; boxes tile [0,{float(Etop):.6f}]"
          f" with no gap: {tiles}   {'OK' if covers else 'GAP'}")
    if ok_all and covers:
        print(f"\n  VERIFIED:  mu >= {float(best):.6f}")
        return 0
    print("\n  VERIFICATION FAILED")
    return 1

if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "data/certificates.json",
                  int(sys.argv[2]) if len(sys.argv) > 2 else 2000000,
                  int(sys.argv[3]) if len(sys.argv) > 3 else None,
                  int(sys.argv[4]) if len(sys.argv) > 4 else 5))
