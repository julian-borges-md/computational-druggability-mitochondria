#!/usr/bin/env python3
"""Spearman correlation, its p value, and the normal quantile, using nothing
but the Python standard library.

preflight.py must run anywhere, including a machine with no scientific stack.
A verification script that cannot be run is worse than no verification script,
because it gets skipped and then believed. Values match scipy.stats.spearmanr
and scipy.stats.norm.ppf to better than 1e-10 on this data.
"""
import math

def _ranks(xs):
    """Average ranks, so ties are handled the way Spearman requires."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r

def _pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return float('nan')
    return sxy / math.sqrt(sxx * syy)

def _betacf(a, b, x, itmax=300, eps=3e-16):
    """Continued fraction for the incomplete beta function, modified Lentz."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < tiny: d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny: d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny: c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny: d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny: c = tiny
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h

def betainc(a, b, x):
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log1p(-x))
    front = math.exp(lbeta)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                          + b * math.log1p(-x) + a * math.log(x)) * _betacf(b, a, 1.0 - x) / b

def t_sf(t, df):
    """Upper tail of Student's t."""
    t = abs(t)
    return 0.5 * betainc(df / 2.0, 0.5, df / (df + t * t))

def spearmanr(x, y):
    """Returns (rho, two sided p). Same definition scipy uses: Pearson on
    average ranks, p from the t approximation with n - 2 degrees of freedom."""
    n = len(x)
    rho = _pearson(_ranks(list(x)), _ranks(list(y)))
    if not (-1.0 < rho < 1.0):
        return rho, 0.0
    t = rho * math.sqrt((n - 2) / ((1.0 + rho) * (1.0 - rho)))
    return rho, 2.0 * t_sf(t, n - 2)

def norm_ppf(p):
    """Inverse standard normal CDF, by bisection on erf. Exact to 1e-14."""
    lo, hi = -40.0, 40.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if 0.5 * (1.0 + math.erf(mid / math.sqrt(2.0))) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0

def fisher_ci(rho, n, conf=0.95):
    """Fisher z confidence interval for a correlation."""
    z = math.atanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    q = norm_ppf(0.5 + conf / 2.0)
    return math.tanh(z - q * se), math.tanh(z + q * se)

def mdes(n, alpha=0.05, power=0.80):
    """Minimum detectable |rho| at this sample size."""
    return math.tanh((norm_ppf(1 - alpha / 2) + norm_ppf(power)) / math.sqrt(n - 3))
