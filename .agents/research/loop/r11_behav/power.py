"""Exact power of a two-sided Fisher exact test (alpha 0.05) by enumeration."""

import itertools
import sys

from scipy.stats import binom, fisher_exact


def power(n1, n2, p1, p2, alpha=0.05):
    total = 0.0
    for a, b in itertools.product(range(n1 + 1), range(n2 + 1)):
        _odds, p = fisher_exact([[a, n1 - a], [b, n2 - b]], alternative="two-sided")
        if p <= alpha:
            total += binom.pmf(a, n1, p1) * binom.pmf(b, n2, p2)
    return total


ns = [int(x) for x in sys.argv[1:]] or [6, 12, 16, 20]
pairs = [(1.0, 0.5), (0.95, 0.5), (0.9, 0.5), (0.9, 0.3), (0.8, 0.3), (0.7, 0.2), (0.5, 0.1), (0.3, 0.0), (0.4, 0.05)]
print("p1 vs p2 | " + " | ".join(f"N={n}" for n in ns))
for p1, p2 in pairs:
    print(f"{p1:.2f} vs {p2:.2f} | " + " | ".join(f"{power(n, n, p1, p2):.2f}" for n in ns))
# the smallest detectable split at N=12 when one arm is at ceiling (12/12)
for k in range(12, -1, -1):
    _o, p = fisher_exact([[12, 0], [k, 12 - k]])
    if p <= 0.05:
        print("N=12: 12/12 vs", k, "/12 is the largest count still significant, p =", round(p, 4))
        break
for k in range(0, 13):
    _o, p = fisher_exact([[0, 12], [k, 12 - k]])
    if p <= 0.05:
        print("N=12: 0/12 vs", k, "/12 is the smallest count significant, p =", round(p, 4))
        break
