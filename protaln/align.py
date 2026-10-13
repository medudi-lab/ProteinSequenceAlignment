"""Semi-global protein alignment (BLOSUM62) against a reference.

End gaps are free so truncated sequences are not penalised in the score; the
truncation itself is reported by the completeness check.
"""

_AA = "ARNDCQEGHILKMFPSTWYV"
_ROWS = """
4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4
"""
BLOSUM62 = {}
for _a, _row in zip(_AA, _ROWS.strip().splitlines()):
    for _b, _v in zip(_AA, _row.split()):
        BLOSUM62[(_a, _b)] = int(_v)

GAP = -6


def score(a, b):
    return BLOSUM62.get((a, b), -1)


def align(ref, query, gap=GAP):
    """Align query to ref. Returns (aligned_ref, aligned_query, score)."""
    n, m = len(ref), len(query)
    # S[i][j]: best score aligning ref[:i] with query[:j]; free leading/trailing
    # gaps at the ends of the *query* (ref end gaps free too, for local-like fit).
    S = [[0] * (m + 1) for _ in range(n + 1)]
    T = [[0] * (m + 1) for _ in range(n + 1)]  # 0 diag, 1 up (gap in query), 2 left (gap in ref)
    for i in range(1, n + 1):
        T[i][0] = 1
    for j in range(1, m + 1):
        T[0][j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = S[i - 1][j - 1] + score(ref[i - 1], query[j - 1])
            u = S[i - 1][j] + (0 if j == m else gap)
            l = S[i][j - 1] + (0 if i == n else gap)
            best = max(d, u, l)
            S[i][j] = best
            T[i][j] = 0 if best == d else (1 if best == u else 2)
    i, j = n, m
    ra, qa = [], []
    while i > 0 or j > 0:
        t = T[i][j]
        if t == 0 and i > 0 and j > 0:
            ra.append(ref[i - 1]); qa.append(query[j - 1]); i -= 1; j -= 1
        elif t == 1 and i > 0:
            ra.append(ref[i - 1]); qa.append("-"); i -= 1
        else:
            ra.append("-"); qa.append(query[j - 1]); j -= 1
    return "".join(reversed(ra)), "".join(reversed(qa)), S[n][m]


def project(aligned_ref, aligned_query):
    """Project onto reference coordinates.

    Returns (residues, insertions): residues[i] is the query residue at reference
    position i ('-' if deleted); insertions maps reference index -> inserted text
    that follows that position (-1 means before the first residue).
    """
    residues, insertions, pos = [], {}, -1
    for r, q in zip(aligned_ref, aligned_query):
        if r == "-":
            insertions[pos] = insertions.get(pos, "") + q
        else:
            pos += 1
            residues.append(q)
    return residues, insertions
