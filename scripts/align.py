"""Pure-python Gotoh alignment (affine gaps) for sequences and profiles; local or semi-global (free end gaps)."""
import numpy as np
from common import B62, AAIDX

NEG = -1e9

def gotoh(S, go, ge, mode="local"):
    """S: (n x m) score matrix. Gap of length k costs go + (k-1)*ge.  mode: 'local' | 'semi' (free end gaps).
    Returns (score, pairs) with pairs = list of (i or None, j or None) in order."""
    n, m = S.shape
    M = np.full((n+1, m+1), NEG); X = np.full((n+1, m+1), NEG); Y = np.full((n+1, m+1), NEG)
    M[0, :] = 0.0; M[:, 0] = 0.0   # free starts (local and semi-global)
    Sl = S.tolist()
    Ml = M.tolist(); Xl = X.tolist(); Yl = Y.tolist()
    best, bi, bj, bs = 0.0, 0, 0, "M"
    local = (mode == "local")
    for i in range(1, n+1):
        Mi, Xi, Yi = Ml[i], Xl[i], Yl[i]; Mp, Xp, Yp = Ml[i-1], Xl[i-1], Yl[i-1]; Si = Sl[i-1]
        for j in range(1, m+1):
            prev = Mp[j-1]
            x = Xp[j-1]; y = Yp[j-1]
            if x > prev: prev = x
            if y > prev: prev = y
            if local and prev < 0: prev = 0.0
            mv = prev + Si[j-1]
            xv = max(Mp[j] - go, Xp[j] - ge, Yp[j] - go)
            yv = max(Mi[j-1] - go, Yi[j-1] - ge, Xi[j-1] - go)
            if local and mv < 0: mv = 0.0
            Mi[j] = mv; Xi[j] = xv; Yi[j] = yv
            if local:
                if mv > best: best, bi, bj, bs = mv, i, j, "M"
    if not local:
        best = NEG
        for j in range(0, m+1):
            for st, mat in (("M", Ml), ("X", Xl), ("Y", Yl)):
                if mat[n][j] > best: best, bi, bj, bs = mat[n][j], n, j, st
        for i in range(0, n+1):
            for st, mat in (("M", Ml), ("X", Xl), ("Y", Yl)):
                if mat[i][m] > best: best, bi, bj, bs = mat[i][m], i, m, st
    # traceback
    pairs = []
    i, j, st = bi, bj, bs
    eps = 1e-7
    if not local:
        # trailing unaligned tail
        for k in range(n, bi, -1): pairs.append((k-1, None))
        for k in range(m, bj, -1): pairs.append((None, k-1))
    while i > 0 and j > 0:
        if st == "M":
            v = Ml[i][j]
            if local and v <= eps: break
            pairs.append((i-1, j-1))
            prev = v - Sl[i-1][j-1]
            if not local and (i == 1 or j == 1) and abs(prev) < eps and (Ml[i-1][j-1] == 0.0):
                i -= 1; j -= 1; break
            if abs(Ml[i-1][j-1] - prev) < eps: st = "M"
            elif abs(Xl[i-1][j-1] - prev) < eps: st = "X"
            elif abs(Yl[i-1][j-1] - prev) < eps: st = "Y"
            else:
                i -= 1; j -= 1; break
            i -= 1; j -= 1
        elif st == "X":
            v = Xl[i][j]; pairs.append((i-1, None))
            if abs(Ml[i-1][j] - go - v) < eps: st = "M"
            elif abs(Xl[i-1][j] - ge - v) < eps: st = "X"
            else: st = "Y"
            i -= 1
        else:
            v = Yl[i][j]; pairs.append((None, j-1))
            if abs(Ml[i][j-1] - go - v) < eps: st = "M"
            elif abs(Yl[i][j-1] - ge - v) < eps: st = "Y"
            else: st = "X"
            j -= 1
    if not local:
        for k in range(i, 0, -1): pairs.append((k-1, None))
        for k in range(j, 0, -1): pairs.append((None, k-1))
    pairs.reverse()
    return best, pairs

def score_matrix(a, b):
    ia = [AAIDX.get(c, -1) for c in a]; ib = [AAIDX.get(c, -1) for c in b]
    S = np.full((len(a), len(b)), -1.0)
    for i, x in enumerate(ia):
        if x < 0: continue
        for j, y in enumerate(ib):
            if y >= 0: S[i, j] = B62[x, y]
    return S

def pair_align(a, b, mode="local", go=11, ge=1):
    sc, pairs = gotoh(score_matrix(a, b), go, ge, mode)
    ra = "".join(a[i] if i is not None else "-" for i, j in pairs)
    rb = "".join(b[j] if j is not None else "-" for i, j in pairs)
    return sc, ra, rb, pairs

def identity(ra, rb):
    cols = [(x, y) for x, y in zip(ra, rb) if x != "-" and y != "-"]
    if not cols: return 0.0
    return sum(x == y for x, y in cols) / len(cols)

def map_positions(a, b, pairs):
    """dict: index in a (0-based) -> index in b for aligned (matched) columns"""
    return {i: j for i, j in pairs if i is not None and j is not None}
