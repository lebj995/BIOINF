#!/usr/bin/env python3
"""
PROCONSUL: probabilistic extension of DIAMOnD-like connectivity significance.
Implements the pseudocode in the provided slides:
- compute connectivity p-values for boundary nodes (DIAMOnD idea)
- transform with Shannon information: -log(p)
- softmax with temperature T to obtain a sampling distribution
- sample ONE gene per iteration, add to seeds, repeat until m genes
- repeat for n_runs, score genes by position (m..1), average, output final ranking
"""

import argparse
import math
import random
from collections import defaultdict

# -----------------------------
# Math helpers: log-combinations and hypergeometric survival
# -----------------------------
def logC(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    # log(n choose k) via log-gamma
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def hypergeom_sf(k_obs_minus_1: int, N: int, K: int, n: int) -> float:
    """
    Survival function: P[X >= k_obs] where X ~ Hypergeom(N, K, n)
    We compute sum_{x=k_obs..min(n,K)} pmf(x) in log-space for stability.
    """
    k_obs = k_obs_minus_1 + 1
    x_min = max(0, n - (N - K))
    x_max = min(n, K)
    if k_obs <= x_min:
        return 1.0
    if k_obs > x_max:
        return 0.0

    # log denominator: log C(N, n)
    log_den = logC(N, n)

    # log-sum-exp
    logs = []
    for x in range(k_obs, x_max + 1):
        log_pmf = logC(K, x) + logC(N - K, n - x) - log_den
        logs.append(log_pmf)

    m = max(logs)
    s = sum(math.exp(v - m) for v in logs)
    return math.exp(m) * s

# -----------------------------
# Graph I/O
# -----------------------------
def read_edge_list(path: str):
    adj = defaultdict(set)
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            if a == b:
                continue
            adj[a].add(b)
            adj[b].add(a)
    return adj

def read_seeds(path: str):
    seeds = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                seeds.append(g)
    # unique keep order
    seen = set()
    out = []
    for g in seeds:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out

# -----------------------------
# PROCONSUL core
# -----------------------------
def softmax_temperature(values, T: float):
    # values are real numbers (here: Shannon info = -log(p))
    # softmax(values / T)
    if T <= 0:
        raise ValueError("Temperature must be > 0.")
    scaled = [v / T for v in values]
    m = max(scaled)
    exps = [math.exp(v - m) for v in scaled]
    Z = sum(exps)
    if Z == 0 or not math.isfinite(Z):
        # fallback: uniform
        return [1.0 / len(values)] * len(values)
    return [e / Z for e in exps]

def one_round(adj, seeds_init, m_pred: int, T: float, rng: random.Random):
    nodes = list(adj.keys())
    node_set = set(nodes)
    S = set(g for g in seeds_init if g in node_set)
    N = len(node_set)
    predicted = []

    # if no seeds in network, impossible
    if len(S) == 0:
        return predicted

    while len(predicted) < m_pred:
        # boundary candidates: neighbors of S not in S
        cand = set()
        for s in S:
            cand.update(adj[s])
        cand.difference_update(S)

        if not cand:
            break

        s_size = len(S)

        # compute p-values for candidates
        cand_list = []
        pvals = []
        infos = []
        for v in cand:
            k = len(adj[v])
            if k == 0:
                continue
            ks = sum((u in S) for u in adj[v])
            if ks == 0:
                continue

            # DIAMOnD connectivity p-value: P[X >= ks], X ~ Hypergeom(N, s_size, k)
            p = hypergeom_sf(ks - 1, N=N, K=s_size, n=k)

            # guard against 0
            if p <= 0.0:
                p = 1e-300
            info = -math.log(p)

            cand_list.append(v)
            pvals.append(p)
            infos.append(info)

        if not cand_list:
            break

        probs = softmax_temperature(infos, T=T)

        # sample one gene
        # (convert to cumulative for speed/compat)
        r = rng.random()
        cum = 0.0
        chosen = cand_list[-1]
        for g, pr in zip(cand_list, probs):
            cum += pr
            if r <= cum:
                chosen = g
                break

        predicted.append(chosen)
        S.add(chosen)

    return predicted

def proconsul(adj, seeds_init, m_pred: int, n_runs: int, T: float, seed: int = 123):
    rng = random.Random(seed)

    # accumulate scores: position-based (top gets m, then m-1, ..., 1)
    score_sum = defaultdict(float)
    count_in_runs = defaultdict(int)

    for run in range(n_runs):
        # new RNG stream per run for reproducibility but different draws
        run_seed = rng.randint(1, 10**9)
        run_rng = random.Random(run_seed)

        pred = one_round(adj, seeds_init, m_pred=m_pred, T=T, rng=run_rng)

        # assign position scores
        # if shorter than m_pred, score based on produced length
        L = len(pred)
        for i, g in enumerate(pred):
            # higher rank -> higher score
            score = float(L - i)
            score_sum[g] += score
            count_in_runs[g] += 1

    # average score among runs (as in slides: average ranks/scores)
    results = []
    for g, s in score_sum.items():
        avg = s / max(1, n_runs)
        results.append((g, avg, count_in_runs[g]))

    # sort by avg score desc, then by frequency desc, then gene
    results.sort(key=lambda x: (-x[1], -x[2], x[0]))
    return results

# -----------------------------
# CLI
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("netfile", help="Edge list file (2 columns, whitespace-separated)")
    ap.add_argument("seedfile", help="Seed genes file (1 per line)")
    ap.add_argument("--m", type=int, default=100, help="Number of predicted genes per run (m)")
    ap.add_argument("--n_runs", type=int, default=50, help="Number of stochastic runs (n)")
    ap.add_argument("--temp", type=float, default=1.0, help="Softmax temperature T")
    ap.add_argument("--out", required=True, help="Output path (tsv)")
    ap.add_argument("--seed", type=int, default=123, help="Random seed (reproducibility)")
    args = ap.parse_args()

    adj = read_edge_list(args.netfile)
    seeds = read_seeds(args.seedfile)

    res = proconsul(
        adj=adj,
        seeds_init=seeds,
        m_pred=args.m,
        n_runs=args.n_runs,
        T=args.temp,
        seed=args.seed,
    )

    # write output
    with open(args.out, "w", encoding="utf-8") as f:
        f.write("rank\tgene\tavg_score\tcount_in_runs\n")
        for i, (g, avg, cnt) in enumerate(res, start=1):
            f.write(f"{i}\t{g}\t{avg:.6f}\t{cnt}\n")

if __name__ == "__main__":
    main()
