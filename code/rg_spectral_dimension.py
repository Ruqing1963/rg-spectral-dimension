"""
Renormalization-group flow of the spectral dimension: bond percolation on the diamond hierarchical lattice (DHL).

Part A  Exact RG map of the edge-retention probability, p' = 1 - (1 - p^2)^2 (two parallel pairs of edges in
        series).  Fixed points 0 (fragmented), 1 (intact DHL, d_s = 2) and p_c = (sqrt 5 - 1)/2 (unstable);
        thermal eigenvalue lambda = f'(p_c) = 4 p_c^2 = 6 - 2 sqrt 5, nu = ln 2 / ln lambda.
Part B  Mass of the cluster attached to a pole at p_c.  The vector (E[M 1_c], E[M 1_{not c}]) obeys an exact
        linear recursion; its leading eigenvalue Lambda gives the growth exponent theta = log2 Lambda and, with
        resistance growing linearly in the length, d_s = 2 theta / (theta + 1).  Checked by a population-
        dynamics simulation that samples the connection indicator exactly at p_c.
Part C  Direct check: heat-kernel traces of sampled spanning clusters of the generation-n DHL at p_c, and of the
        intact DHL; d_s(tau) = 2 tau <lambda>_tau.  Off-critical samples show the flow to d_s = 0 or 2.
Part D  Rewiring is relevant: Newman-Watts shortcuts on a ring destroy the d_s = 1 plateau beyond tau ~ xi^2.

Usage:  python rg_spectral_dimension.py [--show]
Figures go to ./figures (PNG, PDF), data to ./data (CSV); override with RG_FIG_DIR, RG_DATA_DIR.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

FIG_DIR = os.environ.get("RG_FIG_DIR", "figures")
DATA_DIR = os.environ.get("RG_DATA_DIR", "data")
rng = np.random.default_rng(20261007)
PC = (np.sqrt(5) - 1) / 2


def savefig(fig, name):
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=150)


def savedata(name, columns, cols, meta=()):
    os.makedirs(DATA_DIR, exist_ok=True)
    arr = np.column_stack([np.asarray(c, dtype=float) for c in cols])
    with open(os.path.join(DATA_DIR, f"{name}.csv"), "w", encoding="utf-8", newline="\n") as fh:
        for line in meta:
            fh.write(f"# {line}\n")
        fh.write(",".join(columns) + "\n")
        np.savetxt(fh, arr, delimiter=",", fmt="%.10g")


# ================================================================ Part A: RG map
def f(p):
    return 1 - (1 - p * p) ** 2


def part_a():
    lam = 4 * PC ** 2
    nu = np.log(2) / np.log(lam)
    print(f"[A] RG map p' = 1-(1-p^2)^2: p_c = (sqrt5-1)/2 = {PC:.12f}, f(p_c) - p_c = {f(PC) - PC:.1e}, "
          f"lambda = f'(p_c) = 4 p_c^2 = {lam:.12f} (= 6 - 2 sqrt5 = {6 - 2 * np.sqrt(5):.12f}), nu = ln2/ln lambda = {nu:.10f}")
    starts = [PC - 0.05, PC - 0.01, PC - 0.001, PC + 0.001, PC + 0.01, PC + 0.05]
    flows = []
    for p in starts:
        traj = [p]
        for _ in range(16):
            traj.append(f(traj[-1]))
        flows.append(traj)
    savedata("rg_flow_trajectories", ["level"] + [f"p0_{s:.4f}" for s in starts],
             [np.arange(17)] + flows, ["iterates of p' = 1-(1-p^2)^2 from starting points near p_c"])
    return lam, nu, starts, flows


# ================================================================ Part B: mass growth and spectral dimension
def mass_matrix(P):
    """(a, b) = (E[M 1_c], E[M 1_{not c}]) -> (a', b') for one RG step at connection probability P."""
    A = 4 * P + 4 * P * P * (1 - P)
    B = 4 * P * P * (1 + P)
    C = 2 * (1 - P * P) * (1 - P)
    D = 2 * (1 - P * P) * (1 + P)
    return np.array([[A, B], [C, D]])


def part_b(levels=20, pool=400000):
    Lam_closed = 7 - 2 * np.sqrt(5) + np.sqrt(73 - 32 * np.sqrt(5))
    Lam = np.max(np.linalg.eigvals(mass_matrix(PC)).real)
    theta = np.log2(Lam)
    ds = 2 * theta / (theta + 1)
    print(f"[B] mass recursion at p_c: Lambda = {Lam:.12f} (closed form 7-2sqrt5+sqrt(73-32sqrt5) = {Lam_closed:.12f}), "
          f"d_f = log2 Lambda = {theta:.10f} (2 d_f/(d_f+1) = {ds:.4f} would assume zeta = 1; see below)")
    # exact iteration of (a, b) from a single edge
    ab = np.array([PC, 0.0])
    exact = [ab[0] / PC]
    for _ in range(levels):
        ab = mass_matrix(PC) @ ab
        exact.append(ab[0] / PC)
    # population dynamics: conditional pools for connected / disconnected sub-networks; the indicator is drawn
    # exactly with probability p_c at every level, so the unstable fixed point causes no drift
    con = {"M1": np.ones(pool), "M2": np.ones(pool), "R": np.ones(pool)}       # an open edge
    dis = {"M1": np.zeros(pool), "M2": np.zeros(pool)}                         # a closed edge
    stats = []

    def draw(k):
        c = rng.random(k) < PC
        i = rng.integers(0, pool, k)
        M1 = np.where(c, con["M1"][i], dis["M1"][i])
        M2 = np.where(c, con["M2"][i], dis["M2"][i])
        R = np.where(c, con["R"][i], np.inf)
        return c, M1, M2, R

    for n in range(1, levels + 1):
        K = 6 * pool                                                          # oversample, then split by outcome
        parts = [draw(K) for _ in range(4)]
        (c1, a1, b1, r1), (c2, a2, b2, r2), (c3, a3, b3, r3), (c4, a4, b4, r4) = parts

        def series(cx, mx1, mx2, rx, cy, my1, my2, ry):
            c = cx & cy
            M1 = mx1 + np.where(cx, my1, 0.0)
            M2 = my2 + np.where(cy, mx2, 0.0)
            return c, M1, M2, rx + ry

        sA = series(c1, a1, b1, r1, c2, a2, b2, r2)
        sB = series(c3, a3, b3, r3, c4, a4, b4, r4)
        c = sA[0] | sB[0]
        tot = (np.where(sA[0], sA[1], sA[1] + sA[2]) + np.where(sB[0], sB[1], sB[1] + sB[2]))
        M1 = np.where(c, tot, sA[1] + sB[1])
        M2 = np.where(c, tot, sA[2] + sB[2])
        with np.errstate(divide="ignore"):
            R = 1.0 / (1.0 / sA[3] + 1.0 / sB[3])
        pc_emp = c.mean()
        ci, di = np.flatnonzero(c)[:pool], np.flatnonzero(~c)[:pool]
        assert ci.size == pool and di.size == pool
        con = {"M1": M1[ci], "M2": M2[ci], "R": R[ci]}
        dis = {"M1": M1[di], "M2": M2[di]}
        stats.append((n, pc_emp, con["M1"].mean(), exact[n], np.median(con["R"]) / 2 ** n,
                      np.mean(np.log(con["R"])) / np.log(2) / n))
    arr = np.array(stats)
    grow = arr[-1, 2] / arr[-2, 2]
    # resistance growth exponent zeta = log2 rho from the stationary increments of E[log2 R | c] (last 8 levels)
    mlog = arr[:, 5] * arr[:, 0]
    inc = np.diff(mlog)[-8:]
    zeta, zeta_err = inc.mean(), inc.std(ddof=1) / np.sqrt(len(inc))
    ds_einstein = 2 * theta / (theta + zeta)
    print(f"    population dynamics ({pool} per pool, {levels} levels): P(connect) = {arr[-1, 1]:.4f}; "
          f"E[M | c] exact {arr[-1, 3]:.4e} vs sampled {arr[-1, 2]:.4e}; last growth ratio {grow:.4f} (Lambda {Lam:.4f})")
    print(f"    resistance of the spanning cluster: median R_n/2^n = " +
          ", ".join(f"{s[4]:.3f}" for s in stats[-5:]) + f" (levels {levels - 4}..{levels}); "
          f"E[log2 R | c] increments (last 8 levels) give zeta = {zeta:.4f} +- {zeta_err:.4f}")
    print(f"    Einstein relation for the full pole cluster: d_s = 2 theta/(theta + zeta) = {ds_einstein:.4f}")
    savedata("dhl_critical_cluster", ["level", "P_connect", "mean_mass_given_c_sampled", "mean_mass_given_c_exact",
                                      "median_R_over_2n", "mean_log2R_over_n"],
             [arr[:, i] for i in range(6)], [f"bond percolation on the DHL at p_c={PC:.12f}; pool size {pool}"])
    return Lam, theta, zeta, ds_einstein, arr


# ================================================================ Part C: heat kernel of sampled clusters
def dhl_edges(n):
    """Edges of the generation-n DHL (poles 0 and 1)."""
    edges = np.array([[0, 1]])
    nxt = 2
    for _ in range(n):
        m = len(edges)
        u, v = edges[:, 0], edges[:, 1]
        a = np.arange(nxt, nxt + m)
        b = np.arange(nxt + m, nxt + 2 * m)
        nxt += 2 * m
        edges = np.concatenate([np.stack([u, a], 1), np.stack([a, v], 1), np.stack([u, b], 1), np.stack([b, v], 1)])
    return edges, nxt


def pole_cluster_spectrum(edges, nv, p, require_span=True):
    keep = rng.random(len(edges)) < p
    e = edges[keep]
    if len(e) == 0:
        return None
    A = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(nv, nv))
    ncomp, lab = connected_components(A + A.T, directed=False)
    if require_span and lab[0] != lab[1]:
        return None
    nodes = np.flatnonzero(lab == lab[0])
    sub = e[np.isin(e[:, 0], nodes)]
    idx = -np.ones(nv, int)
    idx[nodes] = np.arange(len(nodes))
    k = len(nodes)
    L = np.zeros((k, k))
    i, j = idx[sub[:, 0]], idx[sub[:, 1]]
    np.add.at(L, (i, j), -1.0)
    np.add.at(L, (j, i), -1.0)
    np.add.at(L, (i, i), 1.0)
    np.add.at(L, (j, j), 1.0)
    return np.linalg.eigvalsh(L)


def ds_curve(spectra, tau):
    """Annealed: <P(tau)> and d_s = -2 dln<P>/dln tau.  Quenched: mean over samples of 2 tau <lambda>_tau."""
    P = np.zeros_like(tau)
    Q = np.zeros_like(tau)
    for lam in spectra:
        w = np.exp(-np.outer(tau, np.clip(lam, 0, None)))
        P += w.mean(1)
        Q += 2 * tau * (w * lam).sum(1) / w.sum(1)
    P /= len(spectra)
    Q /= len(spectra)
    return P, -2 * np.gradient(np.log(P), np.log(tau)), Q


def part_c(period, n=6, samples=60, n_big=7, samples_big=24):
    """period = time-scale factor per RG level (Lambda * 2^zeta); averages are taken over two full periods."""
    edges, nv = dhl_edges(n)
    big = dhl_edges(n_big)
    tau = np.logspace(-1, 5, 320)
    out = {}
    cases = [("intact DHL (p = 1)", 1.0, 1, (edges, nv)), ("p_c, spanning clusters", PC, samples, (edges, nv)),
             ("intact DHL, generation 7", 1.0, 1, big), (f"p_c, generation {n_big}", PC, samples_big, big),
             ("p = p_c + 0.08", PC + 0.08, samples // 2, (edges, nv)), ("p = p_c - 0.08", PC - 0.08, samples // 2, (edges, nv))]
    for label, p, ns, (edges_, nv_) in cases:
        spectra, tries = [], 0
        while len(spectra) < ns and tries < 200 * ns:
            tries += 1
            s = pole_cluster_spectrum(edges_, nv_, p, require_span=(p < 1 and abs(p - PC) < 0.01))
            if s is not None and len(s) > 1:
                spectra.append(s)
        P, ds, dq = ds_curve(spectra, tau)
        out[label] = (P, ds, np.mean([len(s) for s in spectra]), dq)
    t0 = 5.0
    win = (tau >= t0) & (tau <= t0 * period ** 2)
    print(f"[C] DHL heat kernels; d_s(tau) = -2 dln<P>/dln tau averaged in ln tau over two periods, "
          f"{t0} <= tau <= {t0 * period ** 2:.0f} (period factor {period:.3f}):")
    wavg = lambda y: np.trapz(y[win], np.log(tau[win])) / np.log(tau[win][-1] / tau[win][0])
    for label, (P, ds, size, dq) in out.items():
        print(f"    {label:>24}: mean cluster size {size:8.1f}, annealed <d_s> = {wavg(ds):.4f} (range [{ds[win].min():.3f}, "
              f"{ds[win].max():.3f}]), quenched <d_s> = {wavg(dq):.4f}")
    savedata("dhl_heat_kernel", ["tau"] + [f"ds_{k.split(' (')[0].replace(' ', '_').replace(',', '')}" for k in out],
             [tau] + [v[1] for v in out.values()], [f"generation-{n} DHL, annealed heat-kernel trace per vertex"])
    return tau, out, (t0, t0 * period ** 2)


# ================================================================ Part D: rewiring (Newman-Watts shortcuts)
def part_d(N=2000, phis=(0.0, 0.002, 0.01, 0.05)):
    tau = np.logspace(-1, 6, 300)
    curves = {}
    for phi in phis:
        i = np.arange(N)
        e = [np.stack([i, (i + 1) % N], 1)]
        k = rng.binomial(N, phi)
        if k:
            a, b = rng.integers(0, N, k), rng.integers(0, N, k)
            ok = a != b
            e.append(np.stack([a[ok], b[ok]], 1))
        e = np.concatenate(e)
        L = np.zeros((N, N))
        np.add.at(L, (e[:, 0], e[:, 1]), -1.0)
        np.add.at(L, (e[:, 1], e[:, 0]), -1.0)
        np.add.at(L, (e[:, 0], e[:, 0]), 1.0)
        np.add.at(L, (e[:, 1], e[:, 1]), 1.0)
        lam = np.linalg.eigvalsh(L)
        P = np.exp(-np.outer(tau, np.clip(lam, 0, None))).mean(1)
        ds = 2 * tau * (np.exp(-np.outer(tau, np.clip(lam, 0, None))) * lam).mean(1) / P   # = 2 tau <lambda>
        curves[phi] = ds
    print(f"[D] Newman-Watts ring, N = {N}: maximum of d_s(tau) = 2 tau <lambda>_tau")
    for phi, ds in curves.items():
        xi = np.inf if phi == 0 else 1 / phi
        print(f"    shortcut density phi = {phi:<6}: xi = 1/phi = {xi:>6}, max d_s = {ds.max():.3f} at tau = {tau[np.argmax(ds)]:.3g}")
    savedata("newman_watts_ds", ["tau"] + [f"ds_phi_{phi}" for phi in phis], [tau] + list(curves.values()),
             [f"ring of N={N} with Poisson(phi N) random shortcuts; d_s = 2 tau <lambda>_tau"])
    return tau, curves


def main():
    plt.rcParams["font.family"] = "DejaVu Sans"
    lam, nu, starts, flows = part_a()
    Lam, theta, zeta, ds_exact, arr = part_b()
    tau_c, out, win = part_c(Lam * 2 ** zeta)
    tau_d, curves = part_d()

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5), constrained_layout=True)
    a = ax[0, 0]
    x = np.linspace(0, 1, 400)
    a.plot(x, f(x), label="p' = 1 - (1 - p^2)^2")
    a.plot(x, x, "k--", lw=0.8)
    for traj in flows[1:5]:
        a.plot(np.repeat(traj[:6], 2)[1:], np.repeat(traj[1:7], 2)[:-1], lw=0.7)
    a.plot([0, PC, 1], [0, PC, 1], "o", color="k")
    a.annotate(f"p_c = {PC:.4f}\nlambda = {lam:.4f}, nu = {nu:.4f}", (PC, PC), (0.05, 0.75), arrowprops=dict(arrowstyle="->"))
    a.set_xlabel("p"); a.set_ylabel("p'"); a.set_title("Exact RG map of the DHL: fixed points 0, p_c, 1")
    a.legend(fontsize=8, loc="lower right")

    a = ax[0, 1]
    lev = arr[:, 0]
    a.semilogy(lev, arr[:, 3], "k-", label=f"exact recursion, growth Lambda = {Lam:.4f}")
    a.semilogy(lev, arr[:, 2], "o", ms=4, label="population dynamics")
    a.semilogy(lev, 2 ** (arr[:, 5] * lev), "s", ms=4, label=f"exp E[ln R | c] (growth 2^zeta, zeta = {zeta:.3f})")
    a.semilogy(lev, 2.0 ** lev, ":", color="grey", label="2^n")
    a.set_xlabel("level n (length 2^n)")
    a.set_title(f"Critical cluster: d_f = log2 Lambda = {theta:.4f}, d_s = 2 d_f/(d_f + zeta) = {ds_exact:.3f}", fontsize=10)
    a.legend(fontsize=8)

    a = ax[1, 0]
    for label, (P, ds, size, dq) in out.items():
        a.semilogx(tau_c, ds, label=label)
    a.axhline(ds_exact, color="k", ls="--", lw=0.8, label=f"Einstein relation {ds_exact:.3f}")
    a.axhline(2, color="grey", ls=":", lw=0.8)
    a.axvspan(*win, color="grey", alpha=0.1)
    a.set_ylim(0, 2.6)
    a.set_xlabel("tau"); a.set_ylabel("d_s(tau) = -2 dln<P>/dln tau")
    a.set_title("Heat kernel of DHL clusters (shaded: two RG periods)")
    a.legend(fontsize=8, loc="upper right")

    a = ax[1, 1]
    for phi, ds in curves.items():
        a.semilogx(tau_d, ds, label=f"shortcut density {phi}" + ("" if phi else " (ring, d_s = 1)"))
    a.axhline(1, color="grey", ls=":", lw=0.8)
    a.set_xlabel("tau"); a.set_ylabel("d_s(tau) = 2 tau <lambda>")
    a.set_title("Rewiring is relevant: shortcuts destroy the plateau beyond tau ~ xi^2")
    a.legend(fontsize=8)
    savefig(fig, "rg_spectral_dimension")
    print(f"figures written to {FIG_DIR}/, data to {DATA_DIR}/")
    if "--show" in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
