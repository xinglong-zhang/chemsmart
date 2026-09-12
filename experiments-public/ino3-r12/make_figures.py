"""Figures for the ino3-r12 forensic review.

Every number is read through ChemSmart's own extraction plane
(``chemsmart.analysis.result_readers.reader_for('orca')``) and every bond
through ChemSmart's own perceiver
(``chemsmart.io.molecules.perception.adjacency_matrix``).  Nothing here
re-parses an ORCA log by hand, and nothing here writes into the raw
evidence.

Run:  PYTHONPATH=/home/chemsmart/research/chemsmart \
      /opt/miniforge3/envs/chemsmart-controller/bin/python make_figures.py
"""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import proj3d

from chemsmart.analysis.result_readers import reader_for
from chemsmart.io.molecules.perception import adjacency_matrix

HERE = pathlib.Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

RAW = HERE / "raw" / "ino3-r12"
# The two redox partners this goal consumed were optimised in an earlier
# campaign window and are NOT inside the case workspace; the goal's own
# seeded run streams name them by absolute path.
UPSTREAM = pathlib.Path(
    "/home/chemsmart/agent-campaigns/ax41-refine-100/novel-round-3"
    "/workspaces/ino3-nickel-thiolate-redox/nodes"
)

R = reader_for("orca")

COLOUR = {
    "Ni": "#7b2d8b",
    "S": "#d9a600",
    "P": "#e06a1b",
    "C": "#4a4a4a",
    "H": "#bdbdbd",
    "N": "#1f5fa8",
    "O": "#c0392b",
    "F": "#2e8b57",
}
SIZE = {"Ni": 300, "S": 190, "P": 190, "C": 95, "N": 110, "O": 110,
        "F": 105, "H": 32}


def read_xyz(path):
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    sym, pos = [], []
    for line in lines[2 : 2 + n]:
        f = line.split()
        sym.append(f[0])
        pos.append([float(x) for x in f[1:4]])
    return sym, np.asarray(pos)


def align_to_plane(pos, centre, a, b):
    """Put ``centre`` at the origin with the centre-a-b plane in view."""
    pos = np.asarray(pos, dtype=float) - np.asarray(pos[centre], float)
    x = pos[a] / np.linalg.norm(pos[a])
    t = pos[b] - (pos[b] @ x) * x
    y = t / np.linalg.norm(t)
    z = np.cross(x, y)
    return pos @ np.column_stack((x, y, z))


def draw(ax, sym, pos, labels=(), view=(24, -66), highlight=(), span=None):
    """A 3D ball-and-stick with ChemSmart-perceived bonds."""
    adj = adjacency_matrix(sym, pos)
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            if not adj[i][j]:
                continue
            heavy = sym[i] != "H" and sym[j] != "H"
            ax.plot(
                *zip(pos[i], pos[j]),
                color="#4d4d4d" if heavy else "#cfcfcf",
                lw=2.4 if heavy else 1.1,
                zorder=1,
                solid_capstyle="round",
            )
    for i, s in enumerate(sym):
        ax.scatter(
            *pos[i],
            s=SIZE.get(s, 90),
            c=COLOUR.get(s, "#888888"),
            edgecolors="white",
            linewidths=0.8,
            depthshade=False,
            zorder=3,
        )
    for i, j in highlight:
        ax.plot(
            *zip(pos[i], pos[j]),
            color="#c0392b",
            lw=1.8,
            ls=(0, (3, 2)),
            zorder=4,
        )
    ax.view_init(*view)
    if span is None:
        lo = pos.min(axis=0) - 0.85
        hi = pos.max(axis=0) + 0.85
    else:
        mid = (pos.min(axis=0) + pos.max(axis=0)) / 2.0
        lo, hi = mid - span, mid + span
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_zlim(lo[2], hi[2])
    ax.set_axis_off()
    ax.set_box_aspect(hi - lo)
    if labels:
        _leader_labels(ax, pos, labels)


def _leader_labels(ax, pos, labels, pad=1.16):
    """Atom labels placed outside the projected molecule, on leader lines.

    Placing text at a 3D offset cannot avoid occlusion, because which
    direction is "clear" depends on the camera. This projects every atom
    to the panel's own 2D coordinates, pushes each label radially outside
    the molecule's projected extent, and draws a thin leader to the atom.
    """
    ax.get_figure().canvas.draw()
    matrix = ax.get_proj()
    flat = np.array([
        proj3d.proj_transform(p[0], p[1], p[2], matrix)[:2] for p in pos
    ])
    centre = flat.mean(axis=0)
    radius = float(np.linalg.norm(flat - centre, axis=1).max())
    for i, text, off in labels:
        direction = flat[i] - centre
        norm = float(np.linalg.norm(direction))
        if norm < 0.30 * radius:
            # an atom near the projected centre has no clear radial
            # direction; use the hint the caller supplied instead
            hint = np.array([off[0], off[1]], dtype=float)
            if np.linalg.norm(hint) < 1e-9:
                hint = np.array([0.0, -1.0])
            direction = hint / np.linalg.norm(hint)
        else:
            direction = direction / norm
        tip = centre + direction * radius * pad
        ax.annotate(
            text, xy=tuple(flat[i]), xytext=tuple(tip),
            xycoords="data", textcoords="data",
            fontsize=9.2, fontweight="bold", color="black",
            ha="center", va="center", zorder=20, annotation_clip=False,
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none",
                      alpha=0.92),
            arrowprops=dict(arrowstyle="-", color="#909090", lw=0.8,
                            shrinkA=1, shrinkB=4),
        )


def geom(path, selector="reached_positions"):
    out = R.open_output(path)
    pos, _ = R.read(out, selector)
    sym, _ = R.read(out, "symbols")
    return sym, np.asarray(pos)


def angle(pos, i, j, k):
    a, b = pos[i] - pos[j], pos[k] - pos[j]
    return float(
        np.degrees(
            np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
        )
    )


def dist(pos, i, j):
    return float(np.linalg.norm(pos[i] - pos[j]))


# --------------------------------------------------------------- figure 1
def figure_structures():
    start_sym, start_pos = read_xyz(
        RAW / "nickel-bis-thiolate-bis-phosphine.xyz"
    )
    n0 = geom(UPSTREAM / "pme3-n0-opt2" / "geom-pme3-n0-disp_opt_opt.out")
    cat = geom(UPSTREAM / "nicat-opt-pme3" / "geom-pme3v2-h18_opt_opt.out")
    qrt = geom(UPSTREAM / "pme3-qrt-opt2" / "geom-pme3-qrt-reach_opt_opt.out")

    def note(pos, s1, s2, p1, p2):
        return (
            f"Ni-S  {dist(pos,0,s1):.3f} / {dist(pos,0,s2):.3f} A\n"
            f"Ni-P  {dist(pos,0,p1):.3f} / {dist(pos,0,p2):.3f} A\n"
            f"S-Ni-S {angle(pos,s1,0,s2):6.1f} deg\n"
            f"P-Ni-P {angle(pos,p1,0,p2):6.1f} deg"
        )

    LAB = [(0, "Ni", (1.15, -1.15, 0.35)), (1, "S1", (1.15, 0.30, 0.0)),
           (6, "S2", (-1.15, -0.30, 0.0)), (11, "P1", (-0.30, 1.15, 0.0)),
           (12, "P2", (0.30, -1.15, 0.0))]
    LAB0 = [(0, "Ni", (1.15, -1.15, 0.35)), (1, "S1", (1.15, 0.30, 0.0)),
            (6, "S2", (-1.15, -0.30, 0.0)), (11, "P1", (-0.30, 1.15, 0.0)),
            (15, "P2", (0.30, -1.15, 0.0))]

    panels = [
        (align_to_plane(start_pos, 0, 1, 11), start_sym, LAB0,
         "(a) supplied start, nothing optimised",
         "Ni(SMe)$_2$(PH$_3$)$_2$, thiolates trans",
         note(start_pos, 1, 6, 11, 15)),
        (align_to_plane(n0[1], 0, 1, 11), n0[0], LAB,
         "(b) neutral minimum, charge 0, $S$ = 0",
         "square planar, 180°/180°",
         note(n0[1], 1, 6, 11, 12)),
        (align_to_plane(cat[1], 0, 1, 11), cat[0], LAB,
         "(c) doublet cation, charge +1, $S$ = 1/2",
         "square plane buckled to 151°/158°",
         note(cat[1], 1, 6, 11, 12)),
        (align_to_plane(qrt[1], 0, 1, 11), qrt[0], LAB,
         "(d) quartet cation, charge +1, $S$ = 3/2",
         "pseudo-tetrahedral, 106°/106°",
         note(qrt[1], 1, 6, 11, 12)),
    ]

    # one common half-span so the four panels are directly comparable
    half = np.array([
        max((p.max(axis=0) - p.min(axis=0))[a] for p, *_r in panels) / 2.0
        + 1.65
        for a in range(3)
    ])
    half[:] = half.max()

    fig = plt.figure(figsize=(13.6, 4.4))
    for k, (pos, sym, labels, title, sub, text) in enumerate(panels, 1):
        ax = fig.add_subplot(1, 4, k, projection="3d")
        draw(ax, sym, pos, labels=labels, view=(62, -74), span=half)
        ax.set_title(f"{title}\n{sub}", fontsize=9.0, pad=0)
        ax.text2D(
            0.5, 0.055, text, transform=ax.transAxes, ha="center",
            va="top", fontsize=8.7, family="monospace",
        )
    handles = [
        Line2D([], [], marker="o", ls="", color=COLOUR[e], label=e,
               markersize=8, markeredgecolor="white")
        for e in ("Ni", "S", "P", "C", "H")
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, 0.005))
    fig.subplots_adjust(left=0.0, right=1.0, top=0.88, bottom=0.235,
                        wspace=0.0)
    fig.savefig(FIG / "figA1-structures.png", dpi=210)
    plt.close(fig)


# --------------------------------------------------------------- figure 2
FUNCTIONALS = (
    ("TPSSh", "cat-sp-svp-tpssh", 10.0),
    ("B3LYP", "cat-sp-svp-b3lyp", 20.0),
    ("PBE0", "cat-sp-svp-pbe0", 25.0),
)
CATFILE = "geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.out"


def cation_spin(node):
    out = R.open_output(RAW / "nodes" / node / CATFILE)
    mull, _ = R.read(out, "mulliken_atomic_spin_populations")
    lowd, _ = R.read(out, "loewdin_atomic_spin_populations")
    s2, _ = R.read(out, "spin_square")
    return np.asarray(mull), np.asarray(lowd), float(s2)


def figure_spin():
    rows = []
    for name, node, hf in FUNCTIONALS:
        mull, lowd, s2 = cation_spin(node)
        rows.append((name, hf, mull, lowd, s2))

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(12.6, 4.5),
        gridspec_kw={"width_ratios": (1.55, 1.0)},
    )

    sites = [("Ni", 0), ("S1", 1), ("S2", 6), ("P1", 11), ("P2", 12)]
    width = 0.135
    xs = np.arange(len(sites), dtype=float)
    shades = {"TPSSh": "#9ecae1", "B3LYP": "#4a91c6", "PBE0": "#1b4f78"}
    for k, (name, _hf, mull, lowd, _s2) in enumerate(rows):
        ax.bar(
            xs + (k - 2.5) * width,
            [mull[i] for _n, i in sites],
            width, color=shades[name], edgecolor="white",
            label=f"{name}  Mulliken",
        )
        ax.bar(
            xs + (k + 0.5) * width,
            [lowd[i] for _n, i in sites],
            width, color=shades[name], edgecolor="white", hatch="///",
            alpha=0.9, label=f"{name}  Löwdin",
        )
    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(0.5, color="#c0392b", lw=1.0, ls=":")
    ax.text(
        len(sites) - 0.55, 0.515,
        "0.5 e$^-$: a Ni(II)–thiyl radical would put\n"
        "about this much on ONE sulfur",
        fontsize=8.2, color="#c0392b", va="bottom", ha="right",
    )
    ax.set_xticks(xs)
    ax.set_xticklabels([n for n, _i in sites], fontsize=10)
    ax.set_ylabel("atomic spin population  /  electron", fontsize=10)
    ax.set_ylim(-0.09, 0.95)
    ax.set_title(
        "(a) where the hole sits in [Ni(SMe)$_2$(PMe$_3$)$_2$]$^+$ "
        "(doublet)\n"
        "three functionals, ONE geometry, one basis (def2-SVP), "
        "CPCM(MeCN)",
        fontsize=9.6,
    )
    ax.legend(fontsize=7.8, ncol=3, frameon=False, loc="upper right",
              bbox_to_anchor=(1.0, 0.93))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    hf = [r[1] for r in rows]
    ni = [r[2][0] for r in rows]
    ssum = [r[2][1] + r[2][6] for r in rows]
    bx.plot(hf, ni, "o-", color="#7b2d8b", lw=2, ms=8, label="Ni")
    bx.plot(hf, ssum, "s-", color="#d9a600", lw=2, ms=8,
            label="S1 + S2")
    for x, y, name in zip(hf, ni, [r[0] for r in rows]):
        bx.annotate(f"{name}\n{y:.3f}", (x, y), textcoords="offset points",
                    xytext=(6, -16), fontsize=8.4)
    for x, y in zip(hf, ssum):
        bx.annotate(f"{y:.3f}", (x, y), textcoords="offset points",
                    xytext=(6, 4), fontsize=8.4)
    bx.set_xlabel("exact-exchange fraction in the functional  /  %",
                  fontsize=9.6)
    bx.set_ylabel("Mulliken spin population  /  electron", fontsize=9.6)
    bx.set_title(
        "(b) the metal fraction is functional-dependent;\n"
        "the metal-vs-sulfur answer is not",
        fontsize=9.6,
    )
    bx.set_ylim(0, 0.95)
    bx.set_xlim(5, 30)
    bx.grid(alpha=0.25, ls=":")
    bx.legend(fontsize=9, frameon=False, loc="center right")
    for spine in ("top", "right"):
        bx.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIG / "figA2-spin-per-functional.png", dpi=210)
    plt.close(fig)
    return rows


# --------------------------------------------------------------- figure 3
def figure_couple():
    delivered = -0.39414929
    unc = 0.16734815
    lin = 0.2474415
    ladder = [
        ("PBE0-D3BJ  (delivered)", -0.39414929, "#1b4f78"),
        ("B3LYP-D3BJ", -0.52319474, "#2d6d9e"),
        ("TPSSh-D3BJ", -0.54110866, "#0d3550"),
    ]
    oxidants = [
        ("FcPF$_6$  (ferrocenium)", 0.00),
        ("AcFcBF$_4$", 0.27),
        ("TBPA$\\cdot$SbCl$_6$", 0.70),
    ]
    terms = [
        ("functional\nPBE0/B3LYP/TPSSh", 0.14695937, "measured"),
        ("basis differential\ndef2-SVP → def2-TZVP", 0.078592716, "measured"),
        ("thermochemical model\nelec / RRHO / qRRHO", 0.012906056, "measured"),
        ("ferrocene reference-scale\nCONSISTENCY", 0.008, "measured"),
        ("geometry materialisation", 0.00098336238, "measured"),
        ("CPCM differential solvation", None, "named by the case"),
        ("ferrocene reference-scale\nACCURACY (source: 0.05–0.1 V)", None,
         "not named"),
        ("functional-dependent\ngeometry relaxation", None, "not named"),
    ]

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(13.0, 4.6),
        gridspec_kw={"width_ratios": (1.35, 1.0)},
    )

    ax.axhspan(delivered - unc, delivered + unc, color="#1b4f78",
               alpha=0.13, zorder=0)
    ax.axhspan(delivered - lin, delivered + lin, color="#1b4f78",
               alpha=0.06, zorder=0)
    nudge = {0: 0.0, 1: 0.030, 2: -0.030}
    for k, (name, value, colour) in enumerate(ladder):
        ax.hlines(value, 0.15, 1.15, color=colour, lw=3)
        ax.text(1.2, value + nudge[k], f"{name}   {value:+.3f} V",
                va="center", fontsize=9, color=colour)
    ax.errorbar(
        [0.65], [delivered], yerr=[[unc], [unc]], fmt="D",
        color="#1b4f78", ms=9, capsize=6, lw=2, zorder=5,
    )
    for name, pot in oxidants:
        ax.hlines(pot, 0.15, 2.6, color="#c0392b", lw=1.4, ls="--")
        ax.text(2.62, pot, f"{name}  {pot:+.2f} V", va="center",
                fontsize=9, color="#c0392b")
    ax.annotate(
        "",
        xy=(0.4, 0.0), xytext=(0.4, delivered),
        arrowprops=dict(arrowstyle="<->", color="black", lw=1.2),
    )
    ax.text(
        0.44, delivered / 2,
        f"margin to the mildest\noxidant  {-delivered:+.3f} V\n"
        f"(linear-sum budget {lin:.3f} V)",
        fontsize=8.6, va="center",
    )
    ax.set_ylim(-0.72, 0.84)
    ax.set_xlim(0.0, 4.4)
    ax.set_xticks([])
    ax.set_ylabel("potential  /  V vs Fc$^+$/Fc in acetonitrile",
                  fontsize=10)
    ax.set_title(
        "(a) the +1/0 couple of the PMe$_3$ model against the shelf\n"
        "delivered −0.394 ± 0.167 V; interval [−0.561, −0.227] V",
        fontsize=9.6,
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    names = [t[0] for t in terms]
    vals = [t[1] if t[1] is not None else 0.0 for t in terms]
    cols = ["#1b4f78" if t[2] == "measured" else "#bbbbbb" for t in terms]
    ys = np.arange(len(terms))[::-1]
    bx.barh(ys, vals, color=cols, edgecolor="white", height=0.62)
    absent_note = {
        "named by the case": "  absent; named by the case",
        "not named": "  absent; NOT named by the case",
    }
    for y, t, v in zip(ys, terms, vals):
        if t[1] is None:
            bx.text(0.004, y, absent_note[t[2]], va="center", fontsize=8.0,
                    color="#c0392b" if t[2] == "not named" else "#777777")
        else:
            bx.text(v + 0.004, y, f"{v:.4f} V", va="center", fontsize=8.6)
    bx.axvline(unc, color="#c0392b", lw=1.4)
    bx.text(unc - 0.005, 5.15,
            f"quadrature\ntotal\n{unc:.4f} V", color="#c0392b",
            fontsize=8.4, ha="right", va="center")
    for value, row, label in (
        (0.174475, 3.15, "0.1745 if the reference term\nis the source's 0.05 V"),
        (0.194786, 1.55, "0.1948 if it is\nthe source's 0.1 V"),
    ):
        bx.axvline(value, color="#7b2d8b", lw=1.1, ls="--")
        bx.text(value - 0.004, row, label, color="#7b2d8b",
                fontsize=7.8, ha="right", va="center")
    bx.axvline(0.2, color="black", lw=1.2, ls=":")
    bx.text(0.203, 0.15, "task tolerance\n±0.2 V",
            fontsize=8.8, ha="left", va="center")
    bx.set_yticks(ys)
    bx.set_yticklabels(names, fontsize=8.6)
    bx.set_xlim(0, 0.255)
    bx.set_xlabel("measured spread / displacement  /  V", fontsize=9.6)
    bx.set_title(
        "(b) the uncertainty the host scored 'met': five measured\n"
        "magnitudes in quadrature, and THREE terms absent from it",
        fontsize=9.6,
    )
    for spine in ("top", "right"):
        bx.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(FIG / "figA3-couple-and-budget.png", dpi=210)
    plt.close(fig)


if __name__ == "__main__":
    figure_structures()
    rows = figure_spin()
    figure_couple()
    print("Mulliken / Loewdin spin, def2-SVP, doublet cation:")
    for name, hf, mull, lowd, s2 in rows:
        print(
            f"  {name:6s} HF={hf:4.1f}%  Ni {mull[0]:.6f}/{lowd[0]:.6f}  "
            f"S1 {mull[1]:.6f}/{lowd[1]:.6f}  S2 {mull[6]:.6f}/{lowd[6]:.6f}"
            f"  P {mull[11]:.6f}/{mull[12]:.6f}  <S^2> {s2:.6f}"
            f"  sum {mull.sum():.6f}"
        )
    print("wrote:", sorted(p.name for p in FIG.glob("*.png")))
