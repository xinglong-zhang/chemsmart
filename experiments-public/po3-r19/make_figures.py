"""Figures for the po3-r19 forensic review.

Energies, coordinates, frequencies and scan surfaces are read through
ChemSmart's own extraction plane
(``chemsmart.analysis.result_readers.reader_for('orca')``); bonds and
atom roles come from ChemSmart's own perceiver
(``chemsmart.io.molecules.perception.adjacency_matrix``).  Per-scan-point
structures are ORCA's own ``.00N.xyz`` sidecars, read as plain XYZ.
Nothing here writes into the raw evidence.

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
RAW = HERE / "raw" / "po3-r19"
NODES = RAW / "nodes"
RUNS = RAW / ".chemsmart-agent" / "runs"

R = reader_for("orca")
HARTREE_KCAL = 627.5094740631

COLOUR = {"C": "#4a4a4a", "H": "#c9c9c9", "N": "#1f5fa8", "O": "#c0392b",
          "F": "#2e8b57"}
SIZE = {"C": 95, "N": 120, "O": 115, "F": 105, "H": 30}


def read_xyz(path):
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    sym, pos = [], []
    for line in lines[2 : 2 + n]:
        f = line.split()
        sym.append(f[0])
        pos.append([float(x) for x in f[1:4]])
    return sym, np.asarray(pos)


def roles(sym, pos):
    """Chemical roles from the perceived graph, never from atom order."""
    adj = np.asarray(adjacency_matrix(sym, pos))
    nb = lambda i: [j for j in range(len(sym)) if adj[i][j]]
    cf3 = [
        i for i, s in enumerate(sym)
        if s == "C" and sum(1 for j in nb(i) if sym[j] == "F") == 3
    ][0]
    carbonyl = [
        i for i, s in enumerate(sym)
        if s == "C" and sum(1 for j in nb(i) if sym[j] == "O") == 2
    ][0]
    c_cf3 = [j for j in nb(cf3) if sym[j] == "C"][0]
    cand = [
        j for j in nb(carbonyl)
        if sym[j] == "C"
        and j != c_cf3
        and sum(1 for k in nb(j) if sym[k] == "H") == 0
    ]
    c_est = cand[0]
    ch2 = [
        i for i, s in enumerate(sym)
        if s == "C"
        and sum(1 for j in nb(i) if sym[j] == "H") == 2
        and any(sym[j] == "N" for j in nb(i))
    ][0]
    n1 = [j for j in nb(ch2) if sym[j] == "N"][0]
    n2 = [j for j in nb(n1) if sym[j] == "N"][0]
    rest = [j for j in nb(n2) if sym[j] == "N" and j != n1]
    n3 = rest[0] if rest else [
        i for i, s in enumerate(sym) if s == "N" and i not in (n1, n2)
    ][0]
    return dict(n1=n1, n2=n2, n3=n3, c_est=c_est, c_cf3=c_cf3, cf3=cf3,
                carbonyl=carbonyl, ch2=ch2)


def dist(pos, i, j):
    return float(np.linalg.norm(pos[i] - pos[j]))


def align(pos, a, b, c):
    """Put the a-b-c plane in the page; a at the origin."""
    pos = np.asarray(pos, float) - pos[a]
    x = (pos[b] - pos[a])
    x = x / np.linalg.norm(x)
    t = (pos[c] - pos[a])
    t = t - (t @ x) * x
    y = t / np.linalg.norm(t)
    return pos @ np.column_stack((x, y, np.cross(x, y)))


def _leader_labels(ax, pos, labels, pad=1.16):
    """Atom labels placed outside the projected molecule, on leader lines.

    A 3D text offset cannot avoid occlusion, because which direction is
    clear depends on the camera. This projects every atom into the panel's
    own 2D coordinates, pushes each label radially outside the molecule's
    projected extent, and draws a thin leader back to the atom. An atom
    near the projected centre has no radial direction, so the caller's
    offset is used as a screen-space hint instead.
    """
    ax.get_figure().canvas.draw()
    matrix = ax.get_proj()
    flat = np.array([
        proj3d.proj_transform(q[0], q[1], q[2], matrix)[:2] for q in pos
    ])
    centre = flat.mean(axis=0)
    radius = float(np.linalg.norm(flat - centre, axis=1).max())
    for entry in labels:
        i, text, off = entry[:3]
        forced = bool(entry[3]) if len(entry) > 3 else False
        direction = flat[i] - centre
        norm = float(np.linalg.norm(direction))
        if forced or norm < 0.30 * radius:
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
            fontsize=8.8, fontweight="bold", color="black",
            ha="center", va="center", zorder=20, annotation_clip=False,
            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none",
                      alpha=0.92),
            arrowprops=dict(arrowstyle="-", color="#909090", lw=0.8,
                            shrinkA=1, shrinkB=4),
        )


def draw(ax, sym, pos, labels=(), view=(24, -70), dashed=(), pad=0.85,
         span=None):
    adj = adjacency_matrix(sym, pos)
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            if not adj[i][j]:
                continue
            heavy = sym[i] != "H" and sym[j] != "H"
            ax.plot(*zip(pos[i], pos[j]),
                    color="#4d4d4d" if heavy else "#d0d0d0",
                    lw=2.2 if heavy else 1.0, zorder=1,
                    solid_capstyle="round")
    for i, s in enumerate(sym):
        ax.scatter(*pos[i], s=SIZE.get(s, 90), c=COLOUR.get(s, "#888888"),
                   edgecolors="white", linewidths=0.7, depthshade=False,
                   zorder=3)
    for entry in dashed:
        i, j, colour, text = entry[:4]
        frac = entry[4] if len(entry) > 4 else 0.5
        ax.plot(*zip(pos[i], pos[j]), color=colour, lw=2.0,
                ls=(0, (2.5, 2)), zorder=4)
        if text:
            spot = pos[i] + frac * (pos[j] - pos[i])
            ax.text(*spot, text, fontsize=8.4, color=colour, zorder=6,
                    ha="center", va="center", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.1", fc="white",
                              ec="none", alpha=0.9))
    ax.view_init(*view)
    if span is None:
        lo, hi = pos.min(axis=0) - pad, pos.max(axis=0) + pad
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


SCANS = (
    ("scan-esterc4-path", "path A", "driven N1$\\cdots$C(CF$_3$)",
     "c_cf3", "#1b4f78"),
    ("scan-esterc5-path", "path B", "driven N1$\\cdots$C(ester)",
     "c_est", "#b8570e"),
)
SCANOUT = "ts-guess-esterc4-complex_scan_scan.out"


def scan_table(node):
    out = R.open_output(NODES / node / SCANOUT)
    cv, _ = R.read(out, "scan_coordinate_values")
    en, _ = R.read(out, "scan_energies")
    sym, _ = R.read(out, "symbols")
    rows = []
    for k in range(1, 8):
        s, p = read_xyz(NODES / node / SCANOUT.replace(".out", ".%03d.xyz" % k))
        r = roles(s, p)
        rows.append(dict(
            point=k, coord=cv[k - 1],
            kcal=(en[k - 1] - min(en)) * HARTREE_KCAL,
            n1_est=dist(p, r["n1"], r["c_est"]),
            n1_cf3=dist(p, r["n1"], r["c_cf3"]),
            n3_est=dist(p, r["n3"], r["c_est"]),
            n3_cf3=dist(p, r["n3"], r["c_cf3"]),
            cc=dist(p, r["c_est"], r["c_cf3"]),
            sym=s, pos=p, roles=r,
        ))
    return rows


# --------------------------------------------------------------- figure 1
def figure_scans(tables):
    fig, axes = plt.subplots(
        2, 2, figsize=(13.2, 7.4),
        gridspec_kw={"height_ratios": (1.2, 1.0), "hspace": 0.30,
                     "wspace": 0.30},
    )
    emax = max(max(r["kcal"] for r in tables[n]) for n, *_r in SCANS)
    for col, (node, tag, driven, _partner, colour) in enumerate(SCANS):
        rows = tables[node]
        x = [r["coord"] for r in rows]
        e = [r["kcal"] for r in rows]
        lo, hi = min(x) - 0.13, max(x) + 0.13

        ax = axes[0][col]
        ax.plot(x, e, "-", color=colour, lw=2, zorder=2)
        ax.plot(x, e, "o", color="white", mec=colour, mew=2, ms=9, zorder=3)
        kmax, kmin = int(np.argmax(e)), int(np.argmin(e))
        ax.plot(x[kmax], e[kmax], "o", color="#c0392b", ms=13, zorder=4)
        ax.plot(x[kmin], e[kmin], "o", color="#2e8b57", ms=11, zorder=4)
        for r in rows:
            ax.annotate(str(r["point"]), (r["coord"], r["kcal"]),
                        textcoords="offset points", xytext=(0, 12),
                        ha="center", fontsize=8.6, color="#333333")
        ax.annotate(
            "host anomaly: scan.extremum_at_grid_boundary\n"
            f"maximum at the grid END, {x[kmax]:.2f} \u00c5, point 7 of 7\n"
            f"last step +{e[kmax] - e[kmax - 1]:.3f} kcal/mol, "
            f"span {max(e) - min(e):.3f} kcal/mol\n"
            "no point on this surface was carried forward",
            xy=(x[kmax], e[kmax]),
            xytext=(x[3], emax * 0.80),
            fontsize=8.4, color="#c0392b", ha="center", va="top",
            arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.1),
        )
        ax.annotate(f"interior minimum, point {kmin + 1}\n"
                    "a loose complex, not a ridge",
                    xy=(x[kmin], e[kmin]),
                    xytext=(x[kmin], emax * 0.30),
                    fontsize=8.4, color="#2e8b57", ha="center",
                    arrowprops=dict(arrowstyle="->", color="#2e8b57",
                                    lw=1.1))
        ax.set_title(
            f"{tag}  ·  node {node}\nrelaxed scan, B3LYP-D3BJ/def2-SVP, "
            f"7 points, {driven}",
            fontsize=9.6,
        )
        ax.set_ylabel("relative energy  /  kcal mol$^{-1}$", fontsize=9.6)
        ax.set_xlabel("driven coordinate  /  Å", fontsize=9.4)
        ax.grid(alpha=0.25, ls=":")
        ax.set_xlim(hi, lo)
        ax.set_ylim(-1.6, emax * 1.16)   # one scale for both panels
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        bx = axes[1][col]
        # the bond that must ALSO form: terminal N to the carbon that
        # becomes ring C4 in this path's own product
        partner_key = "n3_est" if col == 0 else "n3_cf3"
        nearest_key = "n3_cf3" if col == 0 else "n3_est"
        other = [r[partner_key] for r in rows]
        bx.plot(x, other, "s-", color="#7b2d8b", lw=2, ms=8)
        bx.axhline(1.47, color="#888888", lw=1.1, ls="--")
        bx.text(hi - 0.05, 1.36,
                "1.47 \u00c5 \u2014 what this N\u2013C bond is in the product",
                fontsize=8.3, color="#666666", ha="left", va="top")
        bx.annotate(
            f"{other[0]:.2f} \u00c5 at point 1 \u2192 "
            f"{other[-1]:.2f} \u00c5 at point 7"
            + (f"\nit RETREATS by {other[-1] - other[0]:.2f} \u00c5 "
               "as the driven bond closes"
               if other[-1] > other[0] + 0.5
               else f"\nit moves only {other[-1] - other[0]:+.2f} \u00c5: "
                    "it does not close")
            + f"\n(nearest alkyne C at point 7: "
              f"{rows[-1][nearest_key]:.2f} \u00c5)",
            xy=(x[-1], other[-1]), xytext=(x[3], 6.55),
            fontsize=8.4, color="#7b2d8b", ha="center", va="top",
            arrowprops=dict(arrowstyle="->", color="#7b2d8b", lw=1.1),
        )
        cx = bx.twinx()
        cx.plot(x, [r["cc"] for r in rows], "^-", color="#2e8b57", lw=2,
                ms=8)
        saddle_cc = 1.2421 if col == 0 else 1.2413
        cx.axhline(saddle_cc, color="#2e8b57", lw=1.0, ls=":")
        cx.text(hi - 0.05, saddle_cc - 0.0055,
                f"C$\\equiv$C at the converged saddle, "
                f"{saddle_cc:.4f} \u00c5",
                fontsize=8.2, color="#2e8b57", ha="left", va="top")
        cx.set_ylim(1.200, 1.262)
        cx.set_ylabel("alkyne C$\\equiv$C  /  Å", fontsize=9.4,
                      color="#2e8b57")
        cx.tick_params(axis="y", colors="#2e8b57")
        bx.set_ylim(1.1, 6.7)
        bx.set_xlim(hi, lo)
        bx.set_xlabel("driven coordinate  /  Å", fontsize=9.6)
        bx.set_ylabel("terminal N$\\cdots$C4-to-be  /  Å",
                      fontsize=9.4, color="#7b2d8b")
        bx.tick_params(axis="y", colors="#7b2d8b")
        bx.grid(alpha=0.25, ls=":")
        bx.set_title(
            "the OTHER forming bond (purple, left) never closes;\n"
            "the triple bond (green, right) never rehybridises",
            fontsize=9.2,
        )
        for spine in ("top",):
            bx.spines[spine].set_visible(False)
    fig.subplots_adjust(left=0.075, right=0.935, top=0.90, bottom=0.085)
    fig.savefig(FIG / "figB1-scans.png", dpi=210)
    plt.close(fig)


# --------------------------------------------------------------- figure 2
def figure_scan_structures(tables):
    show = (1, 3, 5, 7)
    frames = {}
    for node, tag, driven, partner, _c in SCANS:
        rows = []
        for pt in show:
            r = tables[node][pt - 1]
            rl = r["roles"]
            rows.append((r, rl,
                         align(r["pos"], rl["c_est"], rl["c_cf3"],
                               rl["carbonyl"])))
        frames[node] = rows
    # one scale for the whole figure: the same molecule in every panel
    everything = np.vstack([a for rows in frames.values() for _r, _l, a in rows])
    half = float((everything.max(axis=0) - everything.min(axis=0)).max()) / 2.0
    span = np.array([half + 0.55] * 3)

    fig = plt.figure(figsize=(14.2, 6.6))
    for row, (node, tag, driven, partner, _c) in enumerate(SCANS):
        for col, (r, rl, apos) in enumerate(frames[node]):
            ax = fig.add_subplot(2, 4, row * 4 + col + 1, projection="3d")
            driven_pair = (rl["n1"], rl[partner])
            c4_to_be = rl["c_est"] if row == 0 else rl["c_cf3"]
            other = (rl["n3"], c4_to_be)
            draw(
                ax, r["sym"], apos,
                labels=[
                    (rl["n1"], "N1", (0.0, 1.0, 0.0)),
                    (rl["n3"], "N3", (0.55, 1.0, 0.0)),
                    (rl["c_est"], "C(ester)", (-1.0, -0.55, 0.0), True),
                    (rl["c_cf3"], "C(CF$_3$)", (1.0, -0.55, 0.0), True),
                ],
                dashed=[
                    (*driven_pair, "#c0392b", "", 0.42),
                    (*other, "#7b2d8b", "", 0.68),
                ],
                view=(16, -80), span=span,
            )
            ax.set_title(
                f"{tag}, point {r['point']} of 7   "
                f"$E-E_{{min}}$ = {r['kcal']:+.2f} kcal/mol",
                fontsize=9.2, pad=-4,
            )
            ax.text2D(
                0.5, 0.095,
                f"driven N1\u00b7\u00b7\u00b7C  "
                f"{dist(apos, *driven_pair):.2f} \u00c5\n"
                f"N3\u00b7\u00b7\u00b7C4-to-be  "
                f"{dist(apos, *other):.2f} \u00c5",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=8.6, family="monospace",
            )
    handles = [
        Line2D([], [], marker="o", ls="", color=COLOUR[e], label=e,
               markersize=8, markeredgecolor="white")
        for e in ("C", "N", "O", "F", "H")
    ] + [
        Line2D([], [], color="#c0392b", lw=2, ls=(0, (2.5, 2)),
               label="the DRIVEN contact  N1$\\cdots$C  /  \u00c5"),
        Line2D([], [], color="#7b2d8b", lw=2, ls=(0, (2.5, 2)),
               label="the other forming bond  N3$\\cdots$C4-to-be  /  \u00c5"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False,
               fontsize=9.4, bbox_to_anchor=(0.5, 0.004))
    fig.subplots_adjust(left=0.0, right=1.0, top=0.96, bottom=0.06,
                        wspace=0.0, hspace=0.10)
    fig.savefig(FIG / "figB2-scan-structures.png", dpi=200)
    plt.close(fig)


# --------------------------------------------------------------- figure 3
PRESADDLE = (
    RUNS / "live-20260912T023123806980Z-07481088a0-7cb5f44c" / "artifacts"
)
SADDLES = (
    ("ts-esterc4-restart", "geom-c4-saddle-reached_optts_optts.out",
     "(b) converged ester-at-C4 saddle\ncycle 3, restarted from the "
     "reached geometry"),
    ("ts-esterc5", "presaddle-esterc5_optts_optts.out",
     "(c) converged ester-at-C5 saddle\ncycle 2, converged first time"),
)


def ring_panel(ax, sym, pos, rl, forming, cross, tilt=80.0, span=None):
    """The forming five-ring face-on: C4, C5 and N1 define the page."""
    apos = align(pos, rl["c_est"], rl["c_cf3"], rl["n1"])
    draw(
        ax, sym, apos,
        labels=[
            (rl["n1"], "N1", (-0.7, 1.0, 0.0), True),
            (rl["n2"], "N2", (0.0, 1.0, 0.0), True),
            (rl["n3"], "N3", (0.7, 1.0, 0.0), True),
            (rl["c_est"], "C(ester)", (-1.0, -0.45, 0.0), True),
            (rl["c_cf3"], "C(CF$_3$)", (1.0, -0.45, 0.0), True),
        ],
        dashed=[(i, j, "#c0392b", "", 0.5) for i, j in forming]
               + [(i, j, "#9aa0a6", "", 0.5) for i, j in cross],
        view=(tilt, -90), span=span,
    )
    return apos


def figure_saddles():
    presaddle_sym, presaddle_pos = read_xyz(PRESADDLE / "presaddle-esterc4.xyz")
    presaddle_rl = roles(presaddle_sym, presaddle_pos)
    panels = [(presaddle_sym, presaddle_pos, presaddle_rl)]
    for node, fn, _title in SADDLES:
        out = R.open_output(NODES / node / fn)
        sym, _ = R.read(out, "symbols")
        pos = np.asarray(R.read(out, "reached_positions")[0])
        panels.append((sym, pos, roles(sym, pos)))
    aligned = [align(pos, rl["c_est"], rl["c_cf3"], rl["n1"])
               for _s, pos, rl in panels]
    half = max(float((a.max(axis=0) - a.min(axis=0)).max()) for a in aligned)
    span = np.array([half / 2.0 + 0.55] * 3)

    fig = plt.figure(figsize=(13.6, 5.4))

    sym, pos, rl = panels[0]
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    apos = ring_panel(
        ax, sym, pos, rl,
        forming=[(rl["n3"], rl["c_est"]), (rl["n1"], rl["c_cf3"])],
        cross=[(rl["n3"], rl["c_cf3"]), (rl["n1"], rl["c_est"])],
        span=span,
    )
    ax.set_title(
        "(a) the ester-at-C4 PRE-SADDLE\ncut from the user's own product "
        "and rejoined",
        fontsize=9.2, pad=-4,
    )
    ax.text2D(
        0.5, 0.235,
        f"forming bonds  N3\u00b7\u00b7\u00b7C(ester) "
        f"{dist(apos, rl['n3'], rl['c_est']):.3f}\n"
        f"               N1\u00b7\u00b7\u00b7C(CF3)  "
        f"{dist(apos, rl['n1'], rl['c_cf3']):.3f}\n"
        f"C(ester)-C(CF3) "
        f"{dist(apos, rl['c_est'], rl['c_cf3']):.3f} - rehybridised\n"
        f"no contact closer than the contacts, in \u00c5",
        transform=ax.transAxes, ha="center", va="top", fontsize=8.4,
        family="monospace",
    )

    for k, (node, fn, title) in enumerate(SADDLES, start=2):
        sym, pos, rl = panels[k - 1]
        out = R.open_output(NODES / node / fn)
        freq, _ = R.read(out, "vibrational_frequencies")
        imag = [f for f in freq if f < 0]
        if k == 2:
            forming = [(rl["n3"], rl["c_est"]), (rl["n1"], rl["c_cf3"])]
            cross = [(rl["n3"], rl["c_cf3"]), (rl["n1"], rl["c_est"])]
            names = ("N3\u00b7\u00b7\u00b7C(ester)", "N1\u00b7\u00b7\u00b7C(CF3) ")
        else:
            forming = [(rl["n3"], rl["c_cf3"]), (rl["n1"], rl["c_est"])]
            cross = [(rl["n3"], rl["c_est"]), (rl["n1"], rl["c_cf3"])]
            names = ("N3\u00b7\u00b7\u00b7C(CF3) ", "N1\u00b7\u00b7\u00b7C(ester)")
        ax = fig.add_subplot(1, 3, k, projection="3d")
        apos = ring_panel(ax, sym, pos, rl, forming, cross, span=span)
        ax.set_title(title, fontsize=9.2, pad=-4)
        margin = (min(dist(apos, *q) for q in cross)
                  - max(dist(apos, *q) for q in forming))
        ax.text2D(
            0.5, 0.235,
            f"forming bonds  {names[0]} {dist(apos, *forming[0]):.3f}\n"
            f"               {names[1]} {dist(apos, *forming[1]):.3f}\n"
            f"cross contacts {dist(apos, *cross[0]):.3f} / "
            f"{dist(apos, *cross[1]):.3f}   margin {margin:.3f}\n"
            f"C-C {dist(apos, rl['c_est'], rl['c_cf3']):.3f} \u00c5   "
            f"one imaginary mode {imag[0]:.1f} cm\u207b\u00b9",
            transform=ax.transAxes, ha="center", va="top", fontsize=8.4,
            family="monospace",
        )
    handles = [
        Line2D([], [], marker="o", ls="", color=COLOUR[e], label=e,
               markersize=8, markeredgecolor="white")
        for e in ("C", "N", "O", "F", "H")
    ] + [
        Line2D([], [], color="#c0392b", lw=2, ls=(0, (2.5, 2)),
               label="forming bond"),
        Line2D([], [], color="#9aa0a6", lw=2, ls=(0, (2.5, 2)),
               label="cross contact \u2014 not a bond"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False,
               fontsize=9.2, bbox_to_anchor=(0.5, 0.004))
    fig.subplots_adjust(left=0.0, right=1.0, top=0.93, bottom=0.105,
                        wspace=0.0)
    fig.savefig(FIG / "figB3-saddles.png", dpi=210)
    plt.close(fig)


# --------------------------------------------------------------- figure 4
def figure_levels():
    levels = [
        ("B3LYP-D3BJ/def2-SVP\n(the opt+freq protocol)", 0.32573920,
         -0.45559011, 1.9145201),
        ("B3LYP-D3BJ/def2-TZVP\n(basis axis)", 0.55999321,
         -0.2213361, 1.3709801),
        ("DSD-BLYP-D3BJ/def2-TZVP\n(functional axis, delivered)", 1.4565436,
         0.67521428, 2.6183712),
    ]
    unc = 1.2323116
    floor_terms = [
        ("geometry convergence\n(measured)", 0.15886191, "#1b4f78"),
        ("RRHO vs Grimme qRRHO\n(measured)", 0.35644379, "#1b4f78"),
        ("conformer + neat medium\n(asserted, never sampled)", 0.75,
         "#bbbbbb"),
        ("functional + basis\n(measured in cycle 5)", 0.89655038, "#1b4f78"),
    ]

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(13.6, 5.0),
        gridspec_kw={"width_ratios": (1.12, 1.0)},
    )

    ys = np.arange(len(levels))[::-1]
    for y, (name, _elec, signed, ratio) in zip(ys, levels):
        colour = "#1b4f78" if signed < 0 else "#b8570e"
        ax.barh(y, signed, color=colour, height=0.46, edgecolor="white")
        side = ("ester-at-C4 lower" if signed < 0
                else "ester-at-C5 lower")
        ax.text(
            0.78, y,
            f"{signed:+.4f}    {side}    {ratio:.2f}:1",
            va="center", ha="left",
            fontsize=9.0, color=colour, fontweight="bold",
        )
    ax.axvline(0, color="black", lw=1.0)
    ax.axvspan(-0.5, 0.5, color="#888888", alpha=0.13, zorder=0)
    ax.text(0.0, len(levels) - 0.42,
            "the task's own decision band, ±0.5 kcal/mol",
            ha="center", va="bottom", fontsize=8.8, color="#555555")
    ax.annotate("", xy=(0.67521428, -0.55), xytext=(-0.45559011, -0.55),
                arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.4))
    ax.text(0.11, -0.70,
            "spread across levels 1.1308 kcal/mol, and it brackets zero:\n"
            "the DIRECTION of the preference is not established",
            ha="center", va="top", fontsize=8.8, color="#c0392b")
    ax.set_xlim(-0.80, 2.30)
    ax.set_ylim(-1.25, len(levels) + 0.08)
    ax.set_yticks(ys)
    ax.set_yticklabels([n for n, *_r in levels], fontsize=8.8)
    ax.set_xlabel(
        "signed $\\Delta\\Delta G^{\\ddag}$(353 K) = "
        "TS(ester-at-C4) $-$ TS(ester-at-C5)  /  kcal mol$^{-1}$",
        fontsize=9.6,
    )
    ax.set_title(
        "(a) one pair of geometries, three electronic-structure levels,\n"
        "identical 353 K RRHO thermal corrections — and the sign flips",
        fontsize=9.8,
    )
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ys = np.arange(len(floor_terms))[::-1]
    for y, (name, mag, colour) in zip(ys, floor_terms):
        bx.barh(y, mag, color=colour, height=0.5, edgecolor="white")
        bx.text(mag + 0.015, y, f"{mag:.4f}", va="center", fontsize=9.0)
    bx.set_yticks(ys)
    bx.set_yticklabels([n for n, *_r in floor_terms], fontsize=8.8)
    bx.axvline(0.5, color="black", lw=1.3, ls=":")
    bx.text(0.49, len(floor_terms) - 0.38, "task tolerance\n0.5 kcal/mol",
            ha="right", va="bottom", fontsize=8.8)
    bx.axvline(0.84545212, color="#c0392b", lw=1.7)
    bx.annotate(
        "0.8455 = RSS of the three non-electronic terms, i.e. the floor\n"
        "with the ELECTRONIC term set to zero: 0.3455 over the tolerance",
        xy=(0.84545212, -0.52), xytext=(1.70, -0.90),
        fontsize=8.8, color="#c0392b", ha="right", va="center",
        arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.2),
    )
    bx.axvline(unc, color="#7b2d8b", lw=1.7, ls="--")
    bx.text(unc + 0.02, len(floor_terms) - 0.38,
            f"delivered total\n{unc:.4f}", fontsize=8.8, color="#7b2d8b",
            ha="left", va="bottom")
    bx.set_xlim(0, 1.72)
    bx.set_ylim(-1.35, len(floor_terms) + 0.22)
    bx.set_xlabel("1$\\sigma$ half-width contributed  /  kcal mol$^{-1}$",
                  fontsize=9.6)
    bx.set_title(
        "(b) why the precision is unreachable, not merely unreached:\n"
        "the THREE non-electronic terms alone exceed the tolerance",
        fontsize=9.8,
    )
    for spine in ("top", "right", "left"):
        bx.spines[spine].set_visible(False)
    bx.tick_params(axis="y", length=0)

    fig.tight_layout()
    fig.savefig(FIG / "figB4-level-ladder.png", dpi=210)
    plt.close(fig)


if __name__ == "__main__":
    tables = {node: scan_table(node) for node, *_rest in SCANS}
    for node, rows in tables.items():
        print("##", node)
        print("  pt  driven   dE/kcal   N1-Cest N1-Ccf3 N3-Cest N3-Ccf3  C-C")
        for r in rows:
            print("  %2d  %6.3f  %8.3f   %7.3f %7.3f %7.3f %7.3f  %6.4f"
                  % (r["point"], r["coord"], r["kcal"], r["n1_est"],
                     r["n1_cf3"], r["n3_est"], r["n3_cf3"], r["cc"]))
    figure_scans(tables)
    figure_scan_structures(tables)
    figure_saddles()
    figure_levels()
    print("wrote:", sorted(p.name for p in FIG.glob("*.png")))
