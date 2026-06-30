#!/usr/bin/env python
"""Render the same structure two ways for comparison: stock ASE GUI bond style
vs. this patch. Saves comparison.png.

Left  - Original ASE: single black bond lines, larger balls, multiplicative
        bond rule d < 1.5*(R_i + R_j) (over-bonds metals).
Right - This patch: half-colored sticks, smaller balls, additive Jmol rule
        d < (R_i + R_j) + 0.45 A.

Usage:
    python compare_views.py [structure]   # default: a rutile RuO2 demo cell
"""
import os
import sys

sys.path.insert(0, os.path.expanduser("~/script"))
import plotstyle  # noqa: E402
plotstyle.apply()

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from ase.data import covalent_radii  # noqa: E402
from ase.data.colors import jmol_colors  # noqa: E402
from ase.neighborlist import NeighborList  # noqa: E402

ATOM_SCALE = 0.89          # ASE GUI default radii_scale


def load(argv):
    files = [a for a in argv[1:] if not a.startswith("-")]
    if files and os.path.exists(files[0]):
        from ase.io import read
        return read(files[0])
    from ase.spacegroup import crystal      # self-contained rutile RuO2 demo
    ruo2 = crystal(["Ru", "O"], basis=[(0, 0, 0), (0.3058, 0.3058, 0)],
                   spacegroup=136, cellpar=[4.4968, 4.4968, 3.1049, 90, 90, 90])
    return ruo2 * (2, 2, 3)


def rot(ax_deg, ay_deg):
    ax, ay = np.radians(ax_deg), np.radians(ay_deg)
    Rx = np.array([[1, 0, 0], [0, np.cos(ax), -np.sin(ax)], [0, np.sin(ax), np.cos(ax)]])
    Ry = np.array([[np.cos(ay), 0, np.sin(ay)], [0, 1, 0], [-np.sin(ay), 0, np.cos(ay)]])
    return Rx @ Ry


def bond_halves(atoms, cutoffs):
    """Each (atom i -> midpoint to neighbor), as ASE/GUI does."""
    nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True)
    nl.update(atoms)
    cell, pos = atoms.cell[:], atoms.positions
    halves = []
    for i in range(len(atoms)):
        idx, off = nl.get_neighbors(i)
        for j, o in zip(idx, off):
            end = pos[j] + o @ cell
            halves.append((i, pos[i], 0.5 * (pos[i] + end)))
    return halves


def draw_panel(ax, atoms, P, depth, title, sub, *, half_color, cutoffs, ball):
    colors = [jmol_colors[z] for z in atoms.numbers]
    rad = ball * covalent_radii[atoms.numbers]
    dn = (depth - depth.min()) / (np.ptp(depth) + 1e-9)

    for i, start, mid in bond_halves(atoms, cutoffs):
        s = start @ R.T
        m = mid @ R.T
        z = 1 + 0.5 * (dn[i] + (m[2] - depth.min()) / (np.ptp(depth) + 1e-9))
        if half_color:
            ax.plot([s[0], m[0]], [s[1], m[1]], color="black", lw=4.0,
                    solid_capstyle="round", zorder=z)
            ax.plot([s[0], m[0]], [s[1], m[1]], color=colors[i], lw=2.8,
                    solid_capstyle="round", zorder=z + 0.01)
        else:
            ax.plot([s[0], m[0]], [s[1], m[1]], color="black", lw=1.4,
                    solid_capstyle="round", zorder=z)

    for i in np.argsort(depth):
        ax.add_patch(Circle((P[i, 0], P[i, 1]), rad[i], facecolor=colors[i],
                            edgecolor="black", lw=0.4, zorder=10 + dn[i]))

    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10, pad=4)
    ax.text(0.5, -0.02, sub, transform=ax.transAxes, ha="center", va="top",
            fontsize=7.5)


R = rot(-72, 18)


def main():
    atoms = load(sys.argv)
    atoms.positions -= atoms.positions.mean(0)
    P = atoms.positions @ R.T
    depth = P[:, 2]
    cov = covalent_radii[atoms.numbers]

    # stock ASE: radii = covalent * atom_scale * 1.5  -> sum = 1.335*(Ci+Cj)
    cut_orig = 1.5 * ATOM_SCALE * cov
    # this patch: radii = covalent + 0.225 -> sum = Ci+Cj+0.45  (unscaled)
    cut_new = cov + 0.225

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(7.0, 4.6))
    draw_panel(a0, atoms, P, depth, "Original ASE",
               r"$d < 1.5\,(R_i{+}R_j)$  ·  black lines  ·  larger balls",
               half_color=False, cutoffs=cut_orig, ball=ATOM_SCALE * 0.65)
    draw_panel(a1, atoms, P, depth, "This patch (Jmol-style)",
               r"$d < R_i{+}R_j{+}0.45$ Å  ·  half-colored  ·  smaller balls",
               half_color=True, cutoffs=cut_new, ball=ATOM_SCALE * 0.5)

    fig.suptitle("ASE GUI bond view  (rutile RuO$_2$)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "comparison.png")
    fig.savefig(out, dpi=200)
    print("wrote", out, "| atoms", len(atoms))


if __name__ == "__main__":
    main()
