#!/usr/bin/env python
"""Patch the ASE GUI so Ctrl+B bonds look and behave like Jmol.

What this patch changes (ase/gui/ui.py and ase/gui/view.py):
  1. Bonds are half-colored: each half of a bond takes the color of the atom
     on that end (e.g. red for the O half, gray for the Ru half), instead of a
     single black line.
  2. Each stick gets a thin black outline (a slightly wider black line is drawn
     underneath) and round end caps, for a Jmol ball-and-stick look.
  3. Bond line width is max(1, int(scale * 0.18)).
  4. Atom balls in bond mode shrink from 0.65x to 0.5x the covalent radius, with
     the bond start offset shrunk to match (so sticks meet the ball with no gap)
     -- closer to Jmol's ball-and-stick proportions.
  5. The bond decision is made Jmol-like (see "Bonding rule" below).

Bonding rule after patching (in get_bonds):
    atoms i and j are bonded  <=>  d(i,j) < R_i + R_j + 0.45 A
  - R = ASE Cordero (2008) covalent radii (ase.data.covalent_radii), UNSCALED:
    the decision reads self.images.covalent_radii directly, so it does NOT use
    images.atom_scale (default 0.89) and is independent of the atom-size slider.
  - +0.45 A additive tolerance, implemented as NeighborList(covalent_radii +
    0.225): the neighbor list bonds when d < cutoff_i + cutoff_j, so each atom
    carries half the tolerance (0.225 A).
  - strict less-than (<), skin = 0, self_interaction = False (an atom may still
    bond to its own periodic image). No lower bound (Jmol's 0.4 A
    minBondDistance is NOT applied).
  This replaces ASE's stock multiplicative rule d < 1.5*(R_i + R_j).
  Examples: Ru-Ru if d < 3.37 A; Ru-O if d < 2.57 A; O-O if d < 1.77 A.

  Remaining differences from stock Jmol (left in place on purpose):
  - radius TABLE: ASE uses Cordero; Jmol uses OpenBabel (default) or BODR 2014.
    For catalysis elements (O, Ru, Pt, Ir, Rh, Fe, C, N, H) these agree to
    within ~0.08 A; biggest gaps are alkali / alkaline-earth / noble gases.
  - no minBondDistance lower bound (only matters for near-coincident atoms).

Re-runnable and idempotent: it locates the installed ASE, makes a one-time
.orig backup of each file, and only edits if the patch is not already there.
Run it again after upgrading ASE or on a new machine.

Usage:
    python patch_ase_bond_colors.py            # patch the ASE of THIS python
    python patch_ase_bond_colors.py /path/to/site-packages/ase   # explicit
    python patch_ase_bond_colors.py --revert   # restore .orig backups

IMPORTANT: run it with the SAME interpreter that runs `ase gui`. To find it:
    head -1 $(which ase)
then call that python, e.g.
    /opt/anaconda3/envs/cathub/bin/python patch_ase_bond_colors.py
"""
import os
import shutil
import sys


# --- edits: each is (find_text, replace_text, already_applied_marker) ---------
# The find_text matches pristine upstream ASE (3.21-3.26). If upstream changes
# these code shapes in a future version, the corresponding edit will be skipped
# with a warning instead of corrupting the file.

UI_OLD_SIG = "    def line(self, bbox, width=1):"
UI_NEW_SIG = "    def line(self, bbox, width=1, color='black', capstyle='butt'):"

# the create_line body has two pristine variants (with/without fill='black')
UI_BODY_VARIANTS = [
    "        self.canvas.create_line(*tuple(int(x) for x in bbox), width=width)",
    "        self.canvas.create_line(*tuple(int(x) for x in bbox), width=width,\n"
    "                                fill='black')",
]
UI_BODY_NEW = (
    "        self.canvas.create_line(*tuple(int(x) for x in bbox), width=width,\n"
    "                                fill=color, capstyle=capstyle)"
)

VIEW_XBONDS_ANCHOR = "        self.X_bonds = self.X[natoms + len(B1):]"
VIEW_STORE_BONDS = (
    VIEW_XBONDS_ANCHOR
    + "\n\n        # keep bond atom indices so draw() can color each half by its atom"
      "\n        self.bonds = bonds"
)

VIEW_LW_OLD = "        bond_linewidth = self.scale * 0.15"
VIEW_LW_NEW = "        bond_linewidth = max(1, int(self.scale * 0.18))"

# Jmol-like ball-and-stick: shrink the atom ball shown in bond mode and the
# matching bond start offset from 0.65 to 0.5 of the covalent radius. Both
# must stay equal so sticks meet the ball surface with no gap.
VIEW_BALL_OLD = "            r *= 0.65"
VIEW_BALL_NEW = "            r *= 0.5"
VIEW_OFFSET_OLD = "r = 0.65 * self.get_covalent_radii()"
VIEW_OFFSET_NEW = "r = 0.5 * self.get_covalent_radii()"

# Jmol-style bond decision in get_bonds(): replace ASE's multiplicative
# d < 1.5*(R_i + R_j) with the additive d < (R_i + R_j) + 0.45 A. NeighborList
# sums per-atom radii, so +0.225 each gives a 0.45 A total tolerance.
VIEW_BOND_OLD = "covalent_radii * 1.5"
VIEW_BOND_NEW = "covalent_radii + 0.225"

# Decouple the bond decision from the display atom-size scale: feed full
# covalent radii (self.images.covalent_radii) instead of self.get_covalent_radii,
# which multiplies by images.atom_scale (default 0.89) and is driven by the
# atom-size slider. Keeps d < (R_i + R_j) + 0.45 A stable, like Jmol.
VIEW_RADII_OLD = "self.get_covalent_radii(atoms)"
VIEW_RADII_NEW = "self.images.covalent_radii[atoms.numbers]"

VIEW_DRAW_OLD = (
    "                else:\n"
    "                    line((X1[a, 0], X1[a, 1],\n"
    "                          X2[a, 0], X2[a, 1]),\n"
    "                         width=bond_linewidth)"
)
VIEW_DRAW_NEW = (
    "                else:\n"
    "                    # Split the bond at its midpoint so each half takes the\n"
    "                    # color of the atom on that end. A slightly wider black\n"
    "                    # underlay gives each stick a visible dark edge.\n"
    "                    bond = self.bonds[a - ncell]\n"
    "                    mx = (X1[a, 0] + X2[a, 0]) // 2\n"
    "                    my = (X1[a, 1] + X2[a, 1]) // 2\n"
    "                    line((X1[a, 0], X1[a, 1], X2[a, 0], X2[a, 1]),\n"
    "                         width=bond_linewidth + 2, color='black',\n"
    "                         capstyle='round')\n"
    "                    line((X1[a, 0], X1[a, 1], mx, my),\n"
    "                         width=bond_linewidth, color=colors[bond[0]],\n"
    "                         capstyle='round')\n"
    "                    line((mx, my, X2[a, 0], X2[a, 1]),\n"
    "                         width=bond_linewidth, color=colors[bond[1]],\n"
    "                         capstyle='round')"
)


def backup(path):
    orig = path + ".orig"
    if not os.path.exists(orig):
        shutil.copy2(path, orig)
        print("  backup ->", orig)


def revert(path):
    orig = path + ".orig"
    if os.path.exists(orig):
        shutil.copy2(orig, path)
        print("  reverted", path)
    else:
        print("  no backup for", path)


def patch_ui(path):
    src = open(path).read()
    if "capstyle=capstyle" in src:
        print("  ui.py already patched")
        return
    if UI_OLD_SIG not in src:
        print("  WARN: ui.py line() signature not found; skipped")
        return
    backup(path)
    src = src.replace(UI_OLD_SIG, UI_NEW_SIG)
    for variant in UI_BODY_VARIANTS:
        if variant in src:
            src = src.replace(variant, UI_BODY_NEW)
            break
    else:
        print("  WARN: ui.py create_line body not matched; skipped")
        return
    open(path, "w").write(src)
    print("  ui.py patched")


def patch_view(path):
    src = open(path).read()
    changed = False
    backed = False

    if "self.bonds = bonds" in src:
        print("  view.py: store-bonds already present")
    elif VIEW_XBONDS_ANCHOR in src:
        backup(path); backed = True
        src = src.replace(VIEW_XBONDS_ANCHOR, VIEW_STORE_BONDS, 1)
        changed = True
    else:
        print("  WARN: view.py X_bonds anchor not found; skipped store-bonds")

    if VIEW_LW_NEW in src:
        print("  view.py: linewidth already patched")
    elif VIEW_LW_OLD in src:
        if not backed:
            backup(path); backed = True
        src = src.replace(VIEW_LW_OLD, VIEW_LW_NEW)
        changed = True
    else:
        print("  WARN: view.py bond_linewidth not found; skipped")

    if "self.bonds[a - ncell]" in src:
        print("  view.py: draw block already patched")
    elif VIEW_DRAW_OLD in src:
        if not backed:
            backup(path); backed = True
        src = src.replace(VIEW_DRAW_OLD, VIEW_DRAW_NEW)
        changed = True
    else:
        print("  WARN: view.py draw block not matched; skipped")

    if VIEW_BALL_NEW in src:
        print("  view.py: ball radius already patched")
    elif VIEW_BALL_OLD in src:
        if not backed:
            backup(path); backed = True
        src = src.replace(VIEW_BALL_OLD, VIEW_BALL_NEW)
        changed = True
    else:
        print("  WARN: view.py ball radius (r *= 0.65) not found; skipped")

    if VIEW_OFFSET_NEW in src:
        print("  view.py: bond offset already patched")
    elif VIEW_OFFSET_OLD in src:
        if not backed:
            backup(path); backed = True
        src = src.replace(VIEW_OFFSET_OLD, VIEW_OFFSET_NEW)
        changed = True
    else:
        print("  WARN: view.py bond offset (0.65 * covalent) not found; skipped")

    if VIEW_BOND_NEW in src:
        print("  view.py: bond decision already patched")
    elif VIEW_BOND_OLD in src:
        if not backed:
            backup(path); backed = True
        src = src.replace(VIEW_BOND_OLD, VIEW_BOND_NEW)
        changed = True
    else:
        print("  WARN: view.py bond decision (covalent_radii * 1.5) not found; skipped")

    if VIEW_RADII_NEW in src:
        print("  view.py: bond radii decoupling already patched")
    elif VIEW_RADII_OLD in src:
        if not backed:
            backup(path); backed = True
        src = src.replace(VIEW_RADII_OLD, VIEW_RADII_NEW)
        changed = True
    else:
        print("  WARN: view.py bond radii call (self.get_covalent_radii(atoms)) not found; skipped")

    if changed:
        open(path, "w").write(src)
        print("  view.py patched")


def find_ase(args):
    for a in args:
        if not a.startswith("-"):
            return a
    import ase
    return os.path.dirname(ase.__file__)


def main():
    args = sys.argv[1:]
    ase_dir = find_ase(args)
    ui = os.path.join(ase_dir, "gui", "ui.py")
    view = os.path.join(ase_dir, "gui", "view.py")
    if not (os.path.exists(ui) and os.path.exists(view)):
        sys.exit("Could not find ase/gui/{ui,view}.py under " + ase_dir)

    print("ASE gui dir:", os.path.join(ase_dir, "gui"))
    if "--revert" in args:
        revert(ui)
        revert(view)
        print("Done (reverted). Recompile cache by reopening `ase gui`.")
        return

    patch_ui(ui)
    patch_view(view)
    print("Done. Reopen `ase gui` and press Ctrl+B.")


if __name__ == "__main__":
    main()
