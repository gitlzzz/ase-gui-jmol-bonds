#!/usr/bin/env python
"""Launch the REAL ase gui with Jmol-style bonds, patched IN MEMORY only.

Same result as patch_ase_bond_colors.py (half-colored Jmol bond sticks, the
rule d < R_i + R_j + 0.45 A, smaller balls) but writes NOTHING to disk: it
loads ase/gui/view.py source, applies the same edits in memory, exec's the
patched module back into place, monkey-patches ui.ASEGUIWindow.line, and then
opens the normal ASE GUI. The installed ASE files are never modified.

Re-patches every run, so it survives ASE upgrades and is fully portable --
copy this one file anywhere and run it with the python that has ase.

The VIEW_EDITS below encode the SAME edits as patch_ase_bond_colors.py.

Usage:
    python ase_gui_jmol_launcher.py STRUCTURE        # GUI with bonds shown
    python ase_gui_jmol_launcher.py --check           # verify bond rule, no GUI

Run with the env that has ase, e.g.
    /opt/anaconda3/envs/cathub/bin/python ase_gui_jmol_launcher.py STRUCTURE
"""
import sys

# (old_pristine, new, already_there_marker); old matches pristine ASE 3.21-3.26
VIEW_EDITS = [
    ("        self.X_bonds = self.X[natoms + len(B1):]",
     "        self.X_bonds = self.X[natoms + len(B1):]\n"
     "\n        self.bonds = bonds",
     "self.bonds = bonds"),
    ("        bond_linewidth = self.scale * 0.15",
     "        bond_linewidth = max(1, int(self.scale * 0.18))",
     "max(1, int(self.scale * 0.18))"),
    ("                else:\n"
     "                    line((X1[a, 0], X1[a, 1],\n"
     "                          X2[a, 0], X2[a, 1]),\n"
     "                         width=bond_linewidth)",
     "                else:\n"
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
     "                         capstyle='round')",
     "self.bonds[a - ncell]"),
    ("            r *= 0.65", "            r *= 0.5", "            r *= 0.5"),
    ("r = 0.65 * self.get_covalent_radii()",
     "r = 0.5 * self.get_covalent_radii()",
     "r = 0.5 * self.get_covalent_radii()"),
    ("covalent_radii * 1.5", "covalent_radii + 0.225", "covalent_radii + 0.225"),
    ("self.get_covalent_radii(atoms)",
     "self.images.covalent_radii[atoms.numbers]",
     "self.images.covalent_radii[atoms.numbers]"),
]


def patch_view_in_memory():
    import ase.gui.view as vmod
    src = open(vmod.__file__).read()
    for old, new, marker in VIEW_EDITS:
        if marker in src:
            continue
        if old in src:
            src = src.replace(old, new)
        else:
            print("WARN: view.py pattern not found, skipped:", repr(old[:40]))
    exec(compile(src, vmod.__file__, "exec"), vmod.__dict__)
    return vmod


def patch_ui_line():
    import ase.gui.ui as uimod

    def line(self, bbox, width=1, color="black", capstyle="butt"):
        self.canvas.create_line(*tuple(int(x) for x in bbox), width=width,
                                fill=color, capstyle=capstyle)

    uimod.ASEGUIWindow.line = line


def check(vmod):
    import numpy as np
    from ase import Atoms
    from ase.data import covalent_radii, atomic_numbers
    z = atomic_numbers["Ru"]
    rad = np.array([covalent_radii[z], covalent_radii[z]])
    for d in (3.30, 3.40):                # Jmol Ru-Ru cutoff is 3.37 A
        a = Atoms("Ru2", positions=[[0, 0, 0], [d, 0, 0]])
        n = len(vmod.get_bonds(a, rad))
        print(f"  Ru-Ru at {d:.2f} A -> {'BOND' if n else 'no bond'} "
              f"(nbonds={n}); expected {'BOND' if d < 3.37 else 'no bond'}")


def main():
    args = sys.argv[1:]
    do_check = "--check" in args
    files = [a for a in args if not a.startswith("-")]

    vmod = patch_view_in_memory()
    patch_ui_line()
    print("in-memory patch applied to", vmod.__file__)

    if do_check:
        check(vmod)
        print("check done (no GUI opened).")
        return

    if not files:
        sys.exit("usage: ase_gui_jmol_launcher.py STRUCTURE  (or --check)")

    from ase.gui.images import Images
    from ase.gui.gui import GUI
    images = Images()
    images.read([files[0]])
    gui = GUI(images, "", True, None)     # show_bonds=True at startup
    gui.run()


if __name__ == "__main__":
    main()
