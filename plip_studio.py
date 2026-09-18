#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLIP Studio v2
==============

Protein-ligand interaction analysis with GUARANTEED PyMOL session export.

Pipeline:
  1. Input: protein + ligand (separate files) OR one complex PDB.
  2. Ligand converted to PDB via OpenBabel if needed.
  3. Complex normalized (ligand forced to HETATM, waters/buffers handled).
  4. PLIP analyses the complex  ->  XML report.
  5. A dedicated headless PyMOL script reads the XML and builds a .pse
     showing the ligand + binding site + every interaction as coloured
     dashed lines.  This bypasses PLIP's own PyMOL visualizer which
     crashes on many docking outputs.

Icon:
  Put a file named `icon.png` next to this script.  It is used as the
  window icon, the taskbar icon (on Windows via AppUserModelID) and the
  icon for every dialog and crash window.

Output (in a fresh per-run folder):
    <name>_complex.pdb
    <name>.pse                       <-- open in PyMOL
    <name>.xml / <name>.txt          <-- PLIP reports
    <name>_interactions.csv
    <name>_binding_site_residues.csv
    <name>_ligand_properties.csv
    run_summary.txt
    PLIPStudio_LOG.txt

Requirements (same interpreter that runs this script):
    pip install plip openbabel lxml pymol-open-source
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
import traceback
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------

APP_NAME = "PLIP Studio"
VERSION  = "2.0.1"
APP_ID   = "PLIPStudio.App.2"     # AppUserModelID (Windows taskbar grouping)
ICON_NAME = "icon.png"

LIGAND_EXTS  = (".pdb", ".ent", ".pdbqt", ".mol2", ".sdf", ".mol")
PROTEIN_EXTS = (".pdb", ".ent", ".cif", ".mmcif")

WATER_CODES = {"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3"}

COMMON_BUFFERS = {
    "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "P6G", "DMS", "TRS", "MES",
    "EPE", "CIT", "ACT", "ACY", "FMT", "IMD", "BME", "DTT", "URE", "MPD",
    "IOD", "CL", "BR", "NA", "K", "MG", "CA", "ZN", "FE", "FE2", "CU",
    "MN", "CO", "NI", "CD", "BA", "SR", "LI", "CS", "RB", "ACE", "NH4",
    "NO3", "SCN", "TAR", "MLI", "BEN", "EOH", "MOH", "IPA", "BU1", "CAC",
}

STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
    "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC","PYL",
}
STANDARD_NT = {"DA","DT","DC","DG","DU","A","U","C","G","I"}

INTERACTION_TYPES = [
    ("hydrophobic_interactions",  "Hydrophobic"),
    ("hydrogen_bonds",            "HydrogenBond"),
    ("water_bridges",             "WaterBridge"),
    ("salt_bridges",              "SaltBridge"),
    ("pi_stacks",                 "PiStacking"),
    ("pi_cation_interactions",    "PiCation"),
    ("halogen_bonds",             "HalogenBond"),
    ("metal_complexes",           "MetalComplex"),
]

CRASH_LOG = "PLIPSTUDIO_CRASH.txt"

# ---------------------------------------------------------------------------
#  Icon helpers  (window + taskbar icon from icon.png next to this script)
# ---------------------------------------------------------------------------

_ICON_REF = None   # keep a strong reference so Tk cannot GC the PhotoImage


def _icon_path() -> Optional[Path]:
    """Locate icon.png next to this script (or in cwd as fallback)."""
    here = Path(__file__).resolve().parent
    for cand in (here / ICON_NAME, Path.cwd() / ICON_NAME):
        if cand.is_file():
            return cand
    return None


def _set_windows_appid(appid: str = APP_ID) -> None:
    """
    Windows groups taskbar buttons by AppUserModelID.  Calling this BEFORE
    creating the first tk.Tk() makes the OS use our icon.png for the taskbar
    button instead of the generic python.exe icon.
    """
    if sys.platform != "win32":
        return
    try:
        import ctypes
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(appid)
    except Exception:
        pass


def _apply_icon(window) -> None:
    """
    Set the window icon from icon.png.  Tk >= 8.6 supports PNG directly.
    Keeps the PhotoImage alive in a module-level global so that it is not
    garbage-collected before the window is shown.  Falls back to a sibling
    icon.ico on Windows if PNG loading fails for any reason.
    """
    global _ICON_REF
    p = _icon_path()
    if p is None:
        return
    try:
        img = tk.PhotoImage(file=str(p))
        _ICON_REF = img                     # <- critical: keep a reference
        window.iconphoto(True, img)         # True = default for all Toplevels
        return
    except Exception:
        pass
    ico = p.with_suffix(".ico")
    if ico.is_file():
        try:
            window.iconbitmap(default=str(ico))
        except Exception:
            pass


# ---------------------------------------------------------------------------
#  Headless PyMOL script that builds the .pse from PLIP's XML
# ---------------------------------------------------------------------------

PSE_SCRIPT = r'''# -*- coding: utf-8 -*-
"""Build a PLIP-style PyMOL session from a complex PDB + PLIP XML report."""
import json
import sys
import xml.etree.ElementTree as ET


def _coo(elem):
    return [float(elem.find(t).text) for t in ("x", "y", "z")]


def main(cfg_path):
    cfg = json.load(open(cfg_path, encoding="utf-8"))

    import pymol
    pymol.finish_launching(["pymol", "-c", "-q", "-i"])
    from pymol import cmd

    cmd.load(cfg["complex"], "complex")
    cmd.remove("solvent")
    cmd.hide("everything", "all")

    cmd.set("dash_gap",    0.20)
    cmd.set("dash_width",  2.8)
    cmd.set("dash_length", 0.40)
    cmd.set("dash_radius", 0.06)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("cartoon_transparency",  0.75)

    STYLE = {
        "hydrogen_bonds":           ("HBond",        "blue"),
        "hydrophobic_interactions": ("Hydrophobic",  "gray50"),
        "halogen_bonds":            ("HalogenBond",  "teal"),
        "pi_stacks":                ("PiStack",      "magenta"),
        "pi_cation_interactions":   ("PiCation",     "orange"),
        "salt_bridges":             ("SaltBridge",   "yellow"),
        "water_bridges":            ("WaterBridge",  "lightblue"),
        "metal_complexes":          ("MetalComplex", "violetpurple"),
    }

    root = ET.parse(cfg["xml"]).getroot()

    dash_objs, ref_objs = [], []
    counter = 0

    for bs in root.findall("bindingsite"):
        ident  = bs.find("identifiers")
        hetid  = (ident.findtext("hetid")    or "").strip()
        chain  = (ident.findtext("chain")    or "").strip()
        pos    = (ident.findtext("position") or "").strip()
        if not hetid or not pos:
            continue

        sel = f"resn {hetid} and resi {pos}"
        if chain:
            sel_chain = f"resn {hetid} and chain {chain} and resi {pos}"
            try:
                if cmd.count_atoms(sel_chain) > 0:
                    sel = sel_chain
            except Exception:
                pass
        cmd.select("ligand", sel)
        cmd.show("sticks", "ligand")
        cmd.color("orange", "ligand and elem C")
        for el, col in (("O","red"),("N","blue"),("S","yellow"),("P","orange"),
                        ("F","lightblue"),("CL","green"),("BR","brown"),
                        ("I","purple"),("H","white")):
            cmd.color(col, f"ligand and elem {el}")

        site_sels = []
        for res in bs.findall("./bs_residues/bs_residue"):
            if res.get("contact") != "True":
                continue
            txt = (res.text or "").strip()
            if len(txt) < 2:
                continue
            num, ch = txt[:-1], txt[-1]
            site_sels.append(f"(chain {ch} and resi {num})")
        if site_sels:
            site_sel = " or ".join(site_sels)
        else:
            site_sel = "byres (ligand around 4.5) and not ligand"
        cmd.select("bsite", site_sel)
        cmd.show("sticks", "bsite and not ligand")
        cmd.color("cyan", "bsite and elem C")
        for el, col in (("O","red"),("N","blue"),("S","yellow")):
            cmd.color(col, f"bsite and elem {el}")

        inter = bs.find("interactions")
        if inter is None:
            continue

        for tag, (label, color) in STYLE.items():
            grp = inter.find(tag)
            if grp is None:
                continue
            for item in grp:
                lig   = item.find("ligcoo")
                prot  = item.find("protcoo")
                water = item.find("watercoo")
                if lig is None or prot is None:
                    continue

                la = _coo(lig)
                pa = _coo(prot)
                counter += 1

                if water is not None:
                    wa = _coo(water)
                    n1, n2 = f"r{counter}a", f"r{counter}b"
                    n3, n4 = f"r{counter}c", f"r{counter}d"
                    for nm, p in ((n1, la), (n2, wa), (n3, wa), (n4, pa)):
                        cmd.pseudoatom(nm, pos=p)
                    o1 = f"{label}_{counter:03d}_1"
                    o2 = f"{label}_{counter:03d}_2"
                    cmd.distance(o1, n1, n2)
                    cmd.distance(o2, n3, n4)
                    cmd.color(color, o1)
                    cmd.color(color, o2)
                    cmd.hide("everything", f"{n1} or {n2} or {n3} or {n4}")
                    dash_objs += [o1, o2]
                    ref_objs  += [n1, n2, n3, n4]
                else:
                    n1, n2 = f"r{counter}a", f"r{counter}b"
                    cmd.pseudoatom(n1, pos=la)
                    cmd.pseudoatom(n2, pos=pa)
                    obj = f"{label}_{counter:03d}"
                    cmd.distance(obj, n1, n2)
                    cmd.color(color, obj)
                    cmd.hide("everything", f"{n1} or {n2}")
                    dash_objs.append(obj)
                    ref_objs += [n1, n2]

    if dash_objs:
        cmd.group("Interactions", " ".join(dash_objs))
    if ref_objs:
        cmd.group("RefPoints", " ".join(ref_objs))
        cmd.disable("RefPoints")

    cmd.hide("labels", "all")
    cmd.center("ligand")
    cmd.zoom("ligand", 6)

    cmd.save(cfg["out"])
    print("PSE_SAVED " + cfg["out"])


if __name__ == "__main__":
    main(sys.argv[1])
'''


# ---------------------------------------------------------------------------
#  Small helpers
# ---------------------------------------------------------------------------

def now_ts() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def human_time() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def clock() -> str:
    return datetime.now().strftime("%H:%M:%S")


def _python_exe() -> str:
    """Return python.exe (not pythonw.exe) so subprocesses have a console."""
    exe = sys.executable
    if exe.lower().endswith("pythonw.exe"):
        cand = exe[:-len("w.exe")] + ".exe"
        if os.path.exists(cand):
            return cand
    return exe


# ---------------------------------------------------------------------------
#  PDB parsing
# ---------------------------------------------------------------------------

@dataclass
class HetGroup:
    code: str
    chain: str
    resnum: int
    natoms: int = 0
    nheavy: int = 0
    has_h: bool = False

    @property
    def label(self) -> str:
        return f"{self.code}:{self.chain or '-'}:{self.resnum}"


def _parse_atom(line: str) -> Optional[dict]:
    if len(line) < 54 or not line.startswith(("ATOM", "HETATM")):
        return None
    rec = {
        "record":  line[0:6].strip(),
        "serial":  line[6:11].strip(),
        "name":    line[12:16].strip(),
        "altloc":  line[16:17],
        "resname": line[17:20].strip(),
        "chain":   line[21:22].strip(),
        "resseq":  line[22:26].strip(),
        "x":       line[30:38].strip(),
        "y":       line[38:46].strip(),
        "z":       line[46:54].strip(),
        "element": line[76:78].strip() if len(line) >= 78 else "",
    }
    try:
        rec["resnum"] = int(rec["resseq"])
    except ValueError:
        return None
    return rec


def _element(rec: dict) -> str:
    if rec["element"]:
        return rec["element"].upper()
    m = re.match(r"([A-Za-z]{1,2})", rec["name"])
    return m.group(1).upper() if m else ""


def scan_het_groups(path: Path) -> List[HetGroup]:
    groups: Dict[Tuple[str, str, int], HetGroup] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = _parse_atom(line)
            if rec is None or rec["record"] != "HETATM":
                continue
            key = (rec["resname"], rec["chain"], rec["resnum"])
            g = groups.setdefault(key, HetGroup(*key))
            g.natoms += 1
            if _element(rec) == "H":
                g.has_h = True
            else:
                g.nheavy += 1
    return list(groups.values())


def detect_atom_ligands(path: Path) -> List[str]:
    codes = set()
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = _parse_atom(line)
            if rec and rec["record"] == "ATOM":
                c = rec["resname"].upper()
                if c not in STANDARD_AA and c not in STANDARD_NT and c not in WATER_CODES:
                    codes.add(c)
    return sorted(codes)


# ---------------------------------------------------------------------------
#  Ligand conversion
# ---------------------------------------------------------------------------

def convert_ligand_to_pdb(src: Path, dst: Path, add_h_ph: Optional[float],
                          log: Callable[[str, str], None]) -> None:
    infix = src.suffix.lower().lstrip(".")
    if infix == "ent":
        infix = "pdb"

    try:
        try:
            from openbabel import openbabel as ob
        except ImportError:
            import openbabel as ob  # type: ignore
        conv = ob.OBConversion()
        if not conv.SetInAndOutFormats(infix, "pdb"):
            raise ValueError(f"OpenBabel does not know format '{infix}'")
        mol = ob.OBMol()
        if not conv.ReadFile(mol, str(src)):
            raise ValueError(f"Cannot read ligand file: {src}")
        if add_h_ph:
            op = ob.OBOp.FindType("addpolarH")
            if op is not None:
                op.Do(mol, str(add_h_ph))
                log(f"Added polar hydrogens at pH {add_h_ph}: {mol.NumAtoms()} atoms", "INFO")
        conv.WriteFile(mol, str(dst))
        log(f"Ligand converted to PDB: {dst.name}", "OK")
        return
    except ImportError:
        pass

    exe = shutil.which("obabel") or shutil.which("obabel.exe")
    if not exe:
        raise RuntimeError("OpenBabel not found. Run: pip install openbabel")
    cmd = [exe, str(src), "-O", str(dst)]
    if add_h_ph:
        cmd += ["-p", str(add_h_ph)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    log(f"Ligand converted via obabel CLI: {dst.name}", "OK")


def _has_hydrogens(path: Path) -> bool:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = _parse_atom(line)
            if rec and _element(rec) == "H":
                return True
    return False


# ---------------------------------------------------------------------------
#  Complex normalization
# ---------------------------------------------------------------------------

def normalize_complex(src: Path, dst: Path, remove_waters: bool,
                      log: Callable[[str, str], None]) -> List[HetGroup]:
    out_lines: List[str] = []
    groups: Dict[Tuple[str, str, int], HetGroup] = {}
    n_fix = n_water = n_buf = 0

    with open(src, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            rec = _parse_atom(line)
            if rec is None:
                continue
            code = rec["resname"].upper()

            if code in WATER_CODES:
                if remove_waters:
                    n_water += 1
                    continue
                out_lines.append(line)
                continue

            if remove_waters and code in COMMON_BUFFERS:
                n_buf += 1
                continue

            newline = line
            if code not in STANDARD_AA and code not in STANDARD_NT:
                if rec["record"] == "ATOM":
                    newline = "HETATM" + line[6:]
                    n_fix += 1
                key = (rec["resname"], rec["chain"], rec["resnum"])
                g = groups.setdefault(key, HetGroup(*key))
                g.natoms += 1
                if _element(rec) == "H":
                    g.has_h = True
                else:
                    g.nheavy += 1
            out_lines.append(newline)

    out_lines.append("END")
    dst.write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    if n_fix:
        log(f"Rewrote {n_fix} ligand records from ATOM to HETATM", "OK")
    if n_water or n_buf:
        log(f"Removed {n_water} waters and {n_buf} buffer/ion records", "INFO")
    log(f"Complex normalized: {dst.name}", "OK")
    return list(groups.values())


def build_complex(protein: Path, ligand: Path, out: Path,
                  remove_waters: bool, log: Callable[[str, str], None]) -> List[HetGroup]:
    prot_lines: List[str] = []
    max_resnum: Dict[str, int] = {}
    max_serial = 0
    n_water = n_buf = 0

    with open(protein, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            rec = _parse_atom(line)
            if rec is None:
                continue
            code  = rec["resname"].upper()
            chain = rec["chain"] or "A"
            max_resnum[chain] = max(max_resnum.get(chain, 0), rec["resnum"])
            if rec["serial"].isdigit():
                max_serial = max(max_serial, int(rec["serial"]))

            if rec["record"] == "ATOM":
                prot_lines.append(line)
                continue
            if code in WATER_CODES:
                if remove_waters:
                    n_water += 1
                    continue
                prot_lines.append(line)
                continue
            if remove_waters and code in COMMON_BUFFERS:
                n_buf += 1
                continue
            prot_lines.append(line)

    lig_recs: List[dict] = []
    with open(ligand, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            rec = _parse_atom(line)
            if rec is not None:
                lig_recs.append(rec)
    if not lig_recs:
        raise RuntimeError("Ligand file contains no ATOM/HETATM records.")

    lig_chains = sorted({r["chain"] for r in lig_recs if r["chain"]})
    target_chain = lig_chains[0] if lig_chains else "L"
    if target_chain in max_resnum and max_resnum[target_chain] + len(lig_recs) > 9999:
        for cand in "LXYZ":
            if cand not in max_resnum:
                target_chain = cand
                break

    start = max_resnum.get(target_chain, 0) + 1
    if start > 9999:
        for cand in "LXYZ":
            if cand not in max_resnum:
                target_chain, start = cand, 1
                break
        else:
            raise RuntimeError("No free residue numbering space for the ligand.")

    remap: Dict[Tuple[str, str, int], int] = {}
    nxt = start
    for r in lig_recs:
        key = (r["resname"], r["chain"], r["resnum"])
        if key not in remap:
            remap[key] = nxt
            nxt += 1

    out_lines = list(prot_lines)
    lig_max = max([int(r["serial"]) for r in lig_recs if r["serial"].isdigit()] or [0])
    serial = max(max_serial, lig_max) + 1

    written: Dict[Tuple[str, str, int], HetGroup] = {}
    for r in lig_recs:
        newnum = remap[(r["resname"], r["chain"], r["resnum"])]
        el = _element(r)
        name = r["name"] if len(r["name"]) == 4 else " " + r["name"]
        line = (
            f"HETATM{serial:>5} {name:<4.4s}{r['altloc'] or ' '}"
            f"{r['resname']:>3s} {target_chain}{newnum:>4d}    "
            f"{r['x']:>8s}{r['y']:>8s}{r['z']:>8s}"
            f"{1.0:>6.2f}{0.0:>6.2f}          {el:>2s}"
        )
        out_lines.append(line)
        g = written.setdefault((r["resname"], target_chain, newnum),
                               HetGroup(r["resname"], target_chain, newnum))
        g.natoms += 1
        if el != "H":
            g.nheavy += 1
        else:
            g.has_h = True
        serial += 1

    out_lines.append("END")
    out.write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    extra = f", removed {n_water} waters and {n_buf} buffers" if (n_water or n_buf) else ""
    log(f"Built complex: {out.name} ({len(prot_lines)} protein + "
        f"{len(lig_recs)} ligand records{extra})", "OK")
    log("Ligands in complex: " + ", ".join(g.label for g in written.values()), "INFO")
    return list(written.values())


# ---------------------------------------------------------------------------
#  PLIP XML -> CSV tables
# ---------------------------------------------------------------------------

def _flatten(elem) -> Dict[str, str]:
    row: Dict[str, str] = {}
    for child in elem:
        if len(child):
            row[child.tag] = ";".join((c.text or "").strip() for c in child)
        else:
            row[child.tag] = (child.text or "").strip()
    return row


def xml_to_tables(xml_path: Path, runname: str,
                  log: Callable[[str, str], None]) -> Dict[str, object]:
    import xml.etree.ElementTree as ET

    root = ET.parse(str(xml_path)).getroot()
    sites: List[dict] = []
    inter_rows: List[Dict[str, str]] = []
    bs_rows: List[Dict[str, str]] = []
    lig_rows: List[Dict[str, str]] = []

    for bs in root.findall("bindingsite"):
        ident = bs.find("identifiers")

        def get(tag: str) -> str:
            return (ident.findtext(tag) or "").strip() if ident is not None else ""

        site = {
            "hetid":    get("hetid"),
            "chain":    get("chain"),
            "position": get("position"),
            "longname": get("longname"),
            "smiles":   get("smiles"),
            "inchikey": get("inchikey"),
            "counts":   {},
        }
        props = bs.find("lig_properties")
        prop_map = {p.tag: (p.text or "").strip() for p in props} if props is not None else {}
        site.update(prop_map)
        lig_rows.append({
            "Site":    site["longname"] or site["hetid"],
            "HETID":   site["hetid"],
            "Chain":   site["chain"],
            "ResNum":  site["position"],
            "SMILES":  site["smiles"],
            "InChIKey":site["inchikey"],
            **prop_map,
        })

        for res in bs.findall("./bs_residues/bs_residue"):
            bs_rows.append({
                "Site":        site["longname"] or site["hetid"],
                "ResNum":      (res.text or "").strip(),
                "AA":          res.get("aa", ""),
                "Interacting": res.get("contact", ""),
                "MinDist_A":   res.get("min_dist", ""),
            })

        inter = bs.find("interactions")
        if inter is not None:
            for tag, pretty in INTERACTION_TYPES:
                grp = inter.find(tag)
                items = list(grp) if grp is not None else []
                site["counts"][pretty] = len(items)
                for it in items:
                    row = _flatten(it)
                    for coo in ("ligcoo", "protcoo", "watercoo"):
                        row.pop(coo, None)
                    inter_rows.append({
                        "Site":            site["longname"] or site["hetid"],
                        "InteractionType": pretty,
                        **row,
                    })
        sites.append(site)

    def write_csv(rows, fname):
        if not rows:
            return None
        cols: List[str] = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
        p = xml_path.parent / fname
        with open(p, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        return p

    for rows, fname in ((inter_rows, f"{runname}_interactions.csv"),
                        (bs_rows,    f"{runname}_binding_site_residues.csv"),
                        (lig_rows,   f"{runname}_ligand_properties.csv")):
        f = write_csv(rows, fname)
        if f:
            log(f"CSV written: {f.name}", "OK")
    return {"sites": sites, "n_interactions": len(inter_rows)}


# ---------------------------------------------------------------------------
#  Run log
# ---------------------------------------------------------------------------

class RunLog:
    def __init__(self, gui_log: Callable[[str, str], None]):
        self._gui = gui_log
        self._fh = None
        self._path: Optional[Path] = None
        self._buffer: List[str] = []

    @property
    def path(self): return self._path

    def attach(self, run_dir: Path) -> Path:
        self._path = Path(run_dir) / "PLIPStudio_LOG.txt"
        self._fh = open(self._path, "a", encoding="utf-8")
        for line in self._buffer:
            self._fh.write(line + "\n")
        self._buffer = []
        self._fh.flush()
        return self._path

    def __call__(self, msg: str, level: str = "INFO"):
        line = f"[{clock()}] [{level:<4}] {msg}"
        if self._fh is not None:
            self._fh.write(line + "\n")
            self._fh.flush()
        else:
            self._buffer.append(line)
        self._gui(msg, level)

    def write_error(self, txt: str) -> Path:
        if self._fh is None:
            self._path = Path(__file__).resolve().parent / "PLIPStudio_LOG.txt"
            self._fh = open(self._path, "a", encoding="utf-8")
            for line in self._buffer:
                self._fh.write(line + "\n")
            self._buffer = []
        self._fh.write("\n" + "=" * 70 + f"\nERROR at {human_time()}\n" + txt +
                       "=" * 70 + "\n")
        self._fh.flush()
        return self._path

    def close(self):
        if self._fh is not None:
            try: self._fh.close()
            except OSError: pass
            self._fh = None


# ---------------------------------------------------------------------------
#  Pipeline
# ---------------------------------------------------------------------------

@dataclass
class Options:
    protein: Optional[Path] = None
    ligand:  Optional[Path] = None
    complex: Optional[Path] = None
    out_base: Optional[Path] = None
    runname: str = ""
    remove_waters: bool = True
    add_hydrogens: bool = True
    threads: int = 1
    keepmod: bool = False
    dnareceptor: bool = False


class Pipeline:
    def __init__(self, opts: Options, runlog: RunLog,
                 status: Callable[[str], None],
                 progress: Callable[[float], None],
                 cancel: threading.Event):
        self.opts = opts
        self.log = runlog
        self.status = status
        self.progress = progress
        self.cancel = cancel
        self.run_dir: Optional[Path] = None
        self.results: List[Path] = []
        self.summary: Optional[dict] = None
        self._proc: Optional[subprocess.Popen] = None

    def _check_cancel(self):
        if self.cancel.is_set():
            raise InterruptedError("Cancelled")

    def _run(self, cmd: List[str]) -> int:
        env = dict(os.environ)
        env["PYTHONIOENCODING"] = "utf-8"
        self.log("CMD> " + " ".join(f'"{c}"' if " " in c else c for c in cmd), "CMD")
        flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
        proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True,
                                encoding="utf-8", errors="replace",
                                creationflags=flags)
        self._proc = proc
        assert proc.stdout is not None
        for line in proc.stdout:
            self.log(line.rstrip("\n"), "OUT")
            if self.cancel.is_set():
                proc.kill()
                raise InterruptedError("Cancelled")
        proc.wait()
        return proc.returncode

    def run(self) -> Path:
        o = self.opts
        if not o.out_base:
            raise RuntimeError("No output folder selected.")

        base = Path(o.out_base)
        runname = o.runname.strip() or f"run_{now_ts()}"
        run_dir = base / runname
        n = 2
        while run_dir.exists():
            run_dir = base / f"{runname}_{n}"
            n += 1
        run_dir.mkdir(parents=True, exist_ok=True)
        self.run_dir = run_dir

        self.log.attach(run_dir)
        self.log(f"{APP_NAME} v{VERSION} started at {human_time()}", "INFO")
        self.log(f"python {sys.version.split()[0]} | exe: {_python_exe()}", "INFO")
        self.log(f"Output folder: {run_dir}", "INFO")

        # ---- 1. build complex ------------------------------------------
        self._check_cancel()
        self.status("Step 1/3: preparing complex ...")
        self.progress(0.05)

        if o.complex:
            raw = run_dir / f"{runname}_input_original.pdb"
            shutil.copy2(o.complex, raw)
            complex_pdb = run_dir / f"{runname}_complex.pdb"
            groups = normalize_complex(raw, complex_pdb, o.remove_waters, self.log)
            if not groups:
                raise RuntimeError(
                    "No non-standard residue (ligand) found in the complex file.\n"
                    "Make sure the ligand is stored as HETATM with a code like LIG.")
        else:
            if not o.protein or not o.ligand:
                raise RuntimeError("Both protein and ligand files are required.")
            lig_pdb = o.ligand
            if o.ligand.suffix.lower() not in (".pdb", ".ent"):
                lig_pdb = run_dir / f"{o.ligand.stem}_converted.pdb"
                ph = 7.4 if (o.add_hydrogens and not _has_hydrogens(o.ligand)) else None
                convert_ligand_to_pdb(o.ligand, lig_pdb, ph, self.log)
            elif o.add_hydrogens and not _has_hydrogens(o.ligand):
                lig_pdb = run_dir / f"{o.ligand.stem}_protonated.pdb"
                convert_ligand_to_pdb(o.ligand, lig_pdb, 7.4, self.log)
            complex_pdb = run_dir / f"{runname}_complex.pdb"
            build_complex(o.protein, lig_pdb, complex_pdb, o.remove_waters, self.log)

        # ---- 2. run PLIP for analysis ----------------------------------
        self._check_cancel()
        self.status("Step 2/3: running PLIP analysis ...")
        self.progress(0.30)

        plip_exe = shutil.which("plip") or shutil.which("plip.exe")
        if plip_exe:
            cmd = [plip_exe]
        else:
            cmd = [_python_exe(), "-m", "plip.plipcmd"]

        cmd += ["-f", str(complex_pdb),
                "-o", str(run_dir) + os.sep,
                "--name", runname,
                "-x", "-t",
                "--maxthreads", str(max(1, o.threads))]
        if not o.add_hydrogens: cmd.append("--nohydro")
        if o.keepmod:           cmd.append("--keepmod")
        if o.dnareceptor:       cmd.append("--dnareceptor")

        rc = self._run(cmd)
        if rc != 0:
            raise RuntimeError(f"PLIP failed with exit code {rc}. See log.")

        xmls = sorted(run_dir.glob(f"{runname}.xml")) or sorted(run_dir.glob("*.xml"))
        if not xmls:
            raise RuntimeError("PLIP did not produce an XML report.")
        xml_path = xmls[0]
        self.log(f"PLIP XML report: {xml_path.name}", "OK")

        # ---- 3. build tables + PSE -------------------------------------
        self._check_cancel()
        self.status("Step 3/3: building PyMOL session and tables ...")
        self.progress(0.70)

        summary = xml_to_tables(xml_path, runname, self.log)
        self.summary = summary
        if summary["n_interactions"] == 0:
            self.log("PLIP found 0 interactions.", "WARN")
        else:
            for s in summary["sites"]:
                parts = ", ".join(f"{k}: {v}" for k, v in s["counts"].items() if v)
                self.log(f"Site {s['longname'] or s['hetid']}: {parts}", "OK")

        pse_path = run_dir / f"{runname}.pse"
        self._build_pse(complex_pdb, xml_path, pse_path)

        for pattern in ("plipfixed.*", "plipfixed*"):
            for junk in run_dir.glob(pattern):
                try: junk.unlink()
                except OSError: pass

        self._write_summary(run_dir, runname, complex_pdb, pse_path, summary)
        self.results = sorted(p for p in run_dir.iterdir() if p.is_file())
        self.progress(1.0)
        self.log(f"Done. {len(self.results)} files in {run_dir}", "OK")
        return run_dir

    def _build_pse(self, complex_pdb: Path, xml_path: Path, out_pse: Path) -> None:
        self.log("Building PyMOL session (.pse) from PLIP XML ...", "INFO")
        script = self.run_dir / "_build_pse.py"           # type: ignore[union-attr]
        script.write_text(PSE_SCRIPT, encoding="utf-8")

        cfg = self.run_dir / "_pse_cfg.json"              # type: ignore[union-attr]
        cfg.write_text(json.dumps({
            "complex": str(complex_pdb),
            "xml":     str(xml_path),
            "out":     str(out_pse),
        }, ensure_ascii=False), encoding="utf-8")

        rc = self._run([_python_exe(), str(script), str(cfg)])

        for f in (script, cfg):
            try: f.unlink()
            except OSError: pass

        if rc == 0 and out_pse.exists():
            size = out_pse.stat().st_size
            self.log(f"PyMOL session created: {out_pse.name} ({size:,} bytes)", "OK")
        else:
            raise RuntimeError(
                "Failed to build the PyMOL session (.pse).\n"
                "Check that 'pymol' is installed in this Python:\n"
                f"    {_python_exe()} -m pip install pymol-open-source\n"
                "The XML report and CSV tables were still produced.")

    def _write_summary(self, run_dir: Path, runname: str,
                       complex_pdb: Path, pse_path: Path,
                       summary: Optional[dict]) -> None:
        o = self.opts
        lines = [
            "=" * 78,
            f"{APP_NAME} v{VERSION} - run summary",
            f"Date      : {human_time()}",
            f"Run name  : {runname}",
            f"Complex   : {complex_pdb.name}",
            f"Protein   : {o.protein or o.complex}",
            f"Ligand    : {o.ligand or '(inside complex file)'}",
            f"PyMOL     : {pse_path.name}",
            f"PLIP opts : threads={o.threads} keepmod={o.keepmod} "
            f"dnareceptor={o.dnareceptor} addH={o.add_hydrogens} "
            f"remove_waters={o.remove_waters}",
            "-" * 78,
        ]
        if summary:
            for s in summary["sites"]:
                lines.append(f"Ligand site: {s['longname'] or s['hetid']} "
                             f"({s['hetid']}:{s['chain']}:{s['position']})")
                lines.append(f"  SMILES: {s['smiles']}")
                lines += [f"  {k:<16s}: {v}" for k, v in s["counts"].items()]
        try:
            from plip.basic.config import __citation_information__ as cite  # type: ignore
            lines += ["-" * 78, "Please cite PLIP:", cite]
        except Exception:
            pass
        lines += ["-" * 78, "Files in this run:"]
        lines += [f"  {p.name}  ({p.stat().st_size:,} bytes)"
                  for p in sorted(run_dir.iterdir()) if p.is_file()]
        lines.append("=" * 78)
        (run_dir / "run_summary.txt").write_text("\n".join(lines), encoding="utf-8")
        self.log("Summary written: run_summary.txt", "OK")


# ---------------------------------------------------------------------------
#  GUI
# ---------------------------------------------------------------------------

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    _HAS_TK = True
except Exception:
    _HAS_TK = False


if _HAS_TK:

    class App(tk.Tk):
        BG = "#f4f6fa"
        ACCENT = "#0e7c66"
        TEXT = "#1f2937"
        MUTED = "#6b7280"

        def __init__(self):
            super().__init__()
            self.title(f"{APP_NAME} {VERSION} - Protein-Ligand Interaction Analysis")
            self.geometry("880x720")
            self.minsize(820, 660)
            self.configure(bg=self.BG)

            # --- icon (window + taskbar) ---
            _apply_icon(self)

            self.log_queue: "queue.Queue[Tuple[str, str]]" = queue.Queue()
            self.cancel_evt = threading.Event()
            self.worker: Optional[threading.Thread] = None
            self.pipeline: Optional[Pipeline] = None
            self.last_run_dir: Optional[Path] = None

            self._style()
            self._layout()
            self.after(120, self._drain_log)

        # ----------------------------------------------------------
        def _style(self):
            st = ttk.Style(self)
            st.theme_use("clam")
            st.configure(".", font=("Segoe UI", 10),
                         background=self.BG, foreground=self.TEXT)
            for w in ("TFrame", "TLabelframe", "TCheckbutton", "TRadiobutton", "TLabel"):
                st.configure(w, background=self.BG)
            st.configure("Title.TLabel", background=self.BG, foreground="#123a5c",
                         font=("Segoe UI", 15, "bold"))
            st.configure("Muted.TLabel", background=self.BG, foreground=self.MUTED)
            st.configure("TLabelframe.Label", background=self.BG,
                         foreground="#123a5c", font=("Segoe UI", 11, "bold"))
            st.configure("TButton", font=("Segoe UI", 10), padding=7)
            st.configure("Big.TButton", font=("Segoe UI", 13, "bold"),
                         background=self.ACCENT, foreground="white", padding=12)
            st.map("Big.TButton",
                   background=[("active", "#0a5c4c"), ("disabled", "#a9b8b3")],
                   foreground=[("disabled", "#eef2f0")])
            st.configure("Horizontal.TProgressbar", troughcolor="#e2e8f0",
                         background="#1f6feb", thickness=16)

        # ----------------------------------------------------------
        def _layout(self):
            pad = {"padx": 14, "pady": 6}

            ttk.Label(self, style="Title.TLabel",
                      text=f"{APP_NAME} - Protein-Ligand Interaction Analysis + PyMOL Session"
                      ).pack(anchor="w", padx=18, pady=(14, 2))
            ttk.Label(self, style="Muted.TLabel",
                      text="Pick files, click Analyze, and open the generated .pse in PyMOL."
                      ).pack(anchor="w", padx=18, pady=(0, 8))

            lf1 = ttk.LabelFrame(self, text="  Input  ")
            lf1.pack(fill="x", **pad)
            lf1.columnconfigure(1, weight=1)

            self.mode = tk.StringVar(value="two")
            ttk.Radiobutton(lf1, value="two", variable=self.mode,
                            text="Protein + ligand in two separate files",
                            command=self._mode_changed
                            ).grid(row=0, column=0, columnspan=3, sticky="w",
                                   padx=10, pady=(8, 2))
            ttk.Radiobutton(lf1, value="one", variable=self.mode,
                            text="One ready-made complex file (ligand is HETATM inside)",
                            command=self._mode_changed
                            ).grid(row=1, column=0, columnspan=3, sticky="w",
                                   padx=10, pady=(0, 6))

            self.lbl_p = ttk.Label(lf1, text="Protein:")
            self.ent_p = ttk.Entry(lf1)
            self.btn_p = ttk.Button(lf1, text="Browse...", width=12,
                                    command=lambda: self._pick(self.ent_p, PROTEIN_EXTS, "Protein"))
            self.lbl_l = ttk.Label(lf1, text="Ligand:")
            self.ent_l = ttk.Entry(lf1)
            self.btn_l = ttk.Button(lf1, text="Browse...", width=12,
                                    command=lambda: self._pick(self.ent_l, LIGAND_EXTS, "Ligand"))
            self.lbl_c = ttk.Label(lf1, text="Complex:")
            self.ent_c = ttk.Entry(lf1)
            self.btn_c = ttk.Button(lf1, text="Browse...", width=12,
                                    command=lambda: self._pick(self.ent_c, PROTEIN_EXTS, "Complex"))
            for e in (self.ent_p, self.ent_l, self.ent_c):
                e.bind("<FocusOut>", lambda _e: self._scan())

            self.scan_var = tk.StringVar(value="")
            ttk.Label(lf1, textvariable=self.scan_var, style="Muted.TLabel",
                      wraplength=820, justify="left"
                      ).grid(row=5, column=0, columnspan=3, sticky="w",
                             padx=12, pady=(2, 10))

            lf2 = ttk.LabelFrame(self, text="  Output folder  ")
            lf2.pack(fill="x", **pad)
            lf2.columnconfigure(1, weight=1)
            ttk.Label(lf2, text="Folder:").grid(row=0, column=0, sticky="e",
                                                padx=(10, 6), pady=6)
            self.ent_out = ttk.Entry(lf2)
            self.ent_out.grid(row=0, column=1, sticky="we", padx=4, pady=6)
            ttk.Button(lf2, text="Browse...", width=12, command=self._pick_out
                       ).grid(row=0, column=2, padx=(4, 10))
            ttk.Label(lf2, style="Muted.TLabel",
                      text="Each run creates its own sub-folder; nothing is overwritten."
                      ).grid(row=1, column=1, sticky="w", padx=4, pady=(0, 8))

            self.adv = tk.BooleanVar(value=False)
            self.btn_adv = ttk.Checkbutton(self, text="Advanced options",
                                           variable=self.adv, command=self._toggle_adv)
            self.btn_adv.pack(anchor="w", padx=20)
            self.lf_adv = ttk.LabelFrame(self, text="  Advanced  ")
            self.opt_rm_water = tk.BooleanVar(value=True)
            self.opt_add_h    = tk.BooleanVar(value=True)
            ttk.Checkbutton(self.lf_adv, variable=self.opt_rm_water,
                            text="Remove waters and common buffers from the complex"
                            ).grid(row=0, column=0, sticky="w", padx=16, pady=3)
            ttk.Checkbutton(self.lf_adv, variable=self.opt_add_h,
                            text="Add polar hydrogens to the ligand if missing"
                            ).grid(row=0, column=1, sticky="w", padx=16, pady=3)

            self.btn_run = ttk.Button(self, text="Analyze and save PyMOL session",
                                      style="Big.TButton", command=self._start)
            self.btn_run.pack(fill="x", padx=60, pady=(10, 4))

            bar = ttk.Frame(self)
            bar.pack(fill="x", padx=60, pady=(0, 4))
            self.pbar = ttk.Progressbar(bar, mode="determinate", maximum=100)
            self.pbar.pack(side="right", fill="x", expand=True)
            self.btn_cancel = ttk.Button(bar, text="Cancel", width=8,
                                         command=self._cancel, state="disabled")
            self.btn_cancel.pack(side="left", padx=(0, 8))

            self.status_var = tk.StringVar(value="Ready.")
            ttk.Label(self, textvariable=self.status_var,
                      wraplength=820, justify="left"
                      ).pack(anchor="w", padx=24, pady=(2, 6))

            lf3 = ttk.LabelFrame(self, text="  Results  ")
            lf3.pack(fill="both", expand=True, **pad)
            self.lst = tk.Listbox(lf3, font=("Consolas", 9),
                                  activestyle="none", height=6)
            self.lst.pack(fill="both", expand=True, padx=10, pady=(8, 4))
            bf = ttk.Frame(lf3)
            bf.pack(fill="x", padx=10, pady=(0, 10))
            ttk.Button(bf, text="Open .pse in PyMOL",
                       command=self._open_pymol).pack(side="left", padx=(0, 6))
            ttk.Button(bf, text="Open output folder",
                       command=self._open_folder).pack(side="left", padx=(0, 6))
            ttk.Button(bf, text="Open selected file",
                       command=self._open_selected).pack(side="left")

            ttk.Label(self, text="Log:", style="Muted.TLabel"
                      ).pack(anchor="w", padx=20, pady=(6, 0))
            self.txt = tk.Text(self, height=8, wrap="none",
                               bg="#101418", fg="#d7e2ec",
                               insertbackground="#d7e2ec",
                               font=("Consolas", 9), state="disabled")
            sb = ttk.Scrollbar(self, command=self.txt.yview)
            self.txt.configure(yscrollcommand=sb.set)
            sb.pack(side="right", fill="y", padx=(0, 14), pady=(0, 10))
            self.txt.pack(fill="both", expand=False, padx=(14, 0), pady=(0, 10))
            for tag, col in [("INFO", "#9fc3e8"), ("OK", "#7ee2a8"),
                             ("WARN", "#ffcc7a"), ("ERR", "#ff8f8f"),
                             ("CMD", "#8b98a5"), ("OUT", "#c6d3de")]:
                self.txt.tag_configure(tag, foreground=col)

            self._mode_changed()

        # ----------------------------------------------------------
        def _mode_changed(self):
            two = self.mode.get() == "two"
            for w in (self.lbl_p, self.ent_p, self.btn_p,
                      self.lbl_l, self.ent_l, self.btn_l,
                      self.lbl_c, self.ent_c, self.btn_c):
                w.grid_forget()
            if two:
                self.lbl_p.grid(row=2, column=0, sticky="e", padx=(10, 6), pady=4)
                self.ent_p.grid(row=2, column=1, sticky="we", padx=4, pady=4)
                self.btn_p.grid(row=2, column=2, padx=(4, 10))
                self.lbl_l.grid(row=3, column=0, sticky="e", padx=(10, 6), pady=4)
                self.ent_l.grid(row=3, column=1, sticky="we", padx=4, pady=4)
                self.btn_l.grid(row=3, column=2, padx=(4, 10))
            else:
                self.lbl_c.grid(row=2, column=0, sticky="e", padx=(10, 6), pady=4)
                self.ent_c.grid(row=2, column=1, sticky="we", padx=4, pady=4)
                self.btn_c.grid(row=2, column=2, padx=(4, 10))
            self._scan()

        def _toggle_adv(self):
            if self.adv.get():
                self.lf_adv.pack(fill="x", padx=14, pady=2, after=self.btn_adv)
            else:
                self.lf_adv.pack_forget()

        def _pick(self, entry, exts, title):
            ft = [(f"{title} files", " ".join("*" + e for e in exts)),
                  ("All files", "*.*")]
            p = filedialog.askopenfilename(title=f"Select {title} file", filetypes=ft)
            if p:
                entry.delete(0, "end")
                entry.insert(0, p)
                self._scan()
                if not self.ent_out.get().strip():
                    self.ent_out.insert(0, str(Path(p).parent / "PLIP_Results"))

        def _pick_out(self):
            p = filedialog.askdirectory(title="Select output folder")
            if p:
                self.ent_out.delete(0, "end")
                self.ent_out.insert(0, p)

        def _scan(self):
            try:
                p = (self.ent_c.get().strip() if self.mode.get() == "one"
                     else self.ent_l.get().strip() or self.ent_p.get().strip())
                if not p:
                    self.scan_var.set("")
                    return
                path = Path(p)
                if not path.exists():
                    self.scan_var.set("File does not exist.")
                    return
                if path.suffix.lower() not in LIGAND_EXTS + PROTEIN_EXTS:
                    self.scan_var.set("")
                    return
                groups = [g for g in scan_het_groups(path)
                          if g.code.upper() not in WATER_CODES]
                if groups:
                    txt = "Ligand detected: " + ", ".join(
                        f"{g.label} ({g.nheavy} heavy atoms)" for g in groups[:4])
                    if len(groups) > 4:
                        txt += f" and {len(groups)-4} more"
                    self.scan_var.set(txt)
                else:
                    atom_ligs = detect_atom_ligands(path)
                    if atom_ligs:
                        self.scan_var.set(
                            "Ligand stored as ATOM (" + ", ".join(atom_ligs) +
                            ") - will be rewritten as HETATM automatically.")
                    else:
                        self.scan_var.set("No ligand HETATM group found in this file.")
            except Exception:
                self.scan_var.set("")

        # ----------------------------------------------------------
        def _log(self, msg, level="INFO"):
            self.log_queue.put((msg, level))

        def _drain_log(self):
            try:
                while True:
                    msg, level = self.log_queue.get_nowait()
                    self.txt.configure(state="normal")
                    self.txt.insert("end", f"[{clock()}] {msg}\n", level)
                    self.txt.see("end")
                    self.txt.configure(state="disabled")
            except queue.Empty:
                pass
            self.after(120, self._drain_log)

        def _set_progress(self, frac): self.pbar["value"] = int(frac * 100)
        def _set_status(self, text):   self.status_var.set(text)

        # ----------------------------------------------------------
        def _gather(self) -> Options:
            def pn(s):
                s = s.strip()
                return Path(s) if s else None
            o = Options(
                protein=pn(self.ent_p.get()),
                ligand =pn(self.ent_l.get()),
                complex=pn(self.ent_c.get()),
                out_base=pn(self.ent_out.get()),
                remove_waters=self.opt_rm_water.get(),
                add_hydrogens=self.opt_add_h.get(),
            )
            if self.mode.get() == "two":
                o.complex = None
            else:
                o.protein = o.ligand = None
            return o

        def _start(self):
            if self.worker and self.worker.is_alive():
                return
            opts = self._gather()
            if opts.out_base is None:
                messagebox.showwarning("Output folder",
                                       "Please choose an output folder.")
                return
            if opts.complex is None and (opts.protein is None or opts.ligand is None):
                messagebox.showwarning("Missing input",
                                       "In two-file mode both protein and ligand are required.")
                return
            missing = [str(p) for p in (opts.protein, opts.ligand, opts.complex)
                       if p and not p.exists()]
            if missing:
                messagebox.showerror("File not found",
                                     "These files do not exist:\n" + "\n".join(missing))
                return

            self.lst.delete(0, "end")
            self.txt.configure(state="normal")
            self.txt.delete("1.0", "end")
            self.txt.configure(state="disabled")
            self.cancel_evt.clear()
            self.btn_run.configure(state="disabled")
            self.btn_cancel.configure(state="normal")
            self._set_progress(0)

            def work():
                t0 = time.time()
                runlog = RunLog(self._log)
                try:
                    pipe = Pipeline(opts, runlog,
                                    lambda s: self.after(0, self._set_status, s),
                                    lambda f: self.after(0, self._set_progress, f),
                                    self.cancel_evt)
                    self.pipeline = pipe
                    run_dir = pipe.run()
                    self.last_run_dir = run_dir
                    self.after(0, self._load_results, run_dir)
                    self.after(0, self._set_status,
                               f"Done in {time.time()-t0:.0f}s  ->  {run_dir}")
                    self.after(0, self._open_folder)
                except InterruptedError:
                    runlog("Cancelled.", "WARN")
                    self.after(0, self._set_status, "Cancelled.")
                except Exception as e:
                    tb = traceback.format_exc()
                    msg = str(e)
                    runlog(f"Error: {msg}", "ERR")
                    runlog(tb, "ERR")
                    logpath = runlog.write_error(tb)

                    def show(m=msg, lp=logpath):
                        self._set_status(f"Error: {m} | log: {lp}")
                        messagebox.showerror(
                            "Analysis failed",
                            f"{m}\n\nFull traceback saved to:\n{lp}")
                        self._os_open(lp)

                    self.after(0, show)
                finally:
                    runlog.close()
                    self.after(0, self._run_finished)

            self.worker = threading.Thread(target=work, daemon=True)
            self.worker.start()

        def _cancel(self):
            self.cancel_evt.set()
            proc = getattr(self.pipeline, "_proc", None) if self.pipeline else None
            if proc is not None and proc.poll() is None:
                try: proc.kill()
                except Exception: pass
            self._set_status("Cancelling ...")

        def _run_finished(self):
            self.btn_run.configure(state="normal")
            self.btn_cancel.configure(state="disabled")

        # ----------------------------------------------------------
        def _load_results(self, run_dir: Path):
            self.lst.delete(0, "end")
            for p in sorted(run_dir.iterdir()):
                if p.is_file():
                    self.lst.insert("end", f"{p.name}   ({p.stat().st_size:,} B)")

        def _selected(self) -> Optional[Path]:
            sel = self.lst.curselection()
            if not sel or not self.last_run_dir:
                return None
            return self.last_run_dir / self.lst.get(sel[0]).split("   (")[0]

        def _first_pse(self) -> Optional[Path]:
            if not self.last_run_dir:
                return None
            pses = sorted(self.last_run_dir.glob("*.pse"))
            return pses[0] if pses else None

        def _open_folder(self):
            if self.last_run_dir:
                self._os_open(self.last_run_dir)

        def _open_selected(self):
            p = self._selected()
            if p: self._os_open(p)

        def _open_pymol(self):
            p = self._selected()
            if p is None or p.suffix.lower() != ".pse":
                p = self._first_pse()
            if p is None:
                messagebox.showinfo("PyMOL", "No .pse found in the last run.")
                return
            self._set_status(f"Opening {p.name} in PyMOL ...")
            try:
                flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
                subprocess.Popen([_python_exe(), "-m", "pymol", str(p)],
                                 creationflags=flags)
            except Exception as e:
                messagebox.showerror("PyMOL", str(e))

        @staticmethod
        def _os_open(path: Path):
            try:
                if sys.platform == "win32":
                    os.startfile(str(path))  # type: ignore
                elif sys.platform == "darwin":
                    subprocess.Popen(["open", str(path)])
                else:
                    subprocess.Popen(["xdg-open", str(path)])
            except Exception as e:
                messagebox.showerror("Open", str(e))


# ---------------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------------

def run_headless(args) -> int:
    opts = Options(
        protein=Path(args.protein) if args.protein else None,
        ligand =Path(args.ligand)  if args.ligand  else None,
        complex=Path(args.complex) if args.complex else None,
        out_base=Path(args.out),
        runname=args.name or "",
        remove_waters=not args.keep_waters,
        add_hydrogens=not args.no_h,
        keepmod=args.keepmod,
        dnareceptor=args.dnareceptor,
        threads=args.threads,
    )

    def log(msg, level):
        print(f"[{level}] {msg}")

    runlog = RunLog(log)
    pipe = Pipeline(opts, runlog, lambda s: None, lambda f: None, threading.Event())
    try:
        pipe.run()
        return 0
    except Exception as e:
        tb = traceback.format_exc()
        log(f"Error: {e}", "ERR")
        lp = runlog.write_error(tb)
        log(f"Full traceback saved to: {lp}", "ERR")
        print(tb)
        return 1
    finally:
        runlog.close()


def _deps_missing() -> List[Tuple[str, str]]:
    mods = [("plip", "plip"), ("openbabel", "openbabel"),
            ("lxml", "lxml"), ("pymol", "pymol-open-source")]
    miss = []
    for mod, pip in mods:
        try:
            __import__(mod)
        except Exception:
            miss.append((mod, pip))
    return miss


def _run_quiet(cmd: List[str], timeout: int = 60):
    flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
    return subprocess.run(cmd, capture_output=True, text=True,
                          creationflags=flags, timeout=timeout)


def _interp_has_deps(py: str) -> bool:
    try:
        return _run_quiet([py, "-c",
                           "import plip, openbabel, lxml, pymol"]).returncode == 0
    except Exception:
        return False


def _candidate_interpreters() -> List[str]:
    cands: List[str] = []
    try:
        out = _run_quiet(["py", "-0p"], timeout=15).stdout or ""
        for line in out.splitlines():
            m = re.match(r"\s*-V:\S+\s*\*?\s*(.+)$", line)
            if m:
                cands.append(m.group(1).strip())
    except Exception:
        pass
    for exe in ("python", "pythonw"):
        try:
            out = _run_quiet(["where", exe], timeout=15).stdout or ""
            cands += [ln.strip() for ln in out.splitlines() if ln.strip()]
        except Exception:
            pass
    for base in (Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "Python",
                 Path(os.environ.get("PROGRAMFILES", "")) / "Python"):
        try:
            if base.is_dir():
                for d in sorted(base.iterdir()):
                    for name in ("python.exe", "pythonw.exe"):
                        q = d / name
                        if q.exists():
                            cands.append(str(q))
        except Exception:
            pass
    seen, out = set(), []
    for c in cands:
        k = os.path.normcase(os.path.abspath(c))
        if k not in seen:
            seen.add(k)
            out.append(c)
    return out


def _find_ready_interpreter() -> Optional[str]:
    me = os.path.normcase(os.path.abspath(sys.executable))
    for cand in _candidate_interpreters():
        p = os.path.normcase(os.path.abspath(cand))
        if p == me:
            continue
        tester = cand
        if tester.lower().endswith("pythonw.exe"):
            twin = tester[:-len("w.exe")] + ".exe"
            tester = twin if os.path.exists(twin) else cand
        if _interp_has_deps(tester):
            if tester.lower().endswith("python.exe"):
                w = tester[:-len(".exe")] + "w.exe"
                if os.path.exists(w):
                    return w
            return tester
    return None


def _pip_install_gui(missing) -> bool:
    # make sure the taskbar icon is used before Tk creates any window
    _set_windows_appid()

    pkgs = [pip for _, pip in missing]
    root = tk.Tk(); root.withdraw()
    _apply_icon(root)                       # icon for the hidden root

    yes = messagebox.askyesno(
        "Install dependencies",
        "These packages are missing in the current Python:\n\n"
        + ", ".join(m for m, _ in missing)
        + "\n\nInstall them now with pip? (requires internet)")
    if not yes:
        root.destroy()
        return False

    win = tk.Toplevel(root)
    win.title("Installing packages ...")
    win.geometry("780x440")
    _apply_icon(win)                        # icon for the install window

    txt = tk.Text(win, wrap="word", font=("Consolas", 9),
                  bg="#101418", fg="#d7e2ec")
    txt.pack(fill="both", expand=True, padx=8, pady=8)
    lbl = tk.Label(win, text="Please wait - do not close this window.",
                   bg="#f4f6fa", fg="#1f2937",
                   font=("Segoe UI", 10, "bold"))
    lbl.pack(side="bottom", fill="x", padx=8, pady=(0, 8))
    win.update()

    flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
    proc = subprocess.Popen([_python_exe(), "-m", "pip", "install", *pkgs],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace",
                            creationflags=flags)
    assert proc.stdout is not None
    for line in proc.stdout:
        txt.insert("end", line); txt.see("end"); win.update()
    rc = proc.wait()
    lbl.configure(text="Install finished." if rc == 0 else "Install failed.")
    win.update()
    messagebox.showinfo("pip", "Install finished." if rc == 0 else "Install failed.")
    root.destroy()
    return rc == 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="plip_studio", description=APP_NAME)
    ap.add_argument("--headless", action="store_true")
    ap.add_argument("--protein"); ap.add_argument("--ligand"); ap.add_argument("--complex")
    ap.add_argument("--out", default="./PLIP_Results")
    ap.add_argument("--name", default="")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--keepmod", action="store_true")
    ap.add_argument("--dnareceptor", action="store_true")
    ap.add_argument("--no-h", action="store_true")
    ap.add_argument("--keep-waters", action="store_true")
    args = ap.parse_args(argv)

    # bootstrap: re-exec in another interpreter that has the packages
    if os.environ.get("PLIPSTUDIO_REEXEC") != "1" and _deps_missing():
        good = _find_ready_interpreter()
        if good:
            print(f"Dependencies found in another interpreter: {good}")
            env = dict(os.environ)
            env["PLIPSTUDIO_REEXEC"] = "1"
            flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
            subprocess.Popen([good, os.path.abspath(__file__), *sys.argv[1:]],
                             env=env, creationflags=flags)
            return 0

    missing = _deps_missing()
    if args.headless:
        if missing:
            print("Missing packages: " + ", ".join(m for m, _ in missing))
            print("Install with:  python -m pip install " +
                  " ".join(p for _, p in missing))
            return 1
        return run_headless(args)

    if missing:
        if not _HAS_TK:
            print("Missing packages: " + ", ".join(m for m, _ in missing))
            return 1
        if not _pip_install_gui(missing):
            return 1

    if not _HAS_TK:
        print("Tkinter is not available in this Python; use --headless.")
        return 1

    # set the taskbar AppUserModelID *before* the first Tk window exists
    _set_windows_appid()

    app = App()
    app.mainloop()
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        crash = Path(__file__).resolve().parent / CRASH_LOG
        tb = traceback.format_exc()
        try:
            with open(crash, "a", encoding="utf-8") as fh:
                fh.write(f"\n===== {human_time()} =====\n{tb}")
        except Exception:
            pass
        print(tb)
        if _HAS_TK:
            try:
                _set_windows_appid()
                r = tk.Tk(); r.withdraw()
                _apply_icon(r)
                messagebox.showerror("Unexpected error",
                                     f"See: {crash}")
                r.destroy()
            except Exception:
                pass
        try:
            if sys.platform == "win32":
                os.startfile(str(crash))  # type: ignore
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(crash)])
            else:
                subprocess.Popen(["xdg-open", str(crash)])
        except Exception:
            pass
        sys.exit(1)