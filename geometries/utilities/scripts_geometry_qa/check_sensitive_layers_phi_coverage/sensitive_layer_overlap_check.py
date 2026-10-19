#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GDML sensitive-overlap scan (GDML units: cm; distance outputs: mm):
- Fixed eta
- Scan phi with configurable step
- Count intersections with "sensitive" volumes defined by MATERIAL NAMES

Robust PyROOT navigation (PyROOT-safe):
- FindNextBoundary + GetStep + manual point update + FindNode
- Avoids TGeoManager.Step(double) binding issues

ALL outputs are written into:
  ouput_overlap_analysis_for_<gdml_file_name_without_ext>/

ROOT (+ PNG) extras:
- Export from ROOT to PNG (with gPad->SetGrid()):
    * h_sensitive_freq (lin + logY)
    * h_large_hole_width_mm (lin + logY) [global]
    * h_large_hole_width_mm_L{layer} (lin + logY) [per-layer]

If --splitIntoRadialLayers is enabled:
Counts:
  * TH2D h_layer_phi                    : X=phi, Y=layer, Z=counts
  * TH1D h_layer_phi_L{layer}           : per-layer counts vs phi
  * TH2D h_layer_s                      : X=s=phi*R(mm), Y=layer, Z=counts (shared x bins)
  * TH1D h_layer_s_L{layer}             : per-layer counts vs s(mm)
Holes (large only; hole width in mm > maxAllowedHoleSize):
  * TH2D h_layer_holes_phi              : X=phi, Y=layer, Z=hole flag (0/1)
  * TH1D h_layer_holes_phi_L{layer}     : per-layer hole flag vs phi (red)
  * TH2D h_layer_holes_s                : X=s=phi*R(mm), Y=layer, Z=hole flag (0/1)
  * TH1D h_layer_holes_s_L{layer}       : per-layer hole flag vs s(mm) (red)
Hole-width distributions:
  * TH1D h_large_hole_width_mm          : distribution of hole widths [mm] (global, only > threshold)
  * TH1D h_large_hole_width_mm_L{layer} : per-layer distributions (only > threshold)

IMPORTANT:
- GDML coordinates are in cm.
- We convert radius to mm internally with CM_TO_MM=10.
- --maxAllowedHoleSize is in mm.

Requires: ROOT (PyROOT), numpy, matplotlib
"""

import argparse
import csv
import math
import os
import re
import sys
from collections import Counter
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

CM_TO_MM = 10.0  # GDML is cm; we report distances in mm


# -------------------------
# Geometry / direction utils
# -------------------------

def eta_to_theta(eta: float) -> float:
    return 2.0 * math.atan(math.exp(-eta))


def direction_from_eta_phi(eta: float, phi: float) -> Tuple[float, float, float]:
    theta = eta_to_theta(eta)
    st = math.sin(theta)
    ct = math.cos(theta)
    dx = st * math.cos(phi)
    dy = st * math.sin(phi)
    dz = ct
    n = math.sqrt(dx * dx + dy * dy + dz * dz)
    if n == 0:
        return 0.0, 0.0, 1.0
    return dx / n, dy / n, dz / n


def safe_mat_name(vol) -> str:
    try:
        m = vol.GetMaterial()
        if m:
            return str(m.GetName())
    except Exception:
        pass
    return ""


def safe_vol_name(vol) -> str:
    try:
        return str(vol.GetName())
    except Exception:
        pass
    return ""


def is_sensitive_volume(vol, sensitive_materials_set: Set[str]) -> bool:
    return safe_mat_name(vol) in sensitive_materials_set


# -------------------------
# ROOT / GDML loading
# -------------------------

def load_gdml_into_root(gdml_path: str):
    try:
        import ROOT  # noqa
    except Exception as e:
        print("ERROR: Cannot import ROOT (PyROOT). Make sure ROOT is set up in your env.", file=sys.stderr)
        raise e

    import ROOT
    geom = ROOT.TGeoManager.Import(gdml_path)
    if not geom:
        raise RuntimeError(f"Failed to import GDML: {gdml_path}")
    return ROOT.gGeoManager


# -------------------------
# Reporting helpers
# -------------------------

def print_count_frequencies(counts, decimals: int = 2, top: int = 0):
    counts_list = [int(x) for x in counts]
    total = len(counts_list)
    if total == 0:
        print("No phi points to summarize.")
        return

    freq = Counter(counts_list)
    items = sorted(freq.items(), key=lambda kv: (-kv[1], -kv[0]))
    if top and top > 0:
        items = items[:top]

    parts = []
    for n, c in items:
        pct = 100.0 * c / total
        parts.append(f"{n} - {pct:.{decimals}f}%")

    print("Frequency of N_sensitive over phi scan:")
    print(", ".join(parts))


def mode_of_counts(counts) -> int:
    counts_list = [int(x) for x in counts]
    if not counts_list:
        return 0
    freq = Counter(counts_list)
    return max(freq.items(), key=lambda kv: (kv[1], kv[0]))[0]


# -------------------------
# Layer parsing from path
# -------------------------
_LAYER_RE = re.compile(r"(?:^|/)(TRKLayer(\d+))(?:[^/]*)(?:$|/)")


def layer_from_path(path: str) -> Optional[int]:
    last = None
    for m in _LAYER_RE.finditer(path):
        last = m
    if not last:
        return None
    return int(last.group(2))


# -------------------------
# Ray tracing / scan per phi
# -------------------------

def trace_phi(
    gGeoManager,
    eta: float,
    phi: float,
    sensitive_materials: List[str],
    start_xyz: Tuple[float, float, float],
    eps: float,
    max_steps: int,
    max_track_len: Optional[float],
    unique: bool,
    split_into_layers: bool,
    debug: bool = False,
) -> Tuple[int, Optional[Dict[int, int]], Optional[Dict[int, float]]]:
    """
    Returns:
      total_count
      layer_counts (layer->count) if split_into_layers else None
      layer_r_est_mm (layer->R_mm) if split_into_layers else None
    """
    sens = set(sensitive_materials)

    x, y, z = map(float, start_xyz)  # cm
    dx, dy, dz = direction_from_eta_phi(eta, phi)

    gGeoManager.InitTrack(x, y, z, dx, dy, dz)
    if gGeoManager.IsOutside():
        return 0, ({} if split_into_layers else None), ({} if split_into_layers else None)

    prev_node = gGeoManager.GetCurrentNode()
    travelled = 0.0
    stuck_counter = 0

    total_entries = 0
    seen_unique: Set[str] = set()

    layer_counts: Dict[int, int] = {}
    layer_r_sum_mm: Dict[int, float] = {}
    layer_r_n: Dict[int, int] = {}

    for istep in range(max_steps):
        if max_track_len is not None and travelled >= max_track_len:
            break

        try:
            _ = gGeoManager.FindNextBoundary()
        except Exception as e:
            if debug:
                print(f"[debug] FindNextBoundary exception at phi={phi:.6f}, step={istep}: {e}", file=sys.stderr)
            break

        step = float(gGeoManager.GetStep())
        if (not math.isfinite(step)) or step <= 0.0:
            step = eps
            stuck_counter += 1
        else:
            stuck_counter = 0

        if stuck_counter > 1000:
            if debug:
                print(f"[debug] Stuck >1000 at phi={phi:.6f}; breaking", file=sys.stderr)
            break

        take = step + eps
        x += take * dx
        y += take * dy
        z += take * dz
        travelled += take

        gGeoManager.SetCurrentPoint(x, y, z)
        node = gGeoManager.FindNode()
        if not node or gGeoManager.IsOutside():
            break

        if node != prev_node:
            vol = node.GetVolume()
            if vol and is_sensitive_volume(vol, sens):
                if split_into_layers:
                    try:
                        path = str(gGeoManager.GetPath())
                    except Exception:
                        path = ""
                    lay = layer_from_path(path)
                    if lay is not None:
                        layer_counts[lay] = layer_counts.get(lay, 0) + 1
                        r_here_cm = math.sqrt(x * x + y * y)
                        r_here_mm = r_here_cm * CM_TO_MM
                        layer_r_sum_mm[lay] = layer_r_sum_mm.get(lay, 0.0) + r_here_mm
                        layer_r_n[lay] = layer_r_n.get(lay, 0) + 1

                if unique:
                    key = f"{int(node):x}|{safe_vol_name(vol)}|{safe_mat_name(vol)}"
                    if key not in seen_unique:
                        seen_unique.add(key)
                        total_entries += 1
                else:
                    total_entries += 1

            prev_node = node

    if not split_into_layers:
        return total_entries, None, None

    layer_r_est_mm: Dict[int, float] = {}
    for lay, s in layer_r_sum_mm.items():
        n = layer_r_n.get(lay, 0)
        if n > 0:
            layer_r_est_mm[lay] = s / n

    return total_entries, layer_counts, layer_r_est_mm


# -------------------------
# Holes (gap) logic per layer
# -------------------------

def compute_large_holes_mask_phi(
    counts_per_phi: np.ndarray,
    phi_step: float,
    R_mm: float,
    max_allowed_hole_mm: float,
) -> np.ndarray:
    """
    Hole = contiguous region where counts==0 and width_mm > threshold.
    width_mm = n_bins * phi_step * R_mm
    Returns 0/1 mask per phi bin.
    """
    n = len(counts_per_phi)
    mask = np.zeros(n, dtype=int)
    if n == 0 or R_mm <= 0:
        return mask

    zeros = (counts_per_phi == 0)

    i = 0
    while i < n:
        if not zeros[i]:
            i += 1
            continue
        j = i
        while j < n and zeros[j]:
            j += 1
        n_bins = j - i
        width_mm = (n_bins * phi_step) * R_mm
        if width_mm > max_allowed_hole_mm:
            mask[i:j] = 1
        i = j

    return mask


def extract_large_hole_widths_mm(
    counts_per_phi: np.ndarray,
    phi_step: float,
    R_mm: float,
    max_allowed_hole_mm: float,
) -> List[float]:
    """
    Return list of widths [mm] for large holes only (width > threshold).
    """
    widths: List[float] = []
    n = len(counts_per_phi)
    if n == 0 or R_mm <= 0:
        return widths

    zeros = (counts_per_phi == 0)
    i = 0
    while i < n:
        if not zeros[i]:
            i += 1
            continue
        j = i
        while j < n and zeros[j]:
            j += 1
        n_bins = j - i
        width_mm = (n_bins * phi_step) * R_mm
        if width_mm > max_allowed_hole_mm:
            widths.append(float(width_mm))
        i = j

    return widths


# -------------------------
# ROOT: export histograms to PNG (grid + lin/logY)
# -------------------------

def root_export_hist_png(root_path: str, hist_name: str, png_path: str, logy: bool):
    import ROOT

    ROOT.gROOT.SetBatch(True)

    f = ROOT.TFile.Open(root_path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {root_path}")

    h = f.Get(hist_name)
    if not h:
        f.Close()
        raise RuntimeError(f"Histogram '{hist_name}' not found in {root_path}")

    c = ROOT.TCanvas(f"c_{hist_name}_{'log' if logy else 'lin'}", "", 900, 650)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetLogy(1 if logy else 0)

    h.SetLineWidth(2)
    h.Draw("HIST")

    c.SaveAs(png_path)
    c.Close()
    f.Close()


# -------------------------
# ROOT outputs
# -------------------------

def save_root_histograms_basic(root_filename, phi_vals, counts, unique: bool):
    import ROOT

    counts_int = [int(x) for x in counts]
    if len(counts_int) == 0:
        raise RuntimeError("Empty counts array; nothing to save.")

    fout = ROOT.TFile(root_filename, "RECREATE")

    nbin = len(phi_vals)
    phi_min = float(phi_vals[0])
    phi_max = float(phi_vals[-1])

    htitle_phi = "Sensitive vs phi;phi [rad];N sensitive"
    if unique:
        htitle_phi = "Unique sensitive vs phi;phi [rad];N unique sensitive"

    h_phi = ROOT.TH1D("h_sensitive_phi", htitle_phi, nbin, phi_min, phi_max)
    for i, c in enumerate(counts_int, start=1):
        h_phi.SetBinContent(i, c)
    h_phi.Write()

    freq = Counter(counts_int)
    nmin = min(freq.keys())
    nmax = max(freq.keys())

    htitle_freq = "Frequency of N_{sensitive};N sensitive;Counts"
    if unique:
        htitle_freq = "Frequency of N_{unique sensitive};N unique sensitive;Counts"

    h_freq = ROOT.TH1D("h_sensitive_freq", htitle_freq, (nmax - nmin + 1), nmin - 0.5, nmax + 0.5)
    for n, cnt in freq.items():
        b = h_freq.FindBin(float(n))
        h_freq.SetBinContent(b, cnt)
        h_freq.SetBinError(b, math.sqrt(cnt))
    h_freq.Write()

    m = mode_of_counts(counts_int)
    h_anom = ROOT.TH1D(
        "h_sensitive_anomaly_phi",
        f"Anomaly vs phi (N - mode={m});phi [rad];N - mode",
        nbin, phi_min, phi_max
    )
    for i, c in enumerate(counts_int, start=1):
        h_anom.SetBinContent(i, c - m)
    h_anom.Write()

    try:
        p_mode = ROOT.TParameter(int)("mode_N_sensitive", int(m))
        p_mode.Write()
    except Exception:
        pass

    fout.Close()
    print(f"Wrote ROOT histograms (basic): {root_filename}")


def _write_hole_width_hist(
    ROOT,
    fout,
    hname: str,
    title: str,
    widths_mm: List[float],
    threshold_mm: float,
):
    """
    Writes a TH1D with hole widths. If widths list is empty, writes an empty hist with a sensible range.
    """
    if widths_mm:
        wmax = max(widths_mm)
        xmax = max(wmax * 1.05, threshold_mm * 1.05)
        nb = 100
        h = ROOT.TH1D(hname, title, nb, 0.0, xmax)
        for w in widths_mm:
            h.Fill(float(w))
        h.Write()
    else:
        h = ROOT.TH1D(hname, title, 50, 0.0, threshold_mm * 5.0)
        h.Write()


def save_root_histograms_layers(
    root_filename: str,
    phi_vals: np.ndarray,
    layer_ids_sorted: List[int],
    layer_counts_by_phi: List[Dict[int, int]],
    layer_radii_mm: Dict[int, float],
    phi_step: float,
    max_allowed_hole_mm: float,
    unique: bool,
) -> List[int]:
    """
    Appends layer histograms to existing ROOT file.
    Returns the list of layers for which per-layer hole-width hists were written (same as layer_ids_sorted).
    """
    import ROOT

    fout = ROOT.TFile(root_filename, "UPDATE")
    if not fout or fout.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file for update: {root_filename}")

    nphi = len(phi_vals)
    phi_min = float(phi_vals[0])
    phi_max = float(phi_vals[-1])

    lay_min = min(layer_ids_sorted)
    lay_max = max(layer_ids_sorted)
    nlay_bins = lay_max - lay_min + 1

    lay_to_row = {lay: idx for idx, lay in enumerate(layer_ids_sorted)}
    mat_phi_counts = np.zeros((len(layer_ids_sorted), nphi), dtype=int)
    for j in range(nphi):
        for lay, cnt in layer_counts_by_phi[j].items():
            if lay in lay_to_row:
                mat_phi_counts[lay_to_row[lay], j] += int(cnt)

    hole_mask_by_layer: Dict[int, np.ndarray] = {}
    hole_widths_by_layer: Dict[int, List[float]] = {}
    large_hole_widths_mm_global: List[float] = []

    for lay in layer_ids_sorted:
        R = layer_radii_mm.get(lay, None)
        row = lay_to_row[lay]
        if R is None or R <= 0:
            hole_mask_by_layer[lay] = np.zeros(nphi, dtype=int)
            hole_widths_by_layer[lay] = []
            continue

        hole_mask_by_layer[lay] = compute_large_holes_mask_phi(
            mat_phi_counts[row, :],
            phi_step=phi_step,
            R_mm=float(R),
            max_allowed_hole_mm=float(max_allowed_hole_mm),
        )
        widths_lay = extract_large_hole_widths_mm(
            mat_phi_counts[row, :],
            phi_step=phi_step,
            R_mm=float(R),
            max_allowed_hole_mm=float(max_allowed_hole_mm),
        )
        hole_widths_by_layer[lay] = widths_lay
        large_hole_widths_mm_global.extend(widths_lay)

    # ---- 2D counts vs phi and layer
    title2d = "Sensitive counts vs phi and layer;phi [rad];Layer index;Counts"
    if unique:
        title2d = "Unique-sensitive counts vs phi and layer;phi [rad];Layer index;Counts"
    h2_phi = ROOT.TH2D("h_layer_phi", title2d, nphi, phi_min, phi_max, nlay_bins, lay_min - 0.5, lay_max + 0.5)

    # ---- 2D holes vs phi and layer
    h2_holes_phi = ROOT.TH2D(
        "h_layer_holes_phi",
        f"Holes (large only) vs phi and layer;phi [rad];Layer index;Hole flag (1=large hole)",
        nphi, phi_min, phi_max, nlay_bins, lay_min - 0.5, lay_max + 0.5
    )

    # ---- per-layer 1D counts vs phi & holes vs phi
    h1_phi_by_layer: Dict[int, "ROOT.TH1D"] = {}
    h1_holes_phi_by_layer: Dict[int, "ROOT.TH1D"] = {}

    for lay in layer_ids_sorted:
        h1_phi_by_layer[lay] = ROOT.TH1D(
            f"h_layer_phi_L{lay}",
            (f"Counts vs phi for TRKLayer{lay};phi [rad];Counts"
             if not unique else f"Unique counts vs phi for TRKLayer{lay};phi [rad];Counts"),
            nphi, phi_min, phi_max
        )
        h1_holes_phi_by_layer[lay] = ROOT.TH1D(
            f"h_layer_holes_phi_L{lay}",
            f"Holes (large only) vs phi for TRKLayer{lay};phi [rad];Hole flag (1=large hole)",
            nphi, phi_min, phi_max
        )
        h1_holes_phi_by_layer[lay].SetLineColor(ROOT.kRed)

    for j, phi in enumerate(phi_vals):
        for lay in layer_ids_sorted:
            c = mat_phi_counts[lay_to_row[lay], j]
            if c:
                h2_phi.Fill(float(phi), float(lay), float(c))
                h1_phi_by_layer[lay].SetBinContent(j + 1, float(c))
            hm = int(hole_mask_by_layer[lay][j])
            if hm:
                h2_holes_phi.Fill(float(phi), float(lay), 1.0)
                h1_holes_phi_by_layer[lay].SetBinContent(j + 1, 1.0)

    h2_phi.Write()
    h2_holes_phi.Write()
    for lay in layer_ids_sorted:
        h1_phi_by_layer[lay].Write()
        h1_holes_phi_by_layer[lay].Write()

    # ---- Distance histograms: s=phi*R_mm
    maxR = max(layer_radii_mm.values()) if layer_radii_mm else 0.0
    if maxR > 0:
        s_min = phi_min * maxR
        s_max = phi_max * maxR

        title2d_s = "Sensitive counts vs s=phi*R and layer;s [mm];Layer index;Counts"
        if unique:
            title2d_s = "Unique-sensitive counts vs s=phi*R and layer;s [mm];Layer index;Counts"
        h2_s = ROOT.TH2D("h_layer_s", title2d_s, nphi, s_min, s_max, nlay_bins, lay_min - 0.5, lay_max + 0.5)

        h2_holes_s = ROOT.TH2D(
            "h_layer_holes_s",
            f"Holes (large only) vs s=phi*R and layer;s [mm];Layer index;Hole flag (1=large hole)",
            nphi, s_min, s_max, nlay_bins, lay_min - 0.5, lay_max + 0.5
        )

        h1_s_by_layer: Dict[int, "ROOT.TH1D"] = {}
        h1_holes_s_by_layer: Dict[int, "ROOT.TH1D"] = {}

        for lay in layer_ids_sorted:
            R = layer_radii_mm.get(lay, None)
            if R is None or R <= 0:
                continue
            smin_l = phi_min * R
            smax_l = phi_max * R

            h1_s_by_layer[lay] = ROOT.TH1D(
                f"h_layer_s_L{lay}",
                (f"Counts vs s=phi*R for TRKLayer{lay};s [mm];Counts"
                 if not unique else f"Unique counts vs s=phi*R for TRKLayer{lay};s [mm];Counts"),
                nphi, smin_l, smax_l
            )
            h1_holes_s_by_layer[lay] = ROOT.TH1D(
                f"h_layer_holes_s_L{lay}",
                f"Holes (large only) vs s=phi*R for TRKLayer{lay};s [mm];Hole flag (1=large hole)",
                nphi, smin_l, smax_l
            )
            h1_holes_s_by_layer[lay].SetLineColor(ROOT.kRed)

        for j, phi in enumerate(phi_vals):
            for lay in layer_ids_sorted:
                R = layer_radii_mm.get(lay, None)
                if R is None or R <= 0:
                    continue
                c = mat_phi_counts[lay_to_row[lay], j]
                s = float(phi) * float(R)
                if c:
                    h2_s.Fill(s, float(lay), float(c))
                    if lay in h1_s_by_layer:
                        h1_s_by_layer[lay].SetBinContent(j + 1, float(c))
                hm = int(hole_mask_by_layer[lay][j])
                if hm:
                    h2_holes_s.Fill(s, float(lay), 1.0)
                    if lay in h1_holes_s_by_layer:
                        h1_holes_s_by_layer[lay].SetBinContent(j + 1, 1.0)

        h2_s.Write()
        h2_holes_s.Write()
        for lay in layer_ids_sorted:
            if lay in h1_s_by_layer:
                h1_s_by_layer[lay].Write()
            if lay in h1_holes_s_by_layer:
                h1_holes_s_by_layer[lay].Write()

    # ---- Hole width distributions (global + per-layer)  (NEW per-layer)
    _write_hole_width_hist(
        ROOT, fout,
        "h_large_hole_width_mm",
        f"Large hole width distribution (global);Hole width [mm];Counts (holes)   (threshold>{max_allowed_hole_mm} mm)",
        large_hole_widths_mm_global,
        max_allowed_hole_mm,
    )

    for lay in layer_ids_sorted:
        _write_hole_width_hist(
            ROOT, fout,
            f"h_large_hole_width_mm_L{lay}",
            f"Large hole width distribution (TRKLayer{lay});Hole width [mm];Counts (holes)   (threshold>{max_allowed_hole_mm} mm)",
            hole_widths_by_layer.get(lay, []),
            max_allowed_hole_mm,
        )

    # Store radii + threshold
    for lay, R in layer_radii_mm.items():
        try:
            p = ROOT.TParameter(float)(f"TRKLayer{lay}_R_mm", float(R))
            p.Write()
        except Exception:
            pass
    try:
        p_hole = ROOT.TParameter(float)("maxAllowedHoleSize_mm", float(max_allowed_hole_mm))
        p_hole.Write()
    except Exception:
        pass

    fout.Close()
    print(f"Appended ROOT histograms (layers + holes + hole-width dists) into: {root_filename}")

    return layer_ids_sorted


# -------------------------
# PNG plotting (matplotlib) for layer maps / per-layer
# -------------------------

def _ensure_matplotlib_backend(no_show: bool):
    if no_show:
        import matplotlib
        matplotlib.use("Agg")


def save_pngs_layers(
    out_prefix: str,
    phi_vals: np.ndarray,
    layer_ids_sorted: List[int],
    layer_counts_by_phi: List[Dict[int, int]],
    layer_radii_mm: Dict[int, float],
    phi_step: float,
    max_allowed_hole_mm: float,
    unique: bool,
):
    import matplotlib.pyplot as plt

    lay_to_row = {lay: i for i, lay in enumerate(layer_ids_sorted)}
    nlay = len(layer_ids_sorted)
    nphi = len(phi_vals)

    mat_phi = np.zeros((nlay, nphi), dtype=float)
    for j in range(nphi):
        d = layer_counts_by_phi[j]
        for lay, cnt in d.items():
            if lay in lay_to_row:
                mat_phi[lay_to_row[lay], j] += float(cnt)

    mat_holes_phi = np.zeros((nlay, nphi), dtype=float)
    for lay in layer_ids_sorted:
        R = layer_radii_mm.get(lay, None)
        if R is None or R <= 0:
            continue
        row = lay_to_row[lay]
        mask = compute_large_holes_mask_phi(
            counts_per_phi=mat_phi[row, :].astype(int),
            phi_step=phi_step,
            R_mm=float(R),
            max_allowed_hole_mm=float(max_allowed_hole_mm),
        )
        mat_holes_phi[row, :] = mask.astype(float)

    # counts vs phi (2D)
    plt.figure()
    plt.imshow(
        mat_phi, aspect="auto", origin="lower",
        extent=[float(phi_vals[0]), float(phi_vals[-1]),
                float(layer_ids_sorted[0]) - 0.5, float(layer_ids_sorted[-1]) + 0.5],
    )
    plt.xlabel("phi [rad]")
    plt.ylabel("Layer index")
    plt.title("Counts vs phi and layer" + (" (unique)" if unique else ""))
    plt.colorbar(label="Counts")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_layer_phi_2d.png", dpi=200)
    plt.close()

    # holes vs phi (2D)
    plt.figure()
    plt.imshow(
        mat_holes_phi, aspect="auto", origin="lower",
        extent=[float(phi_vals[0]), float(phi_vals[-1]),
                float(layer_ids_sorted[0]) - 0.5, float(layer_ids_sorted[-1]) + 0.5],
        vmin=0.0, vmax=1.0
    )
    plt.xlabel("phi [rad]")
    plt.ylabel("Layer index")
    plt.title(f"Holes (large only) vs phi and layer (>{max_allowed_hole_mm} mm)")
    plt.colorbar(label="Hole flag (0/1)")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_layer_holes_phi_2d.png", dpi=200)
    plt.close()

    # per-layer plots + holes (red)
    for lay in layer_ids_sorted:
        row = lay_to_row[lay]
        y = mat_phi[row, :]

        plt.figure()
        plt.plot(phi_vals, y)
        plt.xlabel("phi [rad]")
        plt.ylabel("Counts")
        plt.title(f"TRKLayer{lay}: counts vs phi" + (" (unique)" if unique else ""))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_L{lay}_phi.png", dpi=200)
        plt.close()

        plt.figure()
        plt.plot(phi_vals, mat_holes_phi[row, :], color="red")
        plt.ylim(-0.1, 1.1)
        plt.xlabel("phi [rad]")
        plt.ylabel("Hole flag (1=large hole)")
        plt.title(f"TRKLayer{lay}: large holes vs phi (>{max_allowed_hole_mm} mm)")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_L{lay}_holes_phi.png", dpi=200)
        plt.close()

    if not layer_radii_mm:
        return

    maxR = max(layer_radii_mm.values())
    s_extent_min = float(phi_vals[0]) * float(maxR)
    s_extent_max = float(phi_vals[-1]) * float(maxR)
    denom = (s_extent_max - s_extent_min) if (s_extent_max > s_extent_min) else 1.0

    mat_s = np.zeros((nlay, nphi), dtype=float)
    mat_holes_s = np.zeros((nlay, nphi), dtype=float)

    for j, phi in enumerate(phi_vals):
        for lay in layer_ids_sorted:
            R = layer_radii_mm.get(lay, None)
            if R is None or R <= 0:
                continue
            row = lay_to_row[lay]
            s = float(phi) * float(R)
            col = int(round((s - s_extent_min) / denom * (nphi - 1)))
            if 0 <= col < nphi:
                if mat_phi[row, j] != 0.0:
                    mat_s[row, col] += mat_phi[row, j]
                if mat_holes_phi[row, j] != 0.0:
                    mat_holes_s[row, col] = 1.0

    plt.figure()
    plt.imshow(
        mat_s, aspect="auto", origin="lower",
        extent=[s_extent_min, s_extent_max,
                float(layer_ids_sorted[0]) - 0.5, float(layer_ids_sorted[-1]) + 0.5],
    )
    plt.xlabel("s = phi * R_layer [mm]")
    plt.ylabel("Layer index")
    plt.title("Counts vs s=phi*R and layer" + (" (unique)" if unique else ""))
    plt.colorbar(label="Counts")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_layer_s_2d.png", dpi=200)
    plt.close()

    plt.figure()
    plt.imshow(
        mat_holes_s, aspect="auto", origin="lower",
        extent=[s_extent_min, s_extent_max,
                float(layer_ids_sorted[0]) - 0.5, float(layer_ids_sorted[-1]) + 0.5],
        vmin=0.0, vmax=1.0
    )
    plt.xlabel("s = phi * R_layer [mm]")
    plt.ylabel("Layer index")
    plt.title(f"Holes (large only) vs s=phi*R and layer (>{max_allowed_hole_mm} mm)")
    plt.colorbar(label="Hole flag (0/1)")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_layer_holes_s_2d.png", dpi=200)
    plt.close()

    for lay in layer_ids_sorted:
        R = layer_radii_mm.get(lay, None)
        if R is None or R <= 0:
            continue
        row = lay_to_row[lay]
        x_s = phi_vals * float(R)
        y = mat_phi[row, :]

        plt.figure()
        plt.plot(x_s, y)
        plt.xlabel("s = phi * R_layer [mm]")
        plt.ylabel("Counts")
        plt.title(f"TRKLayer{lay}: counts vs s=phi*R" + (" (unique)" if unique else ""))
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_L{lay}_s.png", dpi=200)
        plt.close()

        plt.figure()
        plt.plot(x_s, mat_holes_phi[row, :], color="red")
        plt.ylim(-0.1, 1.1)
        plt.xlabel("s = phi * R_layer [mm]")
        plt.ylabel("Hole flag (1=large hole)")
        plt.title(f"TRKLayer{lay}: large holes vs s (>{max_allowed_hole_mm} mm)")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_L{lay}_holes_s.png", dpi=200)
        plt.close()


# -------------------------
# Output directory helper
# -------------------------

def make_output_dir(gdml_path: str) -> str:
    base = os.path.basename(gdml_path)
    stem = os.path.splitext(base)[0]
    out_dir = f"ouput_overlap_analysis_for_{stem}"
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


# -------------------------
# Main
# -------------------------

def main():
    ap = argparse.ArgumentParser(description="Scan GDML for sensitive-volume intersections vs phi at fixed eta.")
    ap.add_argument("gdml", help="Path to GDML file")

    ap.add_argument("--eta", type=float, required=True, help="Fixed pseudorapidity eta")
    ap.add_argument("--phi-step", type=float, default=0.000625, help="Phi step in radians (default 0.000625)")
    ap.add_argument("--phi-min", type=float, default=0.0, help="Phi min (default 0)")
    ap.add_argument("--phi-max", type=float, default=2.0 * math.pi, help="Phi max (default 2*pi)")

    ap.add_argument(
        "--materials", nargs="+",
        default=["TRK_SILICON", "FT3_SILICON", "FT3_Silicon", "TF3_Silicon"],
        help="Sensitive materials list (space-separated).",
    )

    ap.add_argument(
        "--start", nargs=3, type=float, default=(0.0, 0.0, 0.0), metavar=("X", "Y", "Z"),
        help="Start point x y z in GDML units (cm). Default 0 0 0."
    )
    ap.add_argument("--eps", type=float, default=1e-4, help="Nudge beyond boundary in GDML units (cm). Default 1e-4.")
    ap.add_argument("--max-steps", type=int, default=200000, help="Max navigation iterations per phi (default 200000)")
    ap.add_argument("--max-track-len", type=float, default=None, help="Optional max track length in GDML units (cm)")
    ap.add_argument("--unique", action="store_true", help="Count unique sensitive volumes instead of entry crossings")

    ap.add_argument("--no-show", action="store_true", help="Do not show interactive plot window")
    ap.add_argument("--debug", action="store_true", help="Print debug info to stderr for problematic navigation")

    ap.add_argument("--freq-decimals", type=int, default=2,
                    help="Decimals in printed frequency percentages (default 2)")
    ap.add_argument("--freq-top", type=int, default=0, help="If >0, print only top-N most frequent N values")

    ap.add_argument("--splitIntoRadialLayers", action="store_true",
                    help="Enable layer assignment via TRKLayerX in paths.")
    ap.add_argument("--maxAllowedHoleSize", type=float, default=0.22,
                    help="Max allowed hole size in mm. Default 0.22 mm")

    # Optional filenames (will still be placed inside output dir)
    ap.add_argument("--root", default=None, help="ROOT output filename (default: overlap_scan.root)")
    ap.add_argument("--csv", default=None, help="CSV output filename (default: scan.csv)")
    ap.add_argument("--png", default=None, help="Global N(phi) PNG filename (default: N_vs_phi.png)")

    args = ap.parse_args()

    _ensure_matplotlib_backend(args.no_show)

    out_dir = make_output_dir(args.gdml)

    root_name = os.path.basename(args.root) if args.root else "overlap_scan.root"
    csv_name = os.path.basename(args.csv) if args.csv else "scan.csv"
    png_name = os.path.basename(args.png) if args.png else "N_vs_phi.png"

    root_path = os.path.join(out_dir, root_name)
    csv_path = os.path.join(out_dir, csv_name)
    global_png_path = os.path.join(out_dir, png_name)

    gGeoManager = load_gdml_into_root(args.gdml)

    phi_vals = np.arange(args.phi_min, args.phi_max + 0.5 * args.phi_step, args.phi_step, dtype=float)
    counts = np.zeros_like(phi_vals, dtype=int)

    layer_counts_by_phi: List[Dict[int, int]] = []
    layer_radii_mm_accum: Dict[int, List[float]] = {}

    for i, phi in enumerate(phi_vals):
        total, lay_counts, lay_r_est_mm = trace_phi(
            gGeoManager=gGeoManager,
            eta=args.eta,
            phi=float(phi),
            sensitive_materials=list(args.materials),
            start_xyz=tuple(args.start),
            eps=args.eps,
            max_steps=args.max_steps,
            max_track_len=args.max_track_len,
            unique=args.unique,
            split_into_layers=args.splitIntoRadialLayers,
            debug=args.debug,
        )
        counts[i] = int(total)

        if args.splitIntoRadialLayers:
            layer_counts_by_phi.append(lay_counts or {})
            if lay_r_est_mm:
                for lay, rmm in lay_r_est_mm.items():
                    layer_radii_mm_accum.setdefault(lay, []).append(float(rmm))

    print_count_frequencies(counts, decimals=args.freq_decimals, top=args.freq_top)

    # CSV
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["phi_rad", "n_sensitive"])
        for p, c in zip(phi_vals, counts):
            w.writerow([f"{p:.10f}", int(c)])
    print(f"Wrote CSV: {csv_path}")

    # ROOT basic
    save_root_histograms_basic(root_path, phi_vals, counts, unique=args.unique)

    # Global plot N(phi)
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(phi_vals, counts)
    plt.xlabel("phi [rad]")
    ylabel = "N unique sensitive volumes" if args.unique else "N entries into sensitive volumes"
    plt.ylabel(ylabel)
    plt.title(f"Sensitive intersections vs phi at eta={args.eta:g}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(global_png_path, dpi=200, bbox_inches="tight")
    print(f"Wrote plot: {global_png_path}")
    if not args.no_show and not args.splitIntoRadialLayers:
        plt.show()
    plt.close()

    # Export h_sensitive_freq to PNG (lin/log)
    freq_lin = os.path.join(out_dir, "h_sensitive_freq_lin.png")
    freq_log = os.path.join(out_dir, "h_sensitive_freq_log.png")
    root_export_hist_png(root_path, "h_sensitive_freq", freq_lin, logy=False)
    root_export_hist_png(root_path, "h_sensitive_freq", freq_log, logy=True)
    print(f"Wrote ROOT-exported PNGs: {freq_lin}, {freq_log}")

    # If layers enabled: append layer histograms + layer PNGs + per-layer hole-width PNGs
    layers_written: List[int] = []

    if args.splitIntoRadialLayers:
        layer_ids = sorted({lay for d in layer_counts_by_phi for lay in d.keys()})
        if not layer_ids:
            print('WARNING: No hits assigned to TRKLayerX (check volume/node names in geometry paths).', file=sys.stderr)
        else:
            layer_radii_mm: Dict[int, float] = {}
            for lay, rs in layer_radii_mm_accum.items():
                if rs:
                    layer_radii_mm[lay] = float(sum(rs) / len(rs))

            layers_written = save_root_histograms_layers(
                root_filename=root_path,
                phi_vals=phi_vals,
                layer_ids_sorted=layer_ids,
                layer_counts_by_phi=layer_counts_by_phi,
                layer_radii_mm=layer_radii_mm,
                phi_step=args.phi_step,
                max_allowed_hole_mm=args.maxAllowedHoleSize,
                unique=args.unique,
            )

            layer_prefix = os.path.join(out_dir, "layers")
            save_pngs_layers(
                out_prefix=layer_prefix,
                phi_vals=phi_vals,
                layer_ids_sorted=layer_ids,
                layer_counts_by_phi=layer_counts_by_phi,
                layer_radii_mm=layer_radii_mm,
                phi_step=args.phi_step,
                max_allowed_hole_mm=args.maxAllowedHoleSize,
                unique=args.unique,
            )
            print(f"Wrote layer PNGs with prefix: {layer_prefix}_*.png")

            # Export global hole-width distribution (lin/log)
            holew_lin = os.path.join(out_dir, "h_large_hole_width_mm_lin.png")
            holew_log = os.path.join(out_dir, "h_large_hole_width_mm_log.png")
            root_export_hist_png(root_path, "h_large_hole_width_mm", holew_lin, logy=False)
            root_export_hist_png(root_path, "h_large_hole_width_mm", holew_log, logy=True)
            print(f"Wrote ROOT-exported PNGs: {holew_lin}, {holew_log}")

            # Export per-layer hole-width distributions (lin/log)
            per_layer_dir = os.path.join(out_dir, "hole_width_distributions_per_layer")
            os.makedirs(per_layer_dir, exist_ok=True)

            for lay in layers_written:
                hname = f"h_large_hole_width_mm_L{lay}"
                out_lin = os.path.join(per_layer_dir, f"{hname}_lin.png")
                out_log = os.path.join(per_layer_dir, f"{hname}_log.png")
                # always exists (we write empty if no holes)
                root_export_hist_png(root_path, hname, out_lin, logy=False)
                root_export_hist_png(root_path, hname, out_log, logy=True)

            print(f"Wrote per-layer hole-width PNGs into: {per_layer_dir}")

    print(f"\nAll outputs saved in directory:\n  {out_dir}\n")


if __name__ == "__main__":
    main()
