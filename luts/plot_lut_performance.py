#!/usr/bin/env python3
"""Plot ALICE 3 LUT performance and compare its pT resolution with QA.

The first input is the binary ``lut_pdg*.dat`` format defined by
``ALICE3/Core/FlatLutEntry.h`` (LUTCOVM_VERSION = 20210801). By default the
script also reads ``/tmp/lutMakerQA.root`` and overlays the QA resolution for
the PDG code stored in the LUT header.

The LUT relative transverse-momentum resolution follows
``TrackSmearer::getPtRes``::

    sigma(pT) / pT = sqrt(covm[14]) * pT

where ``covm[14]`` is the variance of q/pT stored in the LUT entry.  The QA
curve is extracted directly from the three-dimensional histogram

    res_qopT_res_vs_pT_pdg<PDG>

using exactly the same eta-bin window selected from the LUT.  In each generated
pT bin the width of Delta(q/pT) is converted to relative pT resolution through

    sigma(pT) / pT = pT * sigma(Delta(q/pT))

for the unit-charge species supported by the LUT maker.

NumPy is used to decode and manipulate the binary LUT payload. PyROOT is used
for reading the QA ROOT file and drawing the plots.
"""

from __future__ import annotations

import argparse
import math
import struct
import sys
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


LUT_VERSION = 20210801
HEADER_FORMAT = "<iiff" + "iff?3x" * 4
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
ENTRY_SIZE = 312

# Native C++ layout of o2::delphes::lutEntry_t. Explicit offsets make the
# reader independent of NumPy's own structure-alignment choices.
ENTRY_DTYPE = np.dtype(
    {
        "names": [
            "nch",
            "eta",
            "pt",
            "valid",
            "eff",
            "eff2",
            "itof",
            "otof",
            "covm",
            "eigval",
            "eigvec",
            "eiginv",
        ],
        "formats": [
            "<f4",
            "<f4",
            "<f4",
            "?",
            "<f4",
            "<f4",
            "<f4",
            "<f4",
            ("<f4", (15,)),
            ("<f4", (5,)),
            ("<f4", (5, 5)),
            ("<f4", (5, 5)),
        ],
        "offsets": [0, 4, 8, 12, 16, 20, 24, 28, 32, 92, 112, 212],
        "itemsize": ENTRY_SIZE,
    }
)


@dataclass(frozen=True)
class AxisMap:
    nbins: int
    minimum: float
    maximum: float
    logarithmic: bool

    def centers(self) -> np.ndarray:
        values = self.minimum + (np.arange(self.nbins, dtype=float) + 0.5) * (
            self.maximum - self.minimum
        ) / self.nbins
        return np.power(10.0, values) if self.logarithmic else values

    def edges(self) -> np.ndarray:
        values = self.minimum + np.arange(self.nbins + 1, dtype=float) * (
            self.maximum - self.minimum
        ) / self.nbins
        return np.power(10.0, values) if self.logarithmic else values


@dataclass(frozen=True)
class LutHeader:
    version: int
    pdg: int
    mass: float
    field: float
    nch: AxisMap
    radius: AxisMap
    eta: AxisMap
    pt: AxisMap

    @property
    def shape(self) -> tuple[int, int, int, int]:
        return (
            self.nch.nbins,
            self.radius.nbins,
            self.eta.nbins,
            self.pt.nbins,
        )

    @property
    def number_of_entries(self) -> int:
        return math.prod(self.shape)


@dataclass(frozen=True)
class LutData:
    header: LutHeader
    entries: np.ndarray


@dataclass(frozen=True)
class QaResolutionCurve:
    source: str
    histogram_name: str
    label: str
    pt: np.ndarray
    resolution: np.ndarray


@dataclass
class RootPlotObjects:
    """Keep ROOT objects alive until the canvas has been saved."""

    canvas: object
    frame: object
    graphs: list[object]
    legend: object | None
    labels: list[object]


_ROOT_OBJECT_COUNTER = 0


def unique_root_name(prefix: str) -> str:
    global _ROOT_OBJECT_COUNTER
    _ROOT_OBJECT_COUNTER += 1
    return f"{prefix}_{_ROOT_OBJECT_COUNTER}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot ALICE 3 LUT performance and compare pT resolution with LutMaker QA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("lut", type=Path, help="Input lut_pdg*.dat file")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("lut_performance"),
        help="Directory for output plots",
    )
    parser.add_argument(
        "--eta-min",
        type=float,
        default=-0.8,
        help="Minimum eta-bin centre to include",
    )
    parser.add_argument(
        "--eta-max",
        type=float,
        default=0.8,
        help="Maximum eta-bin centre to include",
    )
    parser.add_argument(
        "--nch-bin",
        type=int,
        default=0,
        help="Charged-particle-density bin index",
    )
    parser.add_argument(
        "--radius-bin",
        type=int,
        default=0,
        help="Reference-radius bin index",
    )
    parser.add_argument(
        "--efficiency-field",
        choices=("eff", "eff2"),
        default="eff",
        help="Efficiency field to plot",
    )
    parser.add_argument(
        "--pt-min",
        type=float,
        default=None,
        help="Minimum displayed pT",
    )
    parser.add_argument(
        "--pt-max",
        type=float,
        default=None,
        help="Maximum displayed pT",
    )
    parser.add_argument(
        "--linear-x",
        action="store_true",
        help="Use a linear rather than logarithmic pT axis",
    )
    parser.add_argument(
        "--fractional-resolution",
        action="store_true",
        help="Plot sigma(pT)/pT as a fraction rather than percent",
    )
    parser.add_argument(
        "--average-only",
        action="store_true",
        help="Plot only the average over selected eta bins",
    )
    parser.add_argument(
        "--pdf",
        action="store_true",
        help="Also save PDF versions",
    )
    parser.add_argument(
        "-b",
        "--background",
        action="store_true",
        help=(
            "Run in ROOT batch mode and save immediately after all plots have "
            "been drawn, without waiting for Enter"
        ),
    )
    parser.add_argument(
        "--require-resolution",
        action="store_true",
        help="Fail if covm[14] contains no positive resolution values",
    )
    parser.add_argument(
        "--qa-file",
        type=Path,
        default=Path("/tmp/lutMakerQA.root"),
        help="LutMaker QA ROOT file",
    )
    return parser.parse_args()


def import_root(background: bool) -> object:
    try:
        import ROOT  # type: ignore
    except ImportError as error:
        raise RuntimeError(
            "PyROOT is required for plotting. Run the script in an environment "
            "where ROOT is available."
        ) from error

    ROOT.gROOT.SetBatch(background)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFillColor(0)
    ROOT.gStyle.SetEndErrorSize(0)
    return ROOT


def _axis_map(values: tuple[object, ...], offset: int) -> AxisMap:
    return AxisMap(
        nbins=int(values[offset]),
        minimum=float(values[offset + 1]),
        maximum=float(values[offset + 2]),
        logarithmic=bool(values[offset + 3]),
    )


def read_lut(path: Path) -> LutData:
    if not path.is_file():
        raise FileNotFoundError(f"LUT file does not exist: {path}")

    raw = path.read_bytes()
    if len(raw) < HEADER_SIZE:
        raise RuntimeError(
            f"File is too small for a LUT header: {len(raw)} < {HEADER_SIZE} bytes"
        )

    values = struct.unpack_from(HEADER_FORMAT, raw, 0)
    header = LutHeader(
        version=int(values[0]),
        pdg=int(values[1]),
        mass=float(values[2]),
        field=float(values[3]),
        nch=_axis_map(values, 4),
        radius=_axis_map(values, 8),
        eta=_axis_map(values, 12),
        pt=_axis_map(values, 16),
    )

    if header.version != LUT_VERSION:
        raise RuntimeError(
            f"Unsupported LUT version {header.version}; expected {LUT_VERSION}"
        )

    expected_size = HEADER_SIZE + header.number_of_entries * ENTRY_SIZE
    if len(raw) != expected_size:
        raise RuntimeError(
            "Unexpected LUT size: "
            f"got {len(raw)} bytes, expected {expected_size} bytes from the header "
            f"({header.number_of_entries} entries x {ENTRY_SIZE} bytes)"
        )

    entries = np.frombuffer(
        raw,
        dtype=ENTRY_DTYPE,
        count=header.number_of_entries,
        offset=HEADER_SIZE,
    ).reshape(header.shape)
    return LutData(header=header, entries=entries)


def _half_open_axis_bin_range(
    axis: object,
    minimum: float,
    maximum: float,
    label: str,
) -> tuple[int, int]:
    """Return regular ROOT bins implementing minimum <= value < maximum."""
    if not minimum < maximum:
        raise ValueError(f"Invalid {label} interval [{minimum}, {maximum}]")

    axis_minimum = float(axis.GetXmin())
    axis_maximum = float(axis.GetXmax())
    clipped_minimum = max(minimum, axis_minimum)
    clipped_maximum = min(maximum, axis_maximum)
    if not clipped_minimum < clipped_maximum:
        raise ValueError(
            f"Requested {label} interval [{minimum}, {maximum}] lies outside "
            f"the histogram range [{axis_minimum}, {axis_maximum}]"
        )

    scale = max(abs(clipped_minimum), abs(clipped_maximum), 1.0)
    epsilon = 1.0e-10 * scale
    first = int(axis.FindFixBin(clipped_minimum + epsilon))
    last = int(axis.FindFixBin(clipped_maximum - epsilon))
    first = max(1, min(first, int(axis.GetNbins())))
    last = max(1, min(last, int(axis.GetNbins())))
    if first > last:
        raise ValueError(
            f"No regular histogram bins selected for {label} interval "
            f"[{minimum}, {maximum}]"
        )
    return first, last


def _qopt_resolution_from_th3(
    histogram: object,
    eta_min: float,
    eta_max: float,
    pt_min: float | None,
    pt_max: float | None,
    fractional_resolution: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract pT*sigma(Delta(q/pT)) in the requested eta interval.

    The LutMaker QA histogram is expected to have:
      X = generated pT,
      Y = Delta(q/pT),
      Z = generated eta.

    The projection/profile is done in ROOT so this reproduces the QA width
    extraction used by LutMaker, but after restricting the Z axis to the same
    eta-bin edges selected from the LUT.
    """
    if not histogram.InheritsFrom("TH3"):
        raise RuntimeError(
            f"Object '{histogram.GetName()}' is a {histogram.ClassName()}, "
            "while an eta-resolved TH3 is required"
        )

    working = histogram.Clone(unique_root_name("qa_qopt_source"))
    working.SetDirectory(0)

    eta_axis = working.GetZaxis()
    eta_first, eta_last = _half_open_axis_bin_range(
        eta_axis, eta_min, eta_max, "eta"
    )
    eta_axis.SetRange(eta_first, eta_last)

    # X = pT, Y = Delta(q/pT), Z = eta.  Projecting "yx" after restricting Z
    # produces the pT-versus-residual TH2 used for ProfileX.
    projection = working.Project3D("yx")
    if not projection:
        raise RuntimeError(
            f"Project3D('yx') failed for histogram '{histogram.GetName()}'"
        )
    projection.SetName(unique_root_name("qa_qopt_projection"))
    projection.SetDirectory(0)

    profile = projection.ProfileX(
        unique_root_name("qa_qopt_profile"),
        1,
        -1,
        "s",
    )
    if not profile:
        raise RuntimeError(
            f"ProfileX failed for histogram '{histogram.GetName()}'"
        )
    profile.SetDirectory(0)

    pt_values: list[float] = []
    resolutions: list[float] = []
    entries_per_bin: list[float] = []

    for pt_bin in range(1, int(profile.GetNbinsX()) + 1):
        pt = float(profile.GetXaxis().GetBinCenter(pt_bin))
        if not math.isfinite(pt) or pt <= 0.0:
            continue
        if pt_min is not None and pt < pt_min:
            continue
        if pt_max is not None and pt > pt_max:
            continue

        entries = float(profile.GetBinEntries(pt_bin))
        sigma_q_over_pt = float(profile.GetBinError(pt_bin))
        if (
            entries <= 1.0
            or not math.isfinite(sigma_q_over_pt)
            or sigma_q_over_pt <= 0.0
        ):
            continue

        # Invert q/pT to the linearized relative pT resolution:
        # sigma(pT)/pT ~= pT * sigma(Delta(q/pT)).
        relative_pt_resolution = pt * sigma_q_over_pt
        if not fractional_resolution:
            relative_pt_resolution *= 100.0

        if not math.isfinite(relative_pt_resolution) or relative_pt_resolution <= 0.0:
            continue

        pt_values.append(pt)
        resolutions.append(relative_pt_resolution)
        entries_per_bin.append(entries)

    return (
        np.asarray(pt_values, dtype=float),
        np.asarray(resolutions, dtype=float),
        np.asarray(entries_per_bin, dtype=float),
    )


def read_qa_resolution(
    ROOT: object,
    path: Path,
    pdg: int,
    eta_min: float,
    eta_max: float,
    pt_min: float | None,
    pt_max: float | None,
    fractional_resolution: bool,
) -> list[QaResolutionCurve]:
    if not path.is_file():
        raise FileNotFoundError(f"QA ROOT file does not exist: {path}")

    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open QA ROOT file: {path}")

    histogram_name = f"res_qopT_res_vs_pT_pdg{pdg}"
    try:
        histogram = root_file.Get(histogram_name)
        if not histogram:
            root_file.ls()
            raise RuntimeError(
                f"Histogram '{histogram_name}' was not found in {path}"
            )

        pt, resolution, entries = _qopt_resolution_from_th3(
            histogram,
            eta_min,
            eta_max,
            pt_min,
            pt_max,
            fractional_resolution,
        )
    finally:
        root_file.Close()

    if pt.size == 0:
        raise RuntimeError(
            f"No usable q/pT QA resolution points were found in {path} for "
            f"PDG {pdg}, {eta_min:g} <= eta < {eta_max:g}"
        )

    print(
        f"QA eta selection: {eta_min:g} <= eta < {eta_max:g}; "
        f"{pt.size} pT bins, {int(np.sum(entries))} selected entries"
    )

    return [
        QaResolutionCurve(
            source="qopt",
            histogram_name=histogram_name,
            label=(
                "QA, p_{T} #times #sigma(#Delta(q/p_{T})), "
                f"{eta_min:g} #leq #eta < {eta_max:g}"
            ),
            pt=pt,
            resolution=resolution,
        )
    ]


def validate_bin_index(index: int, number_of_bins: int, label: str) -> None:
    if not 0 <= index < number_of_bins:
        raise ValueError(
            f"{label} index {index} is outside [0, {number_of_bins - 1}]"
        )


def selected_eta_indices(header: LutHeader, eta_min: float, eta_max: float) -> np.ndarray:
    if eta_max < eta_min:
        raise ValueError("--eta-max must be greater than or equal to --eta-min")

    eta = header.eta.centers()
    selected = np.flatnonzero((eta >= eta_min) & (eta <= eta_max))
    if selected.size == 0:
        nearest = int(np.argmin(np.abs(eta - 0.5 * (eta_min + eta_max))))
        raise ValueError(
            f"No eta-bin centres lie in [{eta_min:g}, {eta_max:g}]. "
            f"The nearest available centre is {eta[nearest]:g}."
        )
    return selected


def finite_pt_mask(pt: np.ndarray, pt_min: float | None, pt_max: float | None) -> np.ndarray:
    mask = np.isfinite(pt) & (pt > 0.0)
    if pt_min is not None:
        mask &= pt >= pt_min
    if pt_max is not None:
        mask &= pt <= pt_max
    if not np.any(mask):
        raise ValueError("The requested pT range contains no LUT bins")
    return mask


def finite_average(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    counts = np.sum(finite, axis=0)
    sums = np.sum(np.where(finite, values, 0.0), axis=0)
    average = np.full(values.shape[1], np.nan, dtype=float)
    np.divide(sums, counts, out=average, where=counts > 0)
    return average


def graph_from_arrays(
    ROOT: object,
    x_values: np.ndarray,
    y_values: np.ndarray,
    name: str,
    color: int,
    line_width: int,
    line_style: int = 1,
    marker_style: int = 0,
    marker_size: float = 1.0,
) -> object | None:
    mask = np.isfinite(x_values) & np.isfinite(y_values)
    if not np.any(mask):
        return None

    x = array("d", np.asarray(x_values[mask], dtype=np.float64).tolist())
    y = array("d", np.asarray(y_values[mask], dtype=np.float64).tolist())
    graph = ROOT.TGraph(len(x), x, y)
    graph.SetName(unique_root_name(name))
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    graph.SetLineWidth(line_width)
    graph.SetLineStyle(line_style)
    graph.SetMarkerStyle(marker_style)
    graph.SetMarkerSize(marker_size)
    return graph


def curve_colors(ROOT: object) -> Sequence[int]:
    return (
        ROOT.kAzure + 1,
        ROOT.kOrange + 7,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kCyan + 2,
        ROOT.kRed + 1,
        ROOT.kViolet + 1,
        ROOT.kTeal + 3,
        ROOT.kPink + 7,
        ROOT.kGray + 2,
    )


def make_canvas(ROOT: object, stem: str, logarithmic_x: bool) -> object:
    canvas = ROOT.TCanvas(unique_root_name(f"canvas_{stem}"), stem, 850, 680)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.045)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.09)
    canvas.SetTicks(1, 1)
    canvas.SetGrid(1, 1)
    canvas.SetLogx(logarithmic_x)
    return canvas


def make_frame(
    ROOT: object,
    stem: str,
    pt: np.ndarray,
    y_min: float,
    y_max: float,
    y_title: str,
) -> object:
    """Create the axes with a ROOT TH1D frame."""
    x_min = float(np.min(pt))
    x_max = float(np.max(pt))
    if x_max <= x_min:
        x_max = x_min * 1.01 if x_min > 0.0 else x_min + 1.0

    frame = ROOT.TH1D(
        unique_root_name(f"frame_{stem}"),
        f";p_{{T}} (GeV/#it{{c}});{y_title}",
        100,
        x_min,
        x_max,
    )
    frame.SetDirectory(0)
    frame.SetMinimum(y_min)
    frame.SetMaximum(y_max)
    frame.SetStats(False)
    frame.GetXaxis().SetTitleSize(0.052)
    frame.GetYaxis().SetTitleSize(0.052)
    frame.GetXaxis().SetLabelSize(0.044)
    frame.GetYaxis().SetLabelSize(0.044)
    frame.GetXaxis().SetTitleOffset(1.08)
    frame.GetYaxis().SetTitleOffset(1.20)
    frame.Draw("AXIS")
    return frame


def make_title_labels(ROOT: object, header: LutHeader, title: str) -> list[object]:
    title_label = ROOT.TLatex()
    title_label.SetNDC(True)
    title_label.SetTextFont(42)
    title_label.SetTextSize(0.043)
    title_label.DrawLatex(0.15, 0.935, title)

    info_label = ROOT.TLatex()
    info_label.SetNDC(True)
    info_label.SetTextFont(42)
    info_label.SetTextSize(0.034)
    info_label.SetTextAlign(31)
    info_label.DrawLatex(0.95, 0.935, f"PDG {header.pdg}, B = {header.field:g} T")
    return [title_label, info_label]


def draw_curves(
    ROOT: object,
    pt: np.ndarray,
    values: np.ndarray,
    eta_centers: np.ndarray,
    average_only: bool,
    legend: object,
    name_prefix: str,
    legend_prefix: str = "",
) -> list[object]:
    graphs: list[object] = []
    colors = curve_colors(ROOT)

    # With one selected eta bin, drawing both the individual and average curve
    # would duplicate exactly the same information.
    draw_individual = not average_only and eta_centers.size > 1
    if draw_individual:
        for index, (eta, curve) in enumerate(zip(eta_centers, values, strict=True)):
            graph = graph_from_arrays(
                ROOT,
                pt,
                curve,
                f"{name_prefix}_eta_{index}",
                colors[index % len(colors)],
                line_width=2,
            )
            if graph is None:
                continue
            graph.Draw("L SAME")
            legend.AddEntry(graph, f"{legend_prefix}#eta = {eta:g}", "l")
            graphs.append(graph)

    average = finite_average(values)
    average_label = (
        f"{legend_prefix}average, {eta_centers.min():g} #leq #eta #leq {eta_centers.max():g}"
        if eta_centers.size > 1
        else f"{legend_prefix}#eta = {eta_centers[0]:g}"
    )
    average_graph = graph_from_arrays(
        ROOT,
        pt,
        average,
        f"{name_prefix}_average",
        ROOT.kBlack,
        line_width=4,
    )
    if average_graph is not None:
        average_graph.Draw("L SAME")
        legend.AddEntry(average_graph, average_label, "l")
        graphs.append(average_graph)

    return graphs


def update_canvas(ROOT: object, canvas: object) -> None:
    """Finish drawing one canvas without writing any files."""
    canvas.RedrawAxis()
    canvas.Modified()
    canvas.Update()
    if not ROOT.gROOT.IsBatch():
        ROOT.gSystem.ProcessEvents()


def save_canvas(
    canvas: object,
    output_dir: Path,
    stem: str,
    save_pdf: bool,
) -> None:
    """Write an already-drawn canvas to disk."""
    canvas.SaveAs(str(output_dir / f"{stem}.png"))
    if save_pdf:
        canvas.SaveAs(str(output_dir / f"{stem}.pdf"))


def wait_for_confirmation(ROOT: object) -> None:
    """Keep interactive ROOT canvases visible until the user presses Enter."""
    print(
        "All plots have been drawn. Press Enter to save the files and exit...",
        flush=True,
    )
    try:
        input()
    except EOFError as error:
        raise RuntimeError(
            "Interactive confirmation requires standard input. "
            "Run with --background for unattended execution."
        ) from error
    ROOT.gSystem.ProcessEvents()


def plot_efficiency(
    ROOT: object,
    pt: np.ndarray,
    efficiencies: np.ndarray,
    eta_centers: np.ndarray,
    header: LutHeader,
    args: argparse.Namespace,
) -> RootPlotObjects:
    canvas = make_canvas(ROOT, "efficiency_vs_pt", not args.linear_x)
    frame = make_frame(
        ROOT,
        "efficiency_vs_pt",
        pt,
        0.0,
        1.05,
        "Tracking efficiency",
    )

    legend = ROOT.TLegend(0.17, 0.68, 0.55, 0.87)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.032)
    if not args.average_only and eta_centers.size > 5:
        legend.SetNColumns(2)
        legend.SetX2NDC(0.72)

    graphs = draw_curves(
        ROOT,
        pt,
        efficiencies,
        eta_centers,
        args.average_only,
        legend,
        "efficiency",
    )
    labels = make_title_labels(ROOT, header, "LUT tracking efficiency")
    legend.Draw()
    update_canvas(ROOT, canvas)
    return RootPlotObjects(canvas, frame, graphs, legend, labels)


def plot_resolution(
    ROOT: object,
    pt: np.ndarray,
    resolutions: np.ndarray,
    qa_curves: Sequence[QaResolutionCurve],
    eta_centers: np.ndarray,
    header: LutHeader,
    args: argparse.Namespace,
) -> tuple[bool, RootPlotObjects]:
    finite_positive = np.isfinite(resolutions) & (resolutions > 0.0)
    has_lut_resolution = bool(np.any(finite_positive))
    plot_values = np.where(finite_positive, resolutions, np.nan)

    finite_values: list[np.ndarray] = []
    if has_lut_resolution:
        finite_values.append(plot_values[np.isfinite(plot_values)])
    finite_values.extend(
        curve.resolution[np.isfinite(curve.resolution)] for curve in qa_curves
    )
    finite_values = [values for values in finite_values if values.size > 0]

    if finite_values:
        y_max = max(float(np.max(values)) for values in finite_values)
        y_max = max(1.15 * y_max, 0.01 if args.fractional_resolution else 0.1)
    else:
        y_max = 1.0

    frame_pt_parts = [pt]
    frame_pt_parts.extend(curve.pt for curve in qa_curves if curve.pt.size > 0)
    frame_pt = np.concatenate(frame_pt_parts)

    unit = "" if args.fractional_resolution else " (%)"
    canvas = make_canvas(ROOT, "pt_resolution_lut_vs_qa", not args.linear_x)
    frame = make_frame(
        ROOT,
        "pt_resolution_lut_vs_qa",
        frame_pt,
        0.0,
        y_max,
        f"#sigma(p_{{T}})/p_{{T}}{unit}",
    )
    graphs: list[object] = []
    legend = None

    if has_lut_resolution or qa_curves:
        legend = ROOT.TLegend(0.16, 0.63, 0.69, 0.87)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.030)
        if not args.average_only and eta_centers.size > 5:
            legend.SetNColumns(2)
            legend.SetX2NDC(0.87)

    if has_lut_resolution and legend is not None:
        graphs.extend(
            draw_curves(
                ROOT,
                pt,
                plot_values,
                eta_centers,
                args.average_only,
                legend,
                "lut_resolution",
                legend_prefix="LUT, ",
            )
        )

    qa_styles = {
        "qopt": (ROOT.kRed + 1, 2, 24),
    }
    if legend is not None:
        for curve in qa_curves:
            color, line_style, marker_style = qa_styles[curve.source]
            graph = graph_from_arrays(
                ROOT,
                curve.pt,
                curve.resolution,
                f"qa_resolution_{curve.source}",
                color,
                line_width=3,
                line_style=line_style,
                marker_style=marker_style,
                marker_size=0.85,
            )
            if graph is None:
                continue
            graph.Draw("LP SAME")
            legend.AddEntry(graph, curve.label, "lp")
            graphs.append(graph)

    labels = make_title_labels(
        ROOT,
        header,
        "LUT and q/p_{T}-based QA resolution",
    )

    if legend is not None and graphs:
        legend.Draw()
    else:
        warning = ROOT.TLatex()
        warning.SetNDC(True)
        warning.SetTextAlign(22)
        warning.SetTextFont(62)
        warning.SetTextSize(0.045)
        warning.DrawLatex(0.54, 0.55, "No non-zero p_{T}-resolution values")
        labels.append(warning)

    update_canvas(ROOT, canvas)
    return has_lut_resolution, RootPlotObjects(
        canvas, frame, graphs, legend, labels
    )


def print_summary(
    lut: LutData,
    eta_indices: np.ndarray,
    qa_curves: Sequence[QaResolutionCurve],
    args: argparse.Namespace,
) -> None:
    header = lut.header
    eta_centers = header.eta.centers()[eta_indices]
    print(f"LUT:       {args.lut}")
    print(f"QA:        {args.qa_file}")
    print(f"Version:   {header.version}")
    print(f"PDG:       {header.pdg}")
    print(f"Mass:      {header.mass:g} GeV/c^2")
    print(f"Field:     {header.field:g} T")
    print(
        f"Shape:     nch={header.nch.nbins}, radius={header.radius.nbins}, "
        f"eta={header.eta.nbins}, pT={header.pt.nbins}"
    )
    print("Eta bins:  " + ", ".join(f"{value:g}" for value in eta_centers))
    for curve in qa_curves:
        print(
            f"QA curve:  {curve.histogram_name} "
            f"({curve.pt.size} positive bins)"
        )
    print(
        "NOTE: QA uses the eta-bin edges of exactly the LUT bins selected above "
        "and only the q/pT residual, converted as pT*sigma(Delta(q/pT))."
    )


def main() -> int:
    args = parse_args()
    try:
        ROOT = import_root(args.background)
        lut = read_lut(args.lut)
        validate_bin_index(args.nch_bin, lut.header.nch.nbins, "nch-bin")
        validate_bin_index(args.radius_bin, lut.header.radius.nbins, "radius-bin")
        eta_indices = selected_eta_indices(lut.header, args.eta_min, args.eta_max)

        pt_all = lut.header.pt.centers()
        pt_mask = finite_pt_mask(pt_all, args.pt_min, args.pt_max)
        pt = pt_all[pt_mask]

        selected = lut.entries[
            args.nch_bin,
            args.radius_bin,
            eta_indices,
            :,
        ]
        efficiencies = np.asarray(
            selected[args.efficiency_field], dtype=float
        )[:, pt_mask]

        # TrackSmearer::getPtRes returns sqrt(covm[14]) * pT, i.e. the
        # relative pT resolution. Negative covariance entries are invalid.
        q_over_pt_variance = np.asarray(
            selected["covm"][..., 14], dtype=float
        )[:, pt_mask]
        valid_entries = np.asarray(selected["valid"], dtype=bool)[:, pt_mask]
        q_over_pt_variance = np.where(
            valid_entries & (q_over_pt_variance >= 0.0),
            q_over_pt_variance,
            np.nan,
        )
        resolutions = np.sqrt(q_over_pt_variance) * pt[np.newaxis, :]
        if not args.fractional_resolution:
            resolutions *= 100.0

        # Use the physical edges of the selected LUT eta bins for QA, rather
        # than merely reusing the centre-selection limits.  For the default
        # LUT this maps centres {-0.4, +0.4} to the exact window [-0.8, +0.8).
        eta_edges = lut.header.eta.edges()
        qa_eta_min = float(eta_edges[int(eta_indices[0])])
        qa_eta_max = float(eta_edges[int(eta_indices[-1]) + 1])

        qa_curves = read_qa_resolution(
            ROOT,
            args.qa_file,
            lut.header.pdg,
            qa_eta_min,
            qa_eta_max,
            args.pt_min,
            args.pt_max,
            args.fractional_resolution,
        )

        args.output_dir.mkdir(parents=True, exist_ok=True)
        eta_centers = lut.header.eta.centers()[eta_indices]
        print_summary(lut, eta_indices, qa_curves, args)

        # Keep returned ROOT objects referenced until both canvases are saved.
        efficiency_objects = plot_efficiency(
            ROOT,
            pt,
            efficiencies,
            eta_centers,
            lut.header,
            args,
        )
        has_resolution, resolution_objects = plot_resolution(
            ROOT,
            pt,
            resolutions,
            qa_curves,
            eta_centers,
            lut.header,
            args,
        )
        _keep_alive = (efficiency_objects, resolution_objects)

        # Both canvases are fully drawn and updated before any file is written.
        # Interactive runs pause here so the plots can be inspected first.
        if not args.background:
            wait_for_confirmation(ROOT)

        save_canvas(
            efficiency_objects.canvas,
            args.output_dir,
            "efficiency_vs_pt",
            args.pdf,
        )
        save_canvas(
            resolution_objects.canvas,
            args.output_dir,
            "pt_resolution_lut_vs_qa",
            args.pdf,
        )

        print(f"Efficiency: {args.output_dir / 'efficiency_vs_pt.png'}")
        print(f"Comparison: {args.output_dir / 'pt_resolution_lut_vs_qa.png'}")
        if not has_resolution:
            message = (
                "WARNING: the selected LUT entries contain no positive covm[14] "
                "values; the resolution figure is diagnostic only."
            )
            print(message, file=sys.stderr)
            if args.require_resolution:
                return 2
        return 0
    except (FileNotFoundError, RuntimeError, ValueError, struct.error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
