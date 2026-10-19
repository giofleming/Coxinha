#!/usr/bin/env python3
"""Compare ALICE 3 LUT payload, LUT fast tracking, and ACTS performance.

Three performance sources can be overlaid:

* LUT fast-tracking output from ``AnalysisResults.root``;
* the expected efficiency and pT resolution stored in ``lut_pdg*.dat``;
* the full ACTS tracking performance.

The binary LUT performance is evaluated with the same definitions used by
``plot_lut_performance.py``:

    efficiency = entry.eff (or entry.eff2)
    sigma(pT) / pT = sqrt(entry.covm[14]) * pT

The ACTS and AnalysisResults resolutions are extracted as an RMS or Gaussian
sigma from a residual histogram containing

    (pT_reco - pT_gen) / pT_gen

versus generated pT.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path
from typing import Iterable, Sequence


DEFAULT_LUT_FILE = "SimOutput/AnalysisResults.root"
DEFAULT_ACTS_FILE = (
    "../output/test_nEv1000_PID211_nMeasMin7_ckfChi2Meas45_"
    "ckfMeasPerSurf1/performance_merged_ambi_tracks.root"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare LUT and ACTS tracking efficiency and pT resolution.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("lut_file", nargs="?", default=DEFAULT_LUT_FILE,
                        type=Path, help="LUT AnalysisResults ROOT file")
    parser.add_argument("acts_file", nargs="?", default=DEFAULT_ACTS_FILE,
                        type=Path, help="ACTS track-finding performance ROOT file")
    parser.add_argument(
        "--acts-resolution-file",
        type=Path,
        default=None,
        help="Optional ACTS ROOT file containing a pT-residual histogram",
    )
    parser.add_argument("--pdg", type=int, default=211,
                        help="Particle PDG code used for the LUT objects")
    parser.add_argument("--lut-label", default="LUT fast tracking")
    parser.add_argument("--lut-payload-label", default="LUT payload")
    parser.add_argument("--acts-label", default="ACTS CKF")
    parser.add_argument("-o", "--output-dir", type=Path,
                        default=Path("performance_comparison"))
    parser.add_argument("--root-output", type=Path, default=None,
                        help="Output ROOT file; defaults to OUTPUT_DIR/comparison.root")

    parser.add_argument(
        "--lut-efficiency",
        default=None,
        help=(
            "Path of the LUT pT-efficiency TEfficiency. The path may cross "
            "a TCollection, e.g. alice3-efficiency/output/efficiency_pdg211"
        ),
    )
    parser.add_argument("--acts-efficiency", default="trackeff_vs_pT",
                        help="Path of the ACTS pT-efficiency TEfficiency")

    parser.add_argument("--lut-resolution", default=None,
                        help="Full path of the LUT relative-pT residual TH2")
    parser.add_argument(
        "--lut-dat-file",
        type=Path,
        default=None,
        help=(
            "Optional binary lut_pdg*.dat file. If provided, both LUT efficiency "
            "and pT resolution are computed with plot_lut_performance.py and added "
            "to the comparison plots."
        ),
    )
    parser.add_argument(
        "--lut-performance-script",
        type=Path,
        default=None,
        help=(
            "Path to plot_lut_performance.py. When omitted, common locations next "
            "to this script and in ../luts are searched."
        ),
    )
    parser.add_argument(
        "--lut-efficiency-field",
        choices=("eff", "eff2"),
        default="eff",
        help="Efficiency field read from each binary LUT entry",
    )
    parser.add_argument(
        "--lut-average-only",
        action="store_true",
        help="Show only the eta-averaged curves in the standalone LUT payload plots",
    )
    parser.add_argument(
        "--no-lut-payload-plots",
        action="store_true",
        help="Do not write the standalone plots produced from the binary LUT payload",
    )
    parser.add_argument(
        "--lut-eta-min",
        type=float,
        default=-0.8,
        help="Minimum LUT eta-bin centre used for dat-based resolution",
    )
    parser.add_argument(
        "--lut-eta-max",
        type=float,
        default=0.8,
        help="Maximum LUT eta-bin centre used for dat-based resolution",
    )
    parser.add_argument(
        "--lut-nch-bin",
        type=int,
        default=0,
        help="LUT nch-bin index used for dat-based resolution",
    )
    parser.add_argument(
        "--lut-radius-bin",
        type=int,
        default=0,
        help="LUT radius-bin index used for dat-based resolution",
    )
    parser.add_argument("--acts-resolution", default=None,
                        help="Full path of the ACTS relative-pT residual TH2")
    parser.add_argument(
        "--resolution-estimator",
        choices=("rms", "gaus"),
        default="rms",
        help="Estimator used for each residual distribution",
    )
    parser.add_argument("--min-entries", type=int, default=30,
                        help="Minimum entries required in one pT bin")
    parser.add_argument("--gaus-range", type=float, default=2.0,
                        help="Gaussian fit range in units of the initial RMS")
    parser.add_argument(
        "--fractional-resolution",
        action="store_true",
        help="Plot resolution as a fraction instead of percent",
    )

    parser.add_argument("--pt-min", type=float, default=None,
                        help="Lower displayed pT bound")
    parser.add_argument("--pt-max", type=float, default=None,
                        help="Upper displayed pT bound")
    parser.add_argument("--linear-x", action="store_true",
                        help="Use a linear pT axis")
    parser.add_argument("--pdf", action="store_true",
                        help="Also save PDF versions of the plots")
    parser.add_argument("--require-resolution", action="store_true",
                        help="Return an error if an ACTS resolution object is absent")
    return parser.parse_args()


def import_root():
    try:
        import ROOT  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "PyROOT is required. Run this script in an environment where ROOT is loaded."
        ) from exc
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    return ROOT


def open_root(ROOT, path: Path):
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {path}")
    return root_file


def get_object(root_file, path: str, expected: str | None = None):
    """Resolve a slash-separated path through TDirectory and TCollection objects."""
    current = root_file
    traversed = []
    for component in (part for part in path.split("/") if part):
        traversed.append(component)
        if current.InheritsFrom("TDirectory"):
            next_object = current.Get(component)
        elif current.InheritsFrom("TCollection"):
            next_object = current.FindObject(component)
        else:
            parent = "/".join(traversed[:-1]) or "<file>"
            raise KeyError(
                f"Cannot descend into '{parent}': {current.ClassName()} is neither "
                "a TDirectory nor a TCollection"
            )
        if not next_object:
            raise KeyError(f"Object not found: {'/'.join(traversed)}")
        current = next_object

    if expected and not current.InheritsFrom(expected):
        raise TypeError(
            f"Object '{path}' is a {current.ClassName()}, expected {expected}"
        )
    return current


def walk_objects(container, prefix: str = "") -> Iterable[tuple[str, str]]:
    """Yield (full path, class name) through directories and collections."""
    if container.InheritsFrom("TDirectory"):
        keys = container.GetListOfKeys()
        if not keys:
            return
        children = ((key.GetName(), container.Get(key.GetName())) for key in keys)
    elif container.InheritsFrom("TCollection"):
        children = ((obj.GetName(), obj) for obj in container)
    else:
        return

    for name, obj in children:
        if not obj:
            continue
        path = f"{prefix}/{name}" if prefix else name
        if obj.InheritsFrom("TDirectory") or obj.InheritsFrom("TCollection"):
            yield from walk_objects(obj, path)
        else:
            yield path, obj.ClassName()


def pdg_tag(pdg: int) -> str:
    return f"m{abs(pdg)}" if pdg < 0 else str(pdg)


def default_lut_efficiency_path(args: argparse.Namespace) -> str:
    tag = pdg_tag(args.pdg)
    return args.lut_efficiency or f"alice3-efficiency/output/efficiency_pdg{tag}"


def default_lut_resolution_path(args: argparse.Namespace) -> str:
    tag = pdg_tag(args.pdg)
    return args.lut_resolution or (
        f"alice3-tracking-performance/{tag}/ptResolutionVsPt_{tag}"
    )


def efficiency_graph(obj, name: str):
    if obj.InheritsFrom("TEfficiency"):
        graph = obj.CreateGraph()
    elif obj.InheritsFrom("TGraph"):
        graph = obj.Clone()
    else:
        raise TypeError(
            f"Efficiency object '{obj.GetName()}' is a {obj.ClassName()}, "
            "expected TEfficiency or TGraph"
        )
    graph.SetName(name)
    return graph


def graph_x_range(graph) -> tuple[float, float]:
    values = []
    for i in range(graph.GetN()):
        x = float(graph.GetX()[i])
        ex_low = float(graph.GetErrorXlow(i)) if graph.InheritsFrom("TGraphAsymmErrors") else 0.0
        ex_high = float(graph.GetErrorXhigh(i)) if graph.InheritsFrom("TGraphAsymmErrors") else 0.0
        values.extend((x - ex_low, x + ex_high))
    positive = [value for value in values if value > 0.0 and math.isfinite(value)]
    if not positive:
        return 0.1, 10.0
    return min(positive), max(positive)


def common_x_range(graphs, requested_min: float | None,
                   requested_max: float | None) -> tuple[float, float]:
    ranges = [graph_x_range(graph) for graph in graphs]
    x_min = requested_min if requested_min is not None else max(r[0] for r in ranges)
    x_max = requested_max if requested_max is not None else min(r[1] for r in ranges)
    if x_max <= x_min:
        # Fall back to the union when the automatic overlap is empty.
        x_min = requested_min if requested_min is not None else min(r[0] for r in ranges)
        x_max = requested_max if requested_max is not None else max(r[1] for r in ranges)
    if x_max <= x_min:
        raise RuntimeError(f"Invalid pT range: [{x_min}, {x_max}]")
    return x_min, x_max


def style_graph(ROOT, graph, index: int):
    colors = (
        ROOT.kAzure + 2,
        ROOT.kGreen + 2,
        ROOT.kOrange + 7,
        ROOT.kMagenta + 1,
        ROOT.kRed + 1,
    )
    markers = (
        ROOT.kFullCircle,
        ROOT.kOpenDiamond,
        ROOT.kFullSquare,
        ROOT.kOpenTriangleUp,
        ROOT.kFullTriangleDown,
    )
    graph.SetLineColor(colors[index % len(colors)])
    graph.SetMarkerColor(colors[index % len(colors)])
    graph.SetMarkerStyle(markers[index % len(markers)])
    graph.SetMarkerSize(0.9)
    graph.SetLineWidth(2)


def draw_label(ROOT, text: str):
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextFont(42)
    label.SetTextSize(0.040)
    label.DrawLatex(0.16, 0.92, text)
    return label


def save_canvas(canvas, output_dir: Path, stem: str, save_pdf: bool):
    canvas.SaveAs(str(output_dir / f"{stem}.png"))
    if save_pdf:
        canvas.SaveAs(str(output_dir / f"{stem}.pdf"))


def draw_efficiency(
    ROOT,
    graph_series: Sequence[tuple[object, str]],
    args,
    output_file,
):
    if not graph_series:
        raise RuntimeError("No efficiency curves were supplied")

    for index, (graph, _) in enumerate(graph_series):
        style_graph(ROOT, graph, index)

    x_min, x_max = common_x_range(
        tuple(graph for graph, _ in graph_series), args.pt_min, args.pt_max
    )

    canvas = ROOT.TCanvas("c_efficiency", "Tracking efficiency", 900, 750)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)
    if not args.linear_x:
        if x_min <= 0.0:
            raise RuntimeError("A logarithmic pT axis requires --pt-min > 0")
        canvas.SetLogx()

    frame = ROOT.TH1F("efficiency_frame", "", 1, x_min, x_max)
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.08)
    frame.GetXaxis().SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})")
    frame.GetYaxis().SetTitle("Tracking efficiency")
    frame.GetXaxis().SetTitleOffset(1.10)
    frame.GetYaxis().SetTitleOffset(1.30)
    frame.Draw("AXIS")

    for graph, _ in graph_series:
        graph.Draw("PZ SAME")

    legend_height = 0.055 * len(graph_series) + 0.04
    legend = ROOT.TLegend(0.53, 0.18, 0.91, 0.18 + legend_height)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for graph, label_text in graph_series:
        legend.AddEntry(graph, label_text, "lp")
    legend.Draw()
    label = draw_label(ROOT, f"Tracking performance, PDG {args.pdg}")

    canvas.RedrawAxis()
    save_canvas(canvas, args.output_dir, "efficiency_vs_pt", args.pdf)
    output_file.cd()
    frame.Write()
    for graph, _ in graph_series:
        graph.Write()
    canvas.Write()
    return canvas, frame, legend, label


def resolution_from_th2(ROOT, hist, name: str, args):
    graph = ROOT.TGraphErrors()
    graph.SetName(name)
    graph.SetTitle("")
    point = 0
    scale = 1.0 if args.fractional_resolution else 100.0

    for x_bin in range(1, hist.GetNbinsX() + 1):
        projection = hist.ProjectionY(
            f"{name}_projection_{x_bin}", x_bin, x_bin, "e"
        )
        projection.SetDirectory(0)
        if projection.GetEntries() < args.min_entries:
            continue

        sigma = float(projection.GetRMS())
        sigma_error = float(projection.GetRMSError())
        if args.resolution_estimator == "gaus":
            mean = float(projection.GetMean())
            rms = float(projection.GetRMS())
            if not math.isfinite(rms) or rms <= 0.0:
                continue
            fit_min = mean - args.gaus_range * rms
            fit_max = mean + args.gaus_range * rms
            fit = ROOT.TF1(f"{name}_fit_{x_bin}", "gaus", fit_min, fit_max)
            fit.SetParameters(projection.GetMaximum(), mean, rms)
            fit_result = projection.Fit(fit, "QNR0S")
            if int(fit_result) != 0:
                continue
            sigma = abs(float(fit.GetParameter(2)))
            sigma_error = float(fit.GetParError(2))

        if not math.isfinite(sigma) or sigma <= 0.0:
            continue
        x = float(hist.GetXaxis().GetBinCenter(x_bin))
        x_error = 0.5 * float(hist.GetXaxis().GetBinWidth(x_bin))
        graph.SetPoint(point, x, scale * sigma)
        graph.SetPointError(point, x_error, scale * sigma_error)
        point += 1

    if graph.GetN() == 0:
        raise RuntimeError(
            f"No valid resolution points could be extracted from '{hist.GetName()}'. "
            f"Try reducing --min-entries or changing --resolution-estimator."
        )
    return graph


def scale_graph_y(graph, factor: float):
    for i in range(graph.GetN()):
        x = float(graph.GetX()[i])
        y = factor * float(graph.GetY()[i])
        graph.SetPoint(i, x, y)
        if graph.InheritsFrom("TGraphAsymmErrors"):
            graph.SetPointError(
                i,
                graph.GetErrorXlow(i),
                graph.GetErrorXhigh(i),
                factor * graph.GetErrorYlow(i),
                factor * graph.GetErrorYhigh(i),
            )
        elif graph.InheritsFrom("TGraphErrors"):
            graph.SetPointError(
                i, graph.GetErrorX(i), factor * graph.GetErrorY(i)
            )


def resolution_graph(ROOT, obj, name: str, args):
    if obj.InheritsFrom("TH2"):
        return resolution_from_th2(ROOT, obj, name, args)
    if obj.InheritsFrom("TGraph"):
        graph = obj.Clone(name)
    elif obj.InheritsFrom("TH1"):
        graph = ROOT.TGraphErrors(obj)
        graph.SetName(name)
    else:
        raise TypeError(
            f"Resolution object '{obj.GetName()}' is a {obj.ClassName()}, "
            "expected TH2 residual distribution, TH1, or TGraph"
        )

    if not args.fractional_resolution:
        scale_graph_y(graph, 100.0)
    return graph


def _import_numpy():
    try:
        import numpy as np  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "NumPy is required for --lut-dat-file performance extraction"
        ) from exc
    return np


def _import_lut_plot_module(args: argparse.Namespace):
    script_dir = Path(__file__).resolve().parent
    candidates: list[Path] = []
    if args.lut_performance_script is not None:
        candidates.append(args.lut_performance_script)
    candidates.extend(
        (
            script_dir / "plot_lut_performance.py",
            script_dir / "plot_lut_performance(1).py",
            script_dir / "luts" / "plot_lut_performance.py",
            script_dir.parent / "luts" / "plot_lut_performance.py",
        )
    )

    module_path = next((path.resolve() for path in candidates if path.is_file()), None)
    if module_path is None:
        searched = "\n  ".join(str(path) for path in candidates)
        raise RuntimeError(
            "Cannot locate plot_lut_performance.py. Use --lut-performance-script. "
            "Searched:\n  " + searched
        )

    spec = importlib.util.spec_from_file_location("plot_lut_performance", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import LUT plotting module from: {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, module_path


def graph_from_numpy(ROOT, name: str, x_values, y_values):
    graph = ROOT.TGraphErrors()
    graph.SetName(name)
    graph.SetTitle("")
    point = 0
    for x, y in zip(x_values, y_values):
        if not math.isfinite(float(x)) or not math.isfinite(float(y)):
            continue
        graph.SetPoint(point, float(x), float(y))
        graph.SetPointError(point, 0.0, 0.0)
        point += 1
    if graph.GetN() == 0:
        raise RuntimeError(f"No finite points were produced for graph '{name}'")
    return graph


def lut_performance_from_dat(ROOT, lut_dat_file: Path, args):
    """Read binary LUT efficiency/resolution using plot_lut_performance.py."""
    np = _import_numpy()
    lut_plot, module_path = _import_lut_plot_module(args)

    lut = lut_plot.read_lut(lut_dat_file)
    if abs(int(lut.header.pdg)) != abs(int(args.pdg)):
        raise ValueError(
            f"Binary LUT PDG {lut.header.pdg} does not match requested PDG {args.pdg}"
        )

    lut_plot.validate_bin_index(
        args.lut_nch_bin, lut.header.nch.nbins, "lut-nch-bin"
    )
    lut_plot.validate_bin_index(
        args.lut_radius_bin, lut.header.radius.nbins, "lut-radius-bin"
    )
    eta_indices = lut_plot.selected_eta_indices(
        lut.header, args.lut_eta_min, args.lut_eta_max
    )

    pt_all = lut.header.pt.centers()
    pt_mask = lut_plot.finite_pt_mask(pt_all, args.pt_min, args.pt_max)
    pt = pt_all[pt_mask]

    selected = lut.entries[
        args.lut_nch_bin,
        args.lut_radius_bin,
        eta_indices,
        :,
    ]

    efficiencies = np.asarray(
        selected[args.lut_efficiency_field], dtype=float
    )[:, pt_mask]

    q_over_pt_variance = np.asarray(
        selected["covm"][..., 14], dtype=float
    )[:, pt_mask]
    q_over_pt_variance = np.where(
        q_over_pt_variance >= 0.0, q_over_pt_variance, np.nan
    )
    resolutions = np.sqrt(q_over_pt_variance) * pt[np.newaxis, :]
    if not args.fractional_resolution:
        resolutions *= 100.0

    eta_centers = lut.header.eta.centers()[eta_indices]

    # Use the exact finite-average convention from plot_lut_performance.py.
    average_efficiency = lut_plot.finite_average(efficiencies)
    positive_resolutions = np.where(resolutions > 0.0, resolutions, np.nan)
    average_resolution = lut_plot.finite_average(positive_resolutions)

    efficiency_graph_payload = graph_from_numpy(
        ROOT,
        "lut_payload_efficiency_vs_pt",
        pt,
        average_efficiency,
    )
    resolution_graph_payload = graph_from_numpy(
        ROOT,
        "lut_payload_pt_resolution",
        pt,
        average_resolution,
    )

    standalone_objects: list[object] = []
    if not args.no_lut_payload_plots:
        plot_args = argparse.Namespace(**vars(args))
        plot_args.output_dir = args.output_dir / "lut_payload"
        plot_args.output_dir.mkdir(parents=True, exist_ok=True)
        plot_args.average_only = args.lut_average_only
        plot_args.efficiency_field = args.lut_efficiency_field

        efficiency_objects = lut_plot.plot_efficiency(
            ROOT,
            pt,
            efficiencies,
            eta_centers,
            lut.header,
            plot_args,
        )
        _, resolution_objects = lut_plot.plot_resolution(
            ROOT,
            pt,
            resolutions,
            eta_centers,
            lut.header,
            plot_args,
        )
        standalone_objects.extend((efficiency_objects, resolution_objects))

    eta_text = f"{eta_centers.min():g} <= eta <= {eta_centers.max():g}"
    print(f"LUT payload module:     {module_path}")
    print(f"LUT payload selection:  {eta_text}")
    print(f"LUT payload efficiency: entry.{args.lut_efficiency_field}")
    print("LUT payload resolution: sqrt(covm[14]) * pT")

    return (
        efficiency_graph_payload,
        resolution_graph_payload,
        standalone_objects,
    )


def find_acts_resolution(root_file) -> list[tuple[str, str]]:
    """Find plausible *relative pT residual versus pT* objects only."""
    candidates = []
    for path, class_name in walk_objects(root_file):
        leaf = path.rsplit("/", 1)[-1].lower()
        has_pt = "pt" in leaf
        has_resolution = (
            "resolution" in leaf or "residual" in leaf or leaf.startswith("res_")
        )
        excluded = any(token in leaf for token in (
            "dca", "eta", "phi", "theta", "qop", "qoverp", "invpt",
            "pull", "cov", "eff", "fake", "dup",
        ))
        if has_pt and has_resolution and not excluded:
            candidates.append((path, class_name))
    return candidates


def select_acts_resolution(root_file, requested: str | None):
    if requested:
        return requested
    candidates = find_acts_resolution(root_file)
    if len(candidates) == 1:
        return candidates[0][0]
    if not candidates:
        return None
    names = "\n  ".join(f"{path} [{class_name}]" for path, class_name in candidates)
    raise RuntimeError(
        "More than one possible ACTS pT-resolution object was found. "
        "Choose one with --acts-resolution:\n  " + names
    )


def draw_resolution(
    ROOT,
    graph_series: Sequence[tuple[object, str]],
    args,
    output_file,
):
    if not graph_series:
        raise RuntimeError("No resolution curves were supplied")

    for index, (graph, _) in enumerate(graph_series):
        style_graph(ROOT, graph, index)

    x_min, x_max = common_x_range(
        tuple(graph for graph, _ in graph_series), args.pt_min, args.pt_max
    )

    y_values = []
    for graph, _ in graph_series:
        for i in range(graph.GetN()):
            value = float(graph.GetY()[i])
            if math.isfinite(value) and value >= 0.0:
                y_values.append(value)
    y_max = max(1.0, 1.25 * max(y_values)) if y_values else 1.0

    canvas = ROOT.TCanvas("c_pt_resolution", "pT resolution", 900, 750)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)
    if not args.linear_x:
        if x_min <= 0.0:
            raise RuntimeError("A logarithmic pT axis requires --pt-min > 0")
        canvas.SetLogx()

    frame = canvas.DrawFrame(x_min, 0.0, x_max, y_max)
    frame.GetXaxis().SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})")
    unit = " (%)" if not args.fractional_resolution else ""
    frame.GetYaxis().SetTitle("#sigma(#it{p}_{T})/#it{p}_{T}" + unit)
    frame.GetXaxis().SetTitleOffset(1.10)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.Draw("AXIS")

    for graph, _ in graph_series:
        graph.Draw("PZ SAME")

    legend_height = 0.055 * len(graph_series) + 0.04
    legend = ROOT.TLegend(0.50, 0.86 - legend_height, 0.91, 0.86)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for graph, label_text in graph_series:
        legend.AddEntry(graph, label_text, "lp")
    legend.Draw()

    label = draw_label(ROOT, f"Relative #it{{p}}_{{T}} resolution, PDG {args.pdg}")

    canvas.RedrawAxis()
    save_canvas(canvas, args.output_dir, "pt_resolution_vs_pt", args.pdf)
    output_file.cd()
    frame.Write()
    for graph, _ in graph_series:
        graph.Write()
    canvas.Write()
    return canvas, frame, legend, label


def draw_lut_resolution_only(ROOT, lut_graph, args, output_file):
    style_graph(ROOT, lut_graph, 0)
    x_min, x_max = common_x_range((lut_graph,), args.pt_min, args.pt_max)

    y_values = []
    for i in range(lut_graph.GetN()):
        value = float(lut_graph.GetY()[i])
        if math.isfinite(value) and value >= 0.0:
            y_values.append(value)
    y_max = max(1.0, 1.25 * max(y_values)) if y_values else 1.0

    canvas = ROOT.TCanvas(
        "c_pt_resolution_analysisresults",
        "LUT pT resolution from AnalysisResults",
        900,
        750,
    )
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)
    if not args.linear_x:
        if x_min <= 0.0:
            raise RuntimeError("A logarithmic pT axis requires --pt-min > 0")
        canvas.SetLogx()

    frame = ROOT.TH1F("pt_resolution_analysisresults_frame", "", 1, x_min, x_max)
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(y_max)
    frame.GetXaxis().SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})")
    unit = " (%)" if not args.fractional_resolution else ""
    frame.GetYaxis().SetTitle(
        "#sigma((#it{p}_{T}^{reco}-#it{p}_{T}^{gen})/#it{p}_{T}^{gen})" + unit
    )
    frame.GetXaxis().SetTitleOffset(1.10)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.Draw("AXIS")

    lut_graph.Draw("PZ SAME")
    legend = ROOT.TLegend(0.56, 0.72, 0.91, 0.83)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(lut_graph, "LUT from AnalysisResults", "lp")
    legend.Draw()
    estimator = "RMS" if args.resolution_estimator == "rms" else "Gaussian #sigma"
    label = draw_label(ROOT, f"LUT pT resolution from AnalysisResults ({estimator}), PDG {args.pdg}")

    canvas.RedrawAxis()
    save_canvas(canvas, args.output_dir, "pt_resolution_vs_pt_analysisresults", args.pdf)
    output_file.cd()
    frame.Write()
    lut_graph.Write()
    canvas.Write()
    return canvas, frame, legend, label


def print_available_objects(root_file, heading: str):
    print(heading, file=sys.stderr)
    for path, class_name in walk_objects(root_file):
        print(f"  {path} [{class_name}]", file=sys.stderr)


def main() -> int:
    args = parse_args()
    ROOT = import_root()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root_output = args.root_output or args.output_dir / "comparison.root"
    root_output.parent.mkdir(parents=True, exist_ok=True)

    lut_file = open_root(ROOT, args.lut_file)
    acts_file = open_root(ROOT, args.acts_file)
    acts_resolution_file = acts_file
    owns_resolution_file = False
    if args.acts_resolution_file is not None:
        acts_resolution_file = open_root(ROOT, args.acts_resolution_file)
        owns_resolution_file = True

    output_file = ROOT.TFile.Open(str(root_output), "RECREATE")
    if not output_file or output_file.IsZombie():
        raise RuntimeError(f"Cannot create output ROOT file: {root_output}")

    keep_alive: list[object] = []
    try:
        efficiency_series: list[tuple[object, str]] = []
        resolution_series: list[tuple[object, str]] = []

        # Fast-tracking performance measured after applying the LUT.
        lut_efficiency_path = default_lut_efficiency_path(args)
        lut_efficiency = get_object(lut_file, lut_efficiency_path, "TEfficiency")
        lut_eff_graph = efficiency_graph(
            lut_efficiency, "lut_fast_tracking_efficiency_vs_pt"
        )
        efficiency_series.append((lut_eff_graph, args.lut_label))
        keep_alive.append(lut_eff_graph)
        print(f"LUT fast-tracking efficiency: {lut_efficiency_path}")

        lut_resolution_path = default_lut_resolution_path(args)
        lut_resolution = get_object(lut_file, lut_resolution_path)
        lut_resolution_graph_root = resolution_graph(
            ROOT,
            lut_resolution,
            "lut_fast_tracking_pt_resolution",
            args,
        )
        resolution_series.append((lut_resolution_graph_root, args.lut_label))
        keep_alive.append(lut_resolution_graph_root)
        print(f"LUT fast-tracking resolution: {lut_resolution_path}")

        # Expected performance encoded directly in the binary LUT payload.
        if args.lut_dat_file is not None:
            (
                payload_efficiency_graph,
                payload_resolution_graph,
                payload_plot_objects,
            ) = lut_performance_from_dat(ROOT, args.lut_dat_file, args)
            efficiency_series.append(
                (payload_efficiency_graph, args.lut_payload_label)
            )
            resolution_series.append(
                (payload_resolution_graph, args.lut_payload_label)
            )
            keep_alive.extend(
                (
                    payload_efficiency_graph,
                    payload_resolution_graph,
                    *payload_plot_objects,
                )
            )
            print(f"LUT payload: {args.lut_dat_file}")

        # Full ACTS efficiency.
        acts_efficiency = get_object(
            acts_file, args.acts_efficiency, "TEfficiency"
        )
        acts_eff_graph = efficiency_graph(
            acts_efficiency, "acts_efficiency_vs_pt"
        )
        efficiency_series.append((acts_eff_graph, args.acts_label))
        keep_alive.append(acts_eff_graph)
        print(f"ACTS efficiency: {args.acts_efficiency}")

        keep_alive.extend(
            draw_efficiency(ROOT, efficiency_series, args, output_file)
        )

        # Keep the original standalone AnalysisResults resolution figure.
        keep_alive.extend(
            draw_lut_resolution_only(
                ROOT,
                lut_resolution_graph_root,
                args,
                output_file,
            )
        )

        acts_resolution_path = select_acts_resolution(
            acts_resolution_file, args.acts_resolution
        )
        if acts_resolution_path is None:
            message = (
                f"No relative-pT resolution object was found in "
                f"'{acts_resolution_file.GetName()}'. The available LUT resolution "
                "curves are still plotted, but ACTS resolution requires a fitter-"
                "performance file or an explicit --acts-resolution object."
            )
            print(f"WARNING: {message}", file=sys.stderr)
            print_available_objects(
                acts_resolution_file,
                "Objects available in the ACTS file:",
            )
            if args.require_resolution:
                raise RuntimeError(message)
        else:
            acts_resolution = get_object(
                acts_resolution_file, acts_resolution_path
            )
            acts_resolution_graph = resolution_graph(
                ROOT, acts_resolution, "acts_pt_resolution", args
            )
            resolution_series.append((acts_resolution_graph, args.acts_label))
            keep_alive.append(acts_resolution_graph)
            print(f"ACTS resolution: {acts_resolution_path}")

        if len(resolution_series) >= 2:
            keep_alive.extend(
                draw_resolution(ROOT, resolution_series, args, output_file)
            )
        else:
            print(
                "WARNING: fewer than two resolution curves are available; "
                "no comparison resolution plot was produced.",
                file=sys.stderr,
            )

        output_file.Write()
        print(f"Plots written to: {args.output_dir}")
        if args.lut_dat_file is not None and not args.no_lut_payload_plots:
            print(f"LUT payload plots: {args.output_dir / 'lut_payload'}")
        print(f"ROOT output:      {root_output}")
        return 0
    finally:
        output_file.Close()
        if owns_resolution_file:
            acts_resolution_file.Close()
        acts_file.Close()
        lut_file.Close()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, TypeError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
