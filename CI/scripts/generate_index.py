"""
Generate a static index.html gallery page listing all PNG plots
in a given output directory. Designed to run after compare_monitoring.py
has written its plots, as a standalone step in CI.

Usage:
    python generate_index.py --output-dir comparison_output --title "ACTS Alice3 Monitoring"
"""

import argparse
import os
from datetime import datetime, timezone
import subprocess


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
    background: #1e1e1e;
    color: #ddd;
    margin: 0;
    padding: 2rem;
  }}
  h1 {{
    font-size: 1.5rem;
    margin-bottom: 0.25rem;
  }}
  .meta {{
    color: #888;
    font-size: 0.85rem;
    margin-bottom: 2rem;
  }}
  .gallery {{
    display: flex;
    flex-direction: column;
    gap: 2rem;
  }}
  .plot-card {{
    background: #2a2a2a;
    border-radius: 8px;
    padding: 1rem;
    border: 1px solid #3a3a3a;
  }}
  .plot-card h2 {{
    font-size: 1rem;
    font-weight: 500;
    margin: 0 0 0.75rem 0;
    color: #eee;
  }}
  .plot-card img {{
    max-width: 100%;
    height: auto;
    border-radius: 4px;
    background: #fff;
  }}
  .empty {{
    color: #888;
    font-style: italic;
  }}
</style>
</head>
<body>
  <h1>{title}</h1>
  <div class="meta">Generated {timestamp}{commit_info}</div>
  <div class="gallery">
{plot_cards}
  </div>
</body>
</html>
"""

PLOT_CARD_TEMPLATE = """    <div class="plot-card">
      <h2>{name}</h2>
      <img src="{filename}" alt="{name}" loading="lazy">
    </div>"""


def humanize(filename: str) -> str:
    """Turn a filename like 'residual_x_overlay.png' into 'Residual X Overlay'."""
    name = os.path.splitext(filename)[0]
    return name.replace("_", " ").replace("-", " ").title()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--convert-to-png", help="Convert the pdfs to pngs", action="store_true", default=False)
    parser.add_argument("--output-dir", required=True,
                        help="Directory containing PNG plots; index.html is written here")
    parser.add_argument("--title", default="Monitoring Comparison", help="Page title")
    args = parser.parse_args()

    out_dir = args.output_dir
    if not os.path.isdir(out_dir):
        raise SystemExit(f"Output directory does not exist: {out_dir}")

    if args.convert_to_png:

        pdfs = sorted(f for f in os.listdir(out_dir) if f.lower().endswith(".pdf"))
        for pdf in pdfs:
            pdf_path = os.path.join(out_dir, pdf)
            base = os.path.splitext(pdf_path)[0]
            subprocess.run(
                ["pdftoppm", "-png", "-r", "150", "-singlefile", pdf_path, base],
                check=True,
            )

    pngs = sorted(f for f in os.listdir(out_dir) if f.lower().endswith(".png"))

    commit_sha = os.environ.get("CI_COMMIT_SHORT_SHA")
    mr_iid = os.environ.get("CI_MERGE_REQUEST_IID")
    commit_info = ""
    if commit_sha:
        commit_info = f" &middot; commit {commit_sha}"
    if mr_iid:
        commit_info += f" &middot; MR !{mr_iid}"

    if pngs:
        plot_cards = "\n".join(
            PLOT_CARD_TEMPLATE.format(name=humanize(f), filename=f) for f in pngs
        )
    else:
        plot_cards = '    <div class="empty">No plots found in this directory.</div>'

    html = HTML_TEMPLATE.format(
        title=args.title,
        timestamp=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        commit_info=commit_info,
        plot_cards=plot_cards,
    )

    index_path = os.path.join(out_dir, "index.html")
    with open(index_path, "w") as f:
        f.write(html)

    print(f"Wrote {index_path} with {len(pngs)} plot(s).")


if __name__ == "__main__":
    main()
