"""End-to-end analysis: CZI images -> per-cell figures, cell table, group statistics and plot."""
import argparse
from pathlib import Path

import pandas as pd
from scipy import stats
from skimage.measure import regionprops

from . import clustering, plotting, segmentation
from .io import load_czi

COLUMNS = ["Image", "Cell", "Membrane_Mean", "SD", "Background", "Signal", "CV"]


def analyse_image(czi_path, out_dir, receptor_ch, nuclear_ch, band_um):
    name = Path(czi_path).stem
    stack, nuc, px = load_czi(czi_path, receptor_ch, nuclear_ch)
    labels, _ = segmentation.segment_cells(stack, px, nuc, skip_first=2)
    plotting.plot_qc(out_dir / "segmentation_qc" / f"{name}.png", segmentation.projection(stack, 2), labels,
                     f"{name}: {labels.max()} cells (border cells excluded)")
    bg_cache, rows = {}, []
    for r in regionprops(labels):
        mask = labels == r.label
        z = clustering.best_plane(stack, mask, px, nuc)
        plane = clustering.equatorial_plane(stack, z)
        if z not in bg_cache:
            bg_cache[z] = clustering.estimate_background(plane, labels, px)
        dist, prof, outline = clustering.membrane_profile(plane, mask, px, band_um)
        m = clustering.cell_metrics(prof, bg_cache[z])
        plotting.plot_cell(out_dir / "per_cell" / f"{name}_cell{r.label:02d}.png", plane, outline, dist, prof,
                           f"{name}  cell {r.label}  z={z}  CV={m['CV']:.2f}", px)
        rows.append(dict(Image=name, Cell=r.label, **m, Z_plane=z, Perimeter_um=float(dist[-1]),
                         Area_um2=r.area * px**2))
    return rows


def group_stats(df, control):
    ctrl = df.loc[df.Condition == control, "CV"].dropna()
    kw_p = stats.kruskal(*[g.CV.dropna() for _, g in df.groupby("Condition")]).pvalue
    rows = [dict(Condition=c, n=int(g.CV.notna().sum()), median_CV=g.CV.median(),
                 p_vs_control=stats.mannwhitneyu(g.CV.dropna(), ctrl).pvalue if c != control else float("nan"))
            for c, g in df.groupby("Condition")]
    return pd.DataFrame(rows), kw_p


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="data/raw", help="folder with .czi files")
    ap.add_argument("--samples", default="metadata/samples.csv", help="image -> condition / time table")
    ap.add_argument("--out", default="results")
    ap.add_argument("--receptor-channel", type=int, default=0, help="0-based channel of mCherry (R-PE)")
    ap.add_argument("--nuclear-channel", type=int, default=2,
                    help="0-based nuclear (DAPI) channel used to split touching cells; -1 = none")
    ap.add_argument("--band-um", type=float, default=0.8, help="half-width of the membrane band")
    a = ap.parse_args(argv)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(a.samples)
    nuc = None if a.nuclear_channel < 0 else a.nuclear_channel
    rows = []
    for s in samples.itertuples():
        print(f"{s.image} ({s.condition})", flush=True)
        for r in analyse_image(Path(a.data_dir) / f"{s.image}.czi", out, a.receptor_channel, nuc, a.band_um):
            rows.append(dict(r, Condition=s.condition, Time_min=s.time_min))
    df = pd.DataFrame(rows).replace([float("inf")], float("nan")).dropna(subset=["CV"])
    df = df[COLUMNS + [c for c in df.columns if c not in COLUMNS]]
    df.to_csv(out / "cell_table.csv", index=False)
    df.to_excel(out / "cell_table.xlsx", index=False)

    control = samples.sort_values("time_min").condition.iloc[0]
    table, kw_p = group_stats(df, control)
    table.to_csv(out / "group_stats.csv", index=False)
    print(table.to_string(index=False), f"\nKruskal-Wallis p = {kw_p:.3g}")
    plotting.plot_groups(out / "cv_boxplot.png", df, table, kw_p, control)
