"""Figures: per-cell profile, segmentation QC overlay, group boxplot."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import find_boundaries

COLORS = ["#7f8c8d", "#e74c3c", "#e67e22", "#f1c40f", "#16a085", "#3498db", "#9b59b6"]


def plot_cell(path, plane, outline_yx, dist_um, profile, title, pixel_um, pad_um=6):
    """Left: cell crop with the traced outline. Right: intensity along the outline (like an ImageJ plot profile)."""
    pad = int(pad_um / pixel_um)
    y0, x0 = np.maximum(outline_yx.min(0).astype(int) - pad, 0)
    y1, x1 = outline_yx.max(0).astype(int) + pad
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw=dict(width_ratios=[1, 1.4]))
    a.imshow(plane[y0:y1, x0:x1], cmap="gray", vmax=np.percentile(plane[y0:y1, x0:x1], 99.8))
    a.scatter(outline_yx[:, 1] - x0, outline_yx[:, 0] - y0, c=dist_um, cmap="autumn_r", s=6)
    a.set_title(title, fontsize=9)
    a.axis("off")
    b.plot(dist_um, profile, color="k", lw=1)
    b.set(xlabel="Distance (microns)", ylabel="Gray Value", xlim=(0, dist_um[-1]))
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_qc(path, img, labels, title):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(img, cmap="gray", vmax=np.percentile(img, 99.7))
    edge = np.zeros(labels.shape + (4,))
    edge[find_boundaries(labels, mode="inner")] = [1, 0, 0, 1]
    ax.imshow(edge)
    for l in range(1, labels.max() + 1):
        y, x = np.argwhere(labels == l).mean(0)
        ax.text(x, y, str(l), color="cyan", ha="center", va="center", fontsize=11)
    ax.set_title(title, fontsize=10)
    ax.axis("off")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=90)
    plt.close(fig)


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def plot_groups(path, df, stat_table, kruskal_p, control_label):
    """Box + jittered cells per time point, median trend line, significance vs control."""
    order = list(dict.fromkeys(df.sort_values("Time_min").Condition))
    p_vs = stat_table.set_index("Condition").p_vs_control
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(10, 6.6))
    medians, ymax = [], df.CV.max()
    for i, cond in enumerate(order):
        v = df.loc[df.Condition == cond, "CV"].dropna().values
        c = COLORS[i % len(COLORS)]
        ax.boxplot(v, positions=[i], widths=0.5, patch_artist=True, showfliers=False,
                   boxprops=dict(facecolor=c, alpha=0.5, edgecolor="gray"),
                   medianprops=dict(color="k", lw=2), whiskerprops=dict(color="gray"),
                   capprops=dict(color="gray"))
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=45, color=c, edgecolor="white", zorder=3)
        medians.append(np.median(v))
        ax.text(i, 0.02, f"n={len(v)}", transform=ax.get_xaxis_transform(), ha="center", color="gray", fontsize=9)
        s = stars(p_vs.get(cond, 1.0)) if cond != control_label else ""
        if s:
            ax.text(i, v.max() + 0.03 * ymax, s, ha="center", fontsize=14)
    ax.plot(range(len(order)), medians, color="#8a94a0", lw=1.8, zorder=2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_ylabel("CV (membrane intensity)", fontsize=12)
    ax.set_ylim(top=ymax * 1.15)
    ax.grid(axis="y", alpha=0.25)
    ax.set_title("Receptor (mCherry) clustering after ligand stimulation\n"
                 "higher CV = more clustered   |   lower CV = more uniform", fontsize=12)
    foot = (f"Kruskal-Wallis p={kruskal_p:.3g}   |   significance vs {control_label} (Mann-Whitney): "
            f"*p<0.05  **p<0.01  ***p<0.001   |   n={len(df)} cells\n"
            "Cells are pseudo-replicates (not independent): p-values are overestimated (anti-conservative).")
    fig.text(0.01, 0.005, foot, fontsize=8, color="gray", va="bottom")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)
