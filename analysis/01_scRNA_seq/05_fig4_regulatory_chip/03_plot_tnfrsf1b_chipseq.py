#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pyBigWig
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = Path(".")
WINDOW_BED = BASE / "ref" / "Tnfrsf1b.mm10.window.bed"
GENE_BED = BASE / "ref" / "Tnfrsf1b.mm10.bed"

TRACKS = [
    ("BATF Treg ChIP",  BASE / "bigwig" / "BATF_Treg.bw"),
    ("BATF Tconv ChIP", BASE / "bigwig" / "BATF_Tconv.bw"),
    ("Input Treg",      BASE / "bigwig" / "Input_Treg.bw"),
    ("Input Tconv",     BASE / "bigwig" / "Input_Tconv.bw"),
]

PEAKS = [
    ("BATF Treg peaks",  BASE / "peaks" / "BATF_Treg_peaks.narrowPeak"),
    ("BATF Tconv peaks", BASE / "peaks" / "BATF_Tconv_peaks.narrowPeak"),
]


def read_bed(path):
    fields = path.read_text().strip().splitlines()[0].split("\t")
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    name = fields[3] if len(fields) > 3 else ""
    strand = fields[4] if len(fields) > 4 else "."
    return chrom, start, end, name, strand


def same_chrom(a, b):
    return a == b or a.replace("chr", "") == b.replace("chr", "")


def choose_chrom_for_bw(bw, chrom):
    chroms = bw.chroms()
    if chrom in chroms:
        return chrom
    alt = chrom.replace("chr", "") if chrom.startswith("chr") else "chr" + chrom
    if alt in chroms:
        return alt
    raise ValueError(f"{chrom} not found in bigWig. Available examples: {list(chroms)[:5]}")


def binned_signal(bw_path, chrom, start, end, bins=800):
    with pyBigWig.open(str(bw_path)) as bw:
        bw_chrom = choose_chrom_for_bw(bw, chrom)
        edges = np.linspace(start, end, bins + 1).astype(int)
        vals = []
        for i in range(bins):
            s = int(edges[i])
            e = int(edges[i + 1])
            if e <= s:
                e = s + 1
            v = bw.stats(bw_chrom, s, e, type="mean")[0]
            vals.append(0 if v is None or np.isnan(v) else v)
        x = (edges[:-1] + edges[1:]) / 2
        return x, np.array(vals)


def load_peaks(peak_path, chrom, start, end):
    out = []
    if not peak_path.exists():
        return out
    with peak_path.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            p_chrom = fields[0]
            p_start = int(fields[1])
            p_end = int(fields[2])
            if same_chrom(p_chrom, chrom) and p_end > start and p_start < end:
                out.append((max(p_start, start), min(p_end, end)))
    return out


chrom, start, end, _, _ = read_bed(WINDOW_BED)
g_chrom, g_start, g_end, g_name, g_strand = read_bed(GENE_BED)

existing_tracks = [(label, path) for label, path in TRACKS if path.exists()]
if not existing_tracks:
    raise FileNotFoundError("No bigWig files found in bigwig/. Please check your file names.")

nrows = len(existing_tracks) + 2
fig, axes = plt.subplots(
    nrows,
    1,
    figsize=(12, 1.45 * nrows),
    sharex=True,
    gridspec_kw={"height_ratios": [1] * len(existing_tracks) + [0.55, 0.55]},
)

# bigWig signal tracks
all_max = 0
signals = []
for label, path in existing_tracks:
    x, y = binned_signal(path, chrom, start, end)
    signals.append((label, x, y))
    if len(y):
        all_max = max(all_max, float(np.nanmax(y)))

for ax, (label, x, y) in zip(axes[:len(existing_tracks)], signals):
    ax.fill_between(x, y, 0)
    ax.set_ylabel(label, rotation=0, ha="right", va="center")
    ax.set_ylim(0, all_max * 1.08 if all_max > 0 else 1)
    ax.spines[["top", "right"]].set_visible(False)

# peak bars
peak_ax = axes[len(existing_tracks)]
peak_ax.set_ylim(0, len(PEAKS))
peak_ax.set_yticks([])
peak_ax.set_ylabel("MACS2\npeaks", rotation=0, ha="right", va="center")
peak_ax.spines[["top", "right", "left"]].set_visible(False)

for i, (label, peak_path) in enumerate(PEAKS):
    peaks = load_peaks(peak_path, chrom, start, end)
    y = len(PEAKS) - i - 0.65
    peak_ax.text(start, y + 0.08, label, ha="left", va="bottom", fontsize=9)
    for ps, pe in peaks:
        peak_ax.add_patch(Rectangle((ps, y), pe - ps, 0.22))

# gene model
gene_ax = axes[len(existing_tracks) + 1]
gene_ax.set_ylim(0, 1)
gene_ax.set_yticks([])
gene_ax.set_ylabel("Gene", rotation=0, ha="right", va="center")
gene_ax.spines[["top", "right", "left"]].set_visible(False)

gene_ax.plot([g_start, g_end], [0.5, 0.5], linewidth=2)

if g_strand == "-":
    gene_ax.annotate("", xy=(g_start, 0.5), xytext=(g_end, 0.5),
                     arrowprops=dict(arrowstyle="->", linewidth=2))
    tss = g_end
else:
    gene_ax.annotate("", xy=(g_end, 0.5), xytext=(g_start, 0.5),
                     arrowprops=dict(arrowstyle="->", linewidth=2))
    tss = g_start

gene_ax.axvline(tss, linestyle="--", linewidth=1)
gene_ax.text((g_start + g_end) / 2, 0.68, f"{g_name} ({g_strand})", ha="center", va="bottom")

axes[-1].set_xlim(start, end)
axes[-1].set_xlabel(f"{chrom}:{start:,}-{end:,}  mm10 / GRCm38")

fig.suptitle("BATF ChIP-seq signal at Tnfrsf1b / TNFR2 locus", y=0.995)
fig.tight_layout()

out_pdf = BASE / "plots" / "Tnfrsf1b_BATF_chipseq_peakplot.pdf"
out_png = BASE / "plots" / "Tnfrsf1b_BATF_chipseq_peakplot.png"

fig.savefig(out_pdf)
fig.savefig(out_png, dpi=300)

print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
