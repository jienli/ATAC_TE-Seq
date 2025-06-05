#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

def violin_scatter_box(df, val_col, group_col, suffix='', log_scale=False, figure_size_multiplier=1.5):
    """
    Generate a violin plot with overlaid boxplot and density-based scatter points.

    Parameters:
        df (pandas.DataFrame): DataFrame containing the data.
        val_col (str): Column name for the numeric values to plot.
        group_col (str): Column name for the categorical grouping variable.
        suffix (str): Optional filename suffix.
        log_scale (bool): Whether to apply log scale to the y-axis.
        figure_size_multiplier (float): Scale factor for figure size.

    Returns:
        matplotlib.figure.Figure: The created figure object.
    """
    # Drop missing values
    df = df.dropna(subset=[val_col, group_col])
    # If log scale, filter out non-positive counts
    if log_scale:
        df = df[df[val_col] > 0]
    # Prepare data
    groups = sorted(df[group_col].unique())
    data = [df.loc[df[group_col] == g, val_col].values for g in groups]
    # Create figure
    fig = plt.figure(figsize=(max(6, len(groups)*0.5) * figure_size_multiplier,
                              5 * figure_size_multiplier))
    ax = fig.add_subplot(1,1,1)
    # Violin
    parts = ax.violinplot(data, positions=range(1, len(groups)+1), showextrema=False, widths=0.8)
    for body in parts['bodies']:
        body.set_facecolor('lightgray'); body.set_alpha(0.3)
    # Boxplot
    ax.boxplot(data, positions=range(1, len(groups)+1), widths=0.6,
               showfliers=False, patch_artist=True,
               boxprops=dict(facecolor='lightgray', alpha=0.5))
    # Scatter
    width = 0.4
    for i, g in enumerate(groups, start=1):
        y = df.loc[df[group_col] == g, val_col].values
        if len(y) == 0:
            continue
        # always normalize linearly as per your selection
        if log_scale:
            vals = np.log10(y + 1)
        else:
            vals = y
        vals = vals / vals.max()
        counts_hist, bin_edges = np.histogram(vals, bins=50, density=True)
        smoothed = gaussian_filter1d(counts_hist, sigma=1)
        bin_idx = np.digitize(vals, bin_edges) - 1
        bin_idx = np.clip(bin_idx, 0, len(smoothed)-1)
        densities = smoothed[bin_idx]
        scales = densities / densities.max() * width if densities.max() > 0 else np.full_like(densities, width)
        # x = i + np.random.uniform(-scales, scales)
        # Evenly distribute points across [-scale, scale] for each density value
        x = np.empty_like(vals)
        for s in np.unique(scales):
            idxs = np.where(scales == s)[0]
            k = len(idxs)
            if k > 0:
                positions = np.linspace(-s, s, k)
                x[idxs] = i + positions
        ax.scatter(x, y, alpha=0.6, edgecolors='k', linewidths=0.2, s=20 * figure_size_multiplier)
    # Axes styling
    ax.set_xticks(range(1, len(groups)+1))
    ax.set_xticklabels(groups, rotation=45, ha='right')
    if log_scale:
        ax.set_yscale('log')
        ax.set_ylabel('Read Count (log scale)')
    else:
        ax.set_ylabel('Read Count')
    ax.set_xlabel('TE Subfamily')
    ax.set_title('Counts by TE Subfamily')
    fig.tight_layout()
    return fig

def plot_counts_by_subfamily(df, count_col, subfamily_col, outdir, suffix='', log_scale=False):
    # Delegate plotting
    fig = violin_scatter_box(df, count_col, subfamily_col, suffix, log_scale)
    os.makedirs(outdir, exist_ok=True)
    filename = 'counts_by_subfamily' + (f'_{suffix}' if suffix else '') + '.png'
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Plot read counts by TE subfamily")
    parser.add_argument('--input', default='analysis/merged_results_count.csv',
                        help="Path to merged_results_count.csv")
    parser.add_argument('--count-col', default='filtered/SRR13601447.PE.final.bam',
                        help="Column name containing raw counts")
    parser.add_argument('--subfamily-col', default='TE_rte_subfamily',
                        help="Column name for TE subfamily annotation")
    parser.add_argument('--outdir', default='analysis/plots',
                        help="Directory to write the plot")
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.input)
    plot_counts_by_subfamily(df, args.count_col, args.subfamily_col, args.outdir)

    # Also plot specified L1 subfamilies with 'L1' suffix
    l1_subfamilies = ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]
    df_l1 = df[df[args.subfamily_col].isin(l1_subfamilies)]
    if not df_l1.empty:
        plot_counts_by_subfamily(df_l1, args.count_col, args.subfamily_col, args.outdir, 'L1')

    # Also generate log-scale plots
    plot_counts_by_subfamily(df, args.count_col, args.subfamily_col, args.outdir, suffix='log', log_scale=True)
    if not df_l1.empty:
        plot_counts_by_subfamily(df_l1, args.count_col, args.subfamily_col, args.outdir, suffix='L1_log', log_scale=True)

if __name__ == '__main__':
    main()