#!/usr/bin/env python3
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

def summarize_and_plot(df, var, outdir):
    # drop missing
    dfv = df.dropna(subset=[var])
    groups = sorted(dfv[var].unique())

    # Summary statistics
    stats_fc = dfv.groupby(var)['log2FoldChange'] \
                  .agg(['count','mean','median']) \
                  .rename(columns={'count':'n_peaks'}) 
    stats_padj = dfv.groupby(var)['neglog10_padj'] \
                    .agg(['count','mean','median']) \
                    .rename(columns={'count':'n_peaks'})
    
    stats_fc.to_csv(os.path.join(outdir, f"stats_log2FC_by_{var}.csv"))
    stats_padj.to_csv(os.path.join(outdir, f"stats_neglog10_padj_by_{var}.csv"))

    # Define pairwise comparisons for p-value annotation
    comparisons = []
    if var == 'TE_req_integrative':
        comparisons = [('Yng Intact', 'Yng FL')]
    elif var == 'TE_rte_subfamily':
        comparisons = [('L1HS', 'L1PA4'),('L1HS', 'L1PA5'), ('L1HS', 'L1PA6')]
    elif var == 'TE_genic_orientation_loc':
        comparisons = [('Antisense', 'Intergenic'), ('Intergenic', 'Sense')]

    # Boxplot: log2FoldChange
    data_fc = [dfv.loc[dfv[var]==g, 'log2FoldChange'].values for g in groups]
    positions = range(1, len(groups)+1)
    plt.figure(figsize=(max(6, len(groups)*0.5), 4))
    plt.boxplot(data_fc, positions=positions, labels=groups, showfliers=False)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('log2 Fold Change')
    plt.title(f'log2FC by {var}')
    ymin, ymax = plt.ylim()
    offset = (ymax - ymin) * 0.1
    # add p-value annotations
    for idx, (g1, g2) in enumerate(comparisons):
        if g1 in groups and g2 in groups:
            x1, x2 = groups.index(g1)+1, groups.index(g2)+1
            y1 = max(data_fc[groups.index(g1)]) if len(data_fc[groups.index(g1)])>0 else 0
            y2 = max(data_fc[groups.index(g2)]) if len(data_fc[groups.index(g2)])>0 else 0
            base = max(y1, y2)
            y = ymax + offset * (idx)
            pval = ttest_ind(
                dfv.loc[dfv[var]==g1, 'log2FoldChange'],
                dfv.loc[dfv[var]==g2, 'log2FoldChange'],
                nan_policy='omit'
            ).pvalue
            plt.plot([x1, x2], [y, y], linewidth=1)
            plt.text((x1+x2)/2, y, f"p={pval:.2e}", ha='center', va='bottom')
    plt.ylim(ymin, ymax + offset * (len(comparisons) + 0.3))
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"log2FC_by_{var}.png"))
    plt.close()

    # Boxplot: –log10(padj)
    data_p = [dfv.loc[dfv[var]==g, 'neglog10_padj'].values for g in groups]
    plt.figure(figsize=(max(6, len(groups)*0.5), 4))
    plt.boxplot(data_p, positions=positions, labels=groups, showfliers=False)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('-log10(padj)')
    plt.title(f'-log10(padj) by {var}')
    ymin, ymax = plt.ylim()
    offset = (ymax - ymin) * 0.1
    # add p-value annotations
    for idx, (g1, g2) in enumerate(comparisons):
        if g1 in groups and g2 in groups:
            x1, x2 = groups.index(g1)+1, groups.index(g2)+1
            y1 = max(data_p[groups.index(g1)]) if len(data_p[groups.index(g1)])>0 else 0
            y2 = max(data_p[groups.index(g2)]) if len(data_p[groups.index(g2)])>0 else 0
            base = max(y1, y2)
            y = ymax + offset * (idx)

            pval = ttest_ind(
                dfv.loc[dfv[var]==g1, 'neglog10_padj'],
                dfv.loc[dfv[var]==g2, 'neglog10_padj'],
                nan_policy='omit'
            ).pvalue
            plt.plot([x1, x2], [y, y], linewidth=1)
            plt.text((x1+x2)/2, y, f"p={pval:.2e}", ha='center', va='bottom')
    plt.ylim(ymin, ymax + offset * (len(comparisons) + 0.3))
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"neglog10_padj_by_{var}.png"))
    plt.close()


def main():
    p = argparse.ArgumentParser(
        description="Group-level analysis of merged peak enrichment"
    )
    p.add_argument('--input', required=True,
                   help="Path to merged_results.csv")
    p.add_argument('--outdir', default='analysis',
                   help="Directory to write stats and plots")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load data
    df = pd.read_csv(args.input)
    # compute –log10(padj), avoiding log(0)
    df['neglog10_padj'] = -np.log10(df['padj'].replace(0, np.nextafter(0,1)))

    # grouping variables
    grouping_vars = [
        'TE_rte_subfamily',
        'TE_req_integrative',
        'TE_genic_orientation_loc'
    ]

    for var in grouping_vars:
        print(f"Processing grouping by {var}…")
        summarize_and_plot(df, var, args.outdir)

    # Heatmaps comparing TE_req_integrative vs TE_rte_subfamily
    cmap = 'Reds'
    combo_df = df.dropna(subset=['TE_req_integrative', 'TE_rte_subfamily'])
    # Pivot for mean log2FoldChange
    pivot_fc = combo_df.pivot_table(
        index='TE_req_integrative',
        columns='TE_rte_subfamily',
        values='log2FoldChange',
        aggfunc='mean'
    )
    plt.figure(figsize=(8, 6))
    plt.imshow(pivot_fc.values, aspect='auto', cmap=cmap)
    plt.colorbar(label='Mean log2 Fold Change')
    plt.xticks(range(len(pivot_fc.columns)), pivot_fc.columns, rotation=45, ha='right')
    plt.yticks(range(len(pivot_fc.index)), pivot_fc.index)
    plt.title('Mean log2FC by req_integrative vs rte_subfamily')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'heatmap_log2FC_req_vs_rte.png'))
    plt.close()

    # Pivot for mean -log10(padj)
    pivot_p = combo_df.pivot_table(
        index='TE_req_integrative',
        columns='TE_rte_subfamily',
        values='neglog10_padj',
        aggfunc='mean'
    )
    plt.figure(figsize=(8, 6))
    plt.imshow(pivot_p.values, aspect='auto', cmap=cmap)
    plt.colorbar(label='Mean -log10(padj)')
    plt.xticks(range(len(pivot_p.columns)), pivot_p.columns, rotation=45, ha='right')
    plt.yticks(range(len(pivot_p.index)), pivot_p.index)
    plt.title('Mean -log10(padj) by req_integrative vs rte_subfamily')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'heatmap_padj_req_vs_rte.png'))
    plt.close()

    # Heatmaps for L1 subfamily only
    l1_list = ['L1HS','L1PA2','L1PA3','L1PA4','L1PA5','L1PA6']
    combo_l1 = combo_df[combo_df['TE_rte_subfamily'].isin(l1_list)]
    # Mean log2FC
    pivot_fc_l1 = combo_l1.pivot_table(
        index='TE_req_integrative',
        columns='TE_rte_subfamily',
        values='log2FoldChange',
        aggfunc='mean'
    )
    plt.figure(figsize=(8, 6))
    plt.imshow(pivot_fc_l1.values, aspect='auto', cmap=cmap)
    plt.colorbar(label='Mean log2 Fold Change')
    plt.xticks(range(len(pivot_fc_l1.columns)), pivot_fc_l1.columns, rotation=45, ha='right')
    plt.yticks(range(len(pivot_fc_l1.index)), pivot_fc_l1.index)
    plt.title('Mean log2FC (L1 only) by req_integrative vs rte_subfamily')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'heatmap_log2FC_req_vs_rte_L1.png'))
    plt.close()
    # Mean -log10(padj)
    pivot_p_l1 = combo_l1.pivot_table(
        index='TE_req_integrative',
        columns='TE_rte_subfamily',
        values='neglog10_padj',
        aggfunc='mean'
    )
    plt.figure(figsize=(8, 6))
    plt.imshow(pivot_p_l1.values, aspect='auto', cmap=cmap)
    plt.colorbar(label='Mean -log10(padj)')
    plt.xticks(range(len(pivot_p_l1.columns)), pivot_p_l1.columns, rotation=45, ha='right')
    plt.yticks(range(len(pivot_p_l1.index)), pivot_p_l1.index)
    plt.title('Mean -log10(padj) (L1 only) by req_integrative vs rte_subfamily')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'heatmap_padj_req_vs_rte_L1.png'))
    plt.close()

    # Heatmaps comparing TE_genic_orientation_loc vs TE_req_integrative
    combo2 = df.dropna(subset=['TE_genic_orientation_loc', 'TE_req_integrative'])
    # Mean log2FC
    pivot2_fc = combo2.pivot_table(
        index='TE_genic_orientation_loc',
        columns='TE_req_integrative',
        values='log2FoldChange',
        aggfunc='mean'
    )
    plt.figure(figsize=(8, 6))
    plt.imshow(pivot2_fc.values, aspect='auto', cmap=cmap)
    plt.colorbar(label='Mean log2 Fold Change')
    plt.xticks(range(len(pivot2_fc.columns)), pivot2_fc.columns, rotation=45, ha='right')
    plt.yticks(range(len(pivot2_fc.index)), pivot2_fc.index)
    plt.title('Mean log2FC by genic_orientation_loc vs req_integrative')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'heatmap_log2FC_genic_vs_req.png'))
    plt.close()
    # Mean -log10(padj)
    pivot2_p = combo2.pivot_table(
        index='TE_genic_orientation_loc',
        columns='TE_req_integrative',
        values='neglog10_padj',
        aggfunc='mean'
    )
    plt.figure(figsize=(8, 6))
    plt.imshow(pivot2_p.values, aspect='auto', cmap=cmap)
    plt.colorbar(label='Mean -log10(padj)')
    plt.xticks(range(len(pivot2_p.columns)), pivot2_p.columns, rotation=45, ha='right')
    plt.yticks(range(len(pivot2_p.index)), pivot2_p.index)
    plt.title('Mean -log10(padj) by genic_orientation_loc vs req_integrative')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'heatmap_padj_genic_vs_req.png'))
    plt.close()

    # Heatmaps for L1 subfamily only (genic vs req)
    combo2_l1 = combo2[combo2['TE_rte_subfamily'].isin(l1_list)]
    pivot2_fc_l1 = combo2_l1.pivot_table(
        index='TE_genic_orientation_loc',
        columns='TE_req_integrative',
        values='log2FoldChange',
        aggfunc='mean'
    )
    plt.figure(figsize=(8, 6))
    plt.imshow(pivot2_fc_l1.values, aspect='auto', cmap=cmap)
    plt.colorbar(label='Mean log2 Fold Change')
    plt.xticks(range(len(pivot2_fc_l1.columns)), pivot2_fc_l1.columns, rotation=45, ha='right')
    plt.yticks(range(len(pivot2_fc_l1.index)), pivot2_fc_l1.index)
    plt.title('Mean log2FC (L1 only) by genic_orientation_loc vs req_integrative')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'heatmap_log2FC_genic_vs_req_L1.png'))
    plt.close()
    pivot2_p_l1 = combo2_l1.pivot_table(
        index='TE_genic_orientation_loc',
        columns='TE_req_integrative',
        values='neglog10_padj',
        aggfunc='mean'
    )
    plt.figure(figsize=(8, 6))
    plt.imshow(pivot2_p_l1.values, aspect='auto', cmap=cmap)
    plt.colorbar(label='Mean -log10(padj)')
    plt.xticks(range(len(pivot2_p_l1.columns)), pivot2_p_l1.columns, rotation=45, ha='right')
    plt.yticks(range(len(pivot2_p_l1.index)), pivot2_p_l1.index)
    plt.title('Mean -log10(padj) (L1 only) by genic_orientation_loc vs req_integrative')
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, 'heatmap_padj_genic_vs_req_L1.png'))
    plt.close()

    print("Done. Summary tables and plots are in", args.outdir)


if __name__ == '__main__':
    main()