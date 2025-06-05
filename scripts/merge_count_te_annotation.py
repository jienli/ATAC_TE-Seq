#!/usr/bin/env python3
import argparse
import pandas as pd
import pyranges as pr

def main(args):
    # 2) Load counts.tsv, skip the comment line, keep Chr/Start/End/Strand/Length
    counts = pd.read_csv(
        args.counts,
        sep='\t',
        comment='#',
        dtype={'Geneid': str}
    )
    counts = counts.rename(columns={'Geneid': 'peak'})
    # coords = counts[['peak', 'Chr', 'Start', 'End', 'Strand', 'Length']]
    coords = counts

    # 3) Use counts as the base peaks table
    merged = coords.rename(columns={
        'Chr': 'Chromosome',
        'Start': 'Start',
        'End': 'End',
        'Strand': 'Strand'
    })

    # 3) Load TE annotation and prefix its columns
    annot = pd.read_csv(args.annot)
    annot = annot.rename(columns={
        'gene_id': 'TE_gene_id',
        'seqnames': 'TE_seqnames',
        'type': 'TE_type',
        'start': 'TE_start',
        'end': 'TE_end',
        'strand': 'TE_strand',
        'family': 'TE_family',
        'element_start': 'TE_element_start',
        'element_end': 'TE_element_end',
        'pctdiv': 'TE_pctdiv',
        'length': 'TE_length',
        'repeat_superfamily': 'TE_repeat_superfamily',
        'rte_superfamily': 'TE_rte_superfamily',
        'rte_family': 'TE_rte_family',
        'rte_subfamily': 'TE_rte_subfamily',
        'l1_subfamily': 'TE_l1_subfamily',
        'herv_subfamily_limited': 'TE_herv_subfamily_limited',
        'req_integrative': 'TE_req_integrative',
        'genic_orientation_loc': 'TE_genic_orientation_loc',
        'loc_integrative': 'TE_loc_integrative'
    })
    annot = annot.rename(columns={'TE_seqnames': 'Chromosome', 'TE_start': 'Start', 'TE_end': 'End', 'TE_strand': 'Strand'})
    
    # 4) Merge on genomic coordinates by overlap
    # Create PyRanges objects for peaks and TE annotations
    peaks_gr = pr.PyRanges(df=merged)
    te_gr = pr.PyRanges(df=annot)
    # Perform overlap join
    joined = peaks_gr.join(te_gr)
    final_multi = joined.df

    # Compute overlap lengths and pick TE with longest overlap per peak
    df = joined.df.copy()
    df['overlap_len'] = (
        df.apply(
            lambda row: min(row['End'], row['End_b']) - max(row['Start'], row['Start_b']),
            axis=1
        )
        .clip(lower=0)
    )
    # Select only the TE with the longest overlap for each peak
    final = (
        df.sort_values(['peak', 'overlap_len'], ascending=[True, False])
          .drop_duplicates(subset=['peak'], keep='first')
          .drop(columns=['overlap_len'])
    )



    
    # 5) Write out
    final.to_csv(args.out, index=False)
    final_multi.to_csv(args.out.replace('.csv', '_multi.csv'), index=False)
    print(f"Wrote combined data to {args.out} and {args.out.replace('.csv', '_multi.csv')}")
    print(f"Number of peaks with TE annotation: {len(final)}")
    print(f"Number of peaks with multiple TE annotations: {len(final_multi)}") 

if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description="Merge peak counts with TE annotation"
    )
    p.add_argument('--counts', required=True,
                   help="TSV from featureCounts (has Geneid,Chr,Start,End,Strand,Length)")
    p.add_argument('--annot', required=True,
                   help="CSV of TE annotation (has seqnames,start,end, etc.)")
    p.add_argument('--out', default='merged_results.csv',
                   help="Output CSV filename")
    args = p.parse_args()
    main(args)