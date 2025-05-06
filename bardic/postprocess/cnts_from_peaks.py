import pandas as pd
import bioframe as bf
import argparse
import os
from pathlib import Path

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='Merge contact data with peaks')
    parser.add_argument('-p', '--peaks', required=True, help='Path to peaks file')
    parser.add_argument('-c', '--cnts', required=True, help='Path to contacts file')
    parser.add_argument('-o', '--out_dir', required=True, help='Output directory')
    args = parser.parse_args()

    # Create output directory if needed
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    # Define columns to keep
    peak_cols = ['chrom', 'start', 'end', 'name', 'fc', '-log10p', '-log10q']
    contact_cols = {6: 'name', 12: 'chrom', 13: 'start', 14: 'end'}

    # Read and process files
    peaks = bf.read_table(args.peaks, schema='narrowPeak')[peak_cols]
    cnts = pd.read_csv(args.cnts, sep='\t', header=None)
    cnts = cnts.rename(columns=contact_cols)
    cnts.columns = [str(x) for x in cnts.columns]
    cnts['end_tmp'] = cnts['start'] + 1  # Represent contact as point

    # Merge data
    merged = bf.overlap(
        cnts, 
        peaks, 
        how='left', 
        suffixes=('', '_peak'),
        cols1=['chrom', 'start', 'end_tmp'],
        on=['name']
    ).drop(columns=['name_peak', 'end_tmp'])

    # Save output
    out_name = f"{Path(args.cnts).stem.split('.')[0]}.from_peaks.bed"
    merged.to_csv(Path(args.out_dir) / out_name, sep='\t', index=False)

if __name__ == '__main__':
    main()