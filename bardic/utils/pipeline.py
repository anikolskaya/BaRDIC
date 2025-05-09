from typing import Dict, List

import pandas as pd

from ..api.io import read_bedgraph
from .background import make_background_track
from .binsizes import optimize_bin_sizes
from .dnadataset import bed2h5, encode_invalid_rna_names, decode_invalid_rna_names
from .peaks import estimate_significance, fetch_peaks, format_peaks, fetch_top_percent
from .rdc import dnadataset_to_rdc
from .scaling import calculate_scaling_splines


def run_pipeline(dna_parts_fname: str,
                 dna_dataset_fname: str,
                 rdc_fname: str,
                 chromsizes: Dict[str, int],
                 annotation: pd.DataFrame,
                 selection_results_fname: str,
                 bg_fname: str,
                 rna_list: List,
                 bg_binsize: int,
                 qval_threshold: float,
                 qval_type: str,
                 top_percent: float,
                 peaks_output: str,
                 binsize_params: Dict = {},
                 rdc_params: Dict = {},
                 scaling_params: Dict = {},
                 peaks_format_params: Dict = {},
                 make_bg: bool = True,
                 uniform_bg: bool = False,
                 n_cores: int = 1):
    chromdict = chromsizes

    dna_dataset, invalid_rna_names = bed2h5(dna_parts_fname, dna_dataset_fname, chromdict, annotation)

    selection_df = optimize_bin_sizes(dna_dataset, n_cores=n_cores, **binsize_params)
    
    rna_list = [encode_invalid_rna_names(rna) for rna in rna_list]
    if len(invalid_rna_names) > 0:
        selection_df['gene_name'] = selection_df['gene_name'].apply(decode_invalid_rna_names)

    selection_df.to_csv(selection_results_fname, sep='\t', header=True, index=False)

    if make_bg:
        bg_track = make_background_track(dna_dataset, rna_list, bg_binsize, uniform_bg)
        bg_track.to_csv(bg_fname, header=False, index=False, sep='\t')
    else:
        bg_track = read_bedgraph(bg_fname)

    rdc_data = dnadataset_to_rdc(dna_dataset, bg_track, rdc_fname, n_cores=n_cores, **rdc_params)

    calculate_scaling_splines(rdc_data, n_cores=n_cores, **scaling_params)

    estimate_significance(rdc_data, n_cores)
    peaks = fetch_peaks(rdc_data, qval_threshold, qval_type, n_cores)
    fixed_peaks = fetch_top_percent(rdc_data, qval_type, top_percent, n_cores)

    formatted_peaks = format_peaks(peaks, **peaks_format_params)
    fixed_formatted_peaks = format_peaks(fixed_peaks, **peaks_format_params)
    

    if len(invalid_rna_names) > 0:
        formatted_peaks.iloc[:, 3] = formatted_peaks.iloc[:, 3].apply(decode_invalid_rna_names)
        fixed_formatted_peaks.iloc[:, 3] = fixed_formatted_peaks.iloc[:, 3].apply(decode_invalid_rna_names)

    formatted_peaks.to_csv(peaks_output, sep='\t', header=False, index=False)
    fixed_formatted_peaks.to_csv(peaks_output.replace(peaks_format_params['format'], 
                                                      f'top_{int(top_percent)}perc.{peaks_format_params["format"]}'), sep='\t', header=False, index=False)