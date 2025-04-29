from functools import partial
from typing import Dict, Tuple, Union

import numpy as np
import pandas as pd
import scipy.stats as ss
import statsmodels.api as sm
from tqdm.contrib.concurrent import process_map

from ..api.formats import Rdc, DnaDataset
from ..api.mp import adjust_chunksize


def _calculate_pvals_single(rna_name: str, total_contacts: int, rdc_data: Rdc) -> Tuple[str, pd.DataFrame]:
    pixels = rdc_data.read_pixels(rna_name, value_fields=['signal_count', 'bg_prob'])
    pvals = pixels.apply(lambda row: np.nan if row['signal_count'] == 0 else ss.binomtest(k=row['signal_count'], n=total_contacts, p=row['bg_prob'], alternative='greater').pvalue, axis=1)
    pixels['pvalue'] = pvals
    return rna_name, pixels[['chrom', 'pvalue']]


def _calculate_pvals(rdc_data: Rdc,
                     n_cores: int = 1,
                     chunksize: int = 50) -> Dict[str, pd.DataFrame]:
    cis_contacts_num = rdc_data.read_rna_attribute_batch('cis_contacts')
    trans_contacts_num = rdc_data.read_rna_attribute_batch('trans_contacts')
    rna_names = list(cis_contacts_num.keys())
    total_contacts_gen = (cis_contacts_num[rna_name] + trans_contacts_num[rna_name] for rna_name in rna_names)
    func = partial(_calculate_pvals_single, rdc_data=rdc_data)
    chunksize = adjust_chunksize(len(rna_names), n_cores, chunksize)
    results = dict(process_map(func,
                               rna_names,
                               total_contacts_gen,
                               max_workers=n_cores,
                               chunksize=chunksize,
                               desc='Calculating p-values',
                               unit='RNA'))
    return results


def estimate_significance(rdc_data: Rdc, n_cores: int = 1) -> None:
    total_stats = _calculate_pvals(rdc_data, n_cores)

    rna_dfs_handler = list()
    for rna_name, rna_stats in total_stats.items():
        pvals = rna_stats.dropna().reset_index(drop=True)['pvalue']
        rna_df = pd.DataFrame({'rna_name': rna_name, 'pvalue': pvals})
        _, qvals_rna, *_ = sm.stats.multipletests(rna_df['pvalue'].values, method='fdr_bh')
        rna_df['qvalue_rna'] = qvals_rna
        rna_dfs_handler.append(rna_df)
    meaningful_pvals = pd.concat(rna_dfs_handler, ignore_index=True)

    _, qvals, *_ = sm.stats.multipletests(meaningful_pvals['pvalue'].values, method='fdr_bh')

    meaningful_pvals['qvalue_global'] = qvals

    for rna_name, rna_df in meaningful_pvals.groupby('rna_name'):
        rna_stats = total_stats[rna_name]
        rna_stats.loc[~rna_stats['pvalue'].isna(), 'qvalue_global'] = rna_df['qvalue_global'].values
        rna_stats.loc[~rna_stats['pvalue'].isna(), 'qvalue_rna'] = rna_df['qvalue_rna'].values

    rdc_data.write_pixels_column_batch('pvalue', total_stats)
    rdc_data.write_pixels_column_batch('qvalue_global', total_stats)
    rdc_data.write_pixels_column_batch('qvalue_rna', total_stats)
    rdc_data.are_peaks_estimated = True


def _fetch_peaks_single(rna_name: str, rdc_data: Rdc, threshold: float = 0.05, qvalue_type='global') -> pd.DataFrame:
    pixels = rdc_data.read_pixels(rna_name)
    if qvalue_type not in ['global', 'rna']:
        raise ValueError
    qvalue_col = f'qvalue_{qvalue_type}'
    peaks = pixels.query(f'{qvalue_col} < {threshold}').reset_index(drop=True).rename(columns={qvalue_col: 'qvalue'})
    peaks['rna_name'] = rna_name
    return peaks


def fetch_peaks(rdc_data: Rdc, threshold: float = 0.05, qvalue_type='global', n_cores: int = 1, chunksize: int = 50) -> pd.DataFrame:
    if not rdc_data.are_peaks_estimated:
        raise Exception
    if qvalue_type not in ('global', 'rna'):
        raise ValueError
    annotation = rdc_data.annotation
    rna_names = list(annotation.keys())
    func = partial(_fetch_peaks_single, rdc_data=rdc_data, threshold=threshold, qvalue_type=qvalue_type)
    chunksize = adjust_chunksize(len(rna_names), n_cores, chunksize)
    results = process_map(func,
                          rna_names,
                          max_workers=n_cores,
                          chunksize=chunksize,
                          desc='Fetching peaks',
                          unit='RNA')
    return pd.concat(results, ignore_index=True)

def _fetch_top_percent_single(rna_name: str, rdc_data: Rdc, qvalue_type: str = 'global', top_percent: float = 10.0) -> pd.DataFrame:
    pixels = rdc_data.read_pixels(rna_name)
    qvalue_col = f'qvalue_{qvalue_type}'
    pixels = pixels.dropna(subset=[qvalue_col])
    top_percent_contacts = pixels.sort_values(by=qvalue_col, ascending=True).head(int(pixels.shape[0] * top_percent / 100)).reset_index(drop=True).rename(columns={qvalue_col: 'qvalue'})
    top_percent_contacts['rna_name'] = rna_name
    return top_percent_contacts

def fetch_top_percent(rdc_data: Rdc, qvalue_type: str = 'global', top_percent: float = 10.0, n_cores: int = 1, chunksize: int = 50) -> pd.DataFrame:
    if not rdc_data.are_peaks_estimated:
        raise Exception
    rna_names = list(rdc_data.annotation.keys())
    func = partial(_fetch_top_percent_single, rdc_data=rdc_data, qvalue_type=qvalue_type, top_percent=top_percent)
    chunksize = adjust_chunksize(len(rna_names), n_cores, chunksize)
    results = process_map(func,
                          rna_names,
                          max_workers=n_cores,
                          chunksize=chunksize,
                          desc='Fetching top percent',
                          unit='RNA')
    return pd.concat(results, ignore_index=True)

def format_peaks(peaks_df: pd.DataFrame, format: str = 'narrowPeak', score: Union[str, int] = 0) -> pd.DataFrame:
    if isinstance(score, int):
        peaks_df['score'] = score
    elif isinstance(score, str):
        if score in peaks_df:
            peaks_df['score'] = peaks_df[score]
        else:
            raise ValueError
    else:
        raise ValueError

    peaks_df['strand'] = '.'

    if format == 'narrowPeak':
        peaks_df['peak'] = -1
        peaks_df['-log10pval'] = -np.log10(peaks_df['pvalue'])
        peaks_df['-log10qval'] = -np.log10(peaks_df['qvalue'])
        return peaks_df[['chrom', 'start', 'end', 'rna_name', 'score', 'strand', 'fc', '-log10pval', '-log10qval', 'peak']]
    elif format == 'bed':
        return peaks_df[['chrom', 'start', 'end', 'rna_name', 'fc', 'strand']]
    else:
        raise ValueError


