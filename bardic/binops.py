from typing import Dict
import bioframe as bf
import numpy as np
import pandas as pd


def make_geom_bins(length, start_size, factor):
    """Makes geometrically increasing bins.

    Args:
        length (int): chromosome lenth.
        start_size (int): start bin size, base of the geometric progression.
        factor (float, int): a common ratio between two successive bin sizes.
            Must be not less than 1.

    Returns:
        np.array: bins edges.

    Raises:
        ValueError: in case factor value is less than 1.
    """
    if factor < 1:
        raise ValueError(f"factor value {factor} is less than 1.")
    elif factor == 1:
        bin_edges = np.arange(0, length, start_size, dtype='int')
    else:
        final_index = round(np.log(1 - length * (1 - factor) / start_size) / np.log(factor) - 1)  # from sum of geometric series
        bin_edges = np.cumsum(start_size * factor ** np.arange(0, final_index + 1)).astype('int')
        bin_edges = np.insert(bin_edges, 0, 0)
    bin_edges[-1] = length
    return bin_edges


def prune_geom_bins(geom_bins, max_linear_size):
    """Prunes geometric bins to maximal linear size.

    Given output from `make_geom_bins` removes all bins
    that are longer than `max_linear_size` and fills the
    remaining chromosome length with bins of `max_linear_size`
    in length.

    Args:
        geom_bins (np.array): bins edges from `make_geom_bins`.
        max_linear_size (int): maximal linear size.
    """
    chrom_length = geom_bins[-1]
    scaled_bins = geom_bins[:np.sum(np.diff(geom_bins) <= max_linear_size) + 1]
    if len(scaled_bins) == len(geom_bins):
        return geom_bins
    else:
        linear_bins = np.array(range(scaled_bins[-1], chrom_length, max_linear_size))
        linear_bins[-1] = chrom_length
        return np.concatenate((scaled_bins, linear_bins[1:]))


def make_cis_bins(factor, start_size, chrom_name, chrom_length, gene_start, gene_end, max_linear_size=None, fillgene=False):
    """Makes cis bins of increasing size from gene borders in geometric progression.

    Bin i is of size start_size * factor ** i.

    Args:
        factor (int, float): geometric progression factor.
        start_size (int): start bin size.
        chrom_name (str): chromosome name to put in the dataframe.
        chrom_length (int): length of the chromosome.
        gene_start (int): start coordinate of the gene.
        gene_end (int): end coordinate of the gene.
        min_bin_size (int): minimal linear bin size to
            preserve the bin.
        max_linear_size (int): maximal linear bin size to stop
            geometric growth of bin size and use given bin size
            hereafter (default: None).
        fillgene (bool): if True, will add a single bin of gene
            coordinates, otherwise will exclude gene from binning
            (default: False).

    Returns:
        pd.DataFrame: cis bins as a dataframe with 3 columns:
            chrom, start, end.

    Raises:
        ValueError: in case factor < 1.
    """
    if factor < 1:
        raise ValueError(f"factor value {factor} is less than 1.")
    upstream_size = gene_start
    downstream_size = chrom_length - gene_end

    upstream_edges = make_geom_bins(upstream_size, start_size, factor)
    if max_linear_size is not None:
        upstream_edges = prune_geom_bins(upstream_edges, max_linear_size)
    upstream_edges = gene_start - upstream_edges[::-1]
    upstream_df = pd.DataFrame({'chrom': chrom_name, 'start': upstream_edges[:-1], 'end': upstream_edges[1:]})

    downstream_edges = make_geom_bins(downstream_size, start_size, factor)
    if max_linear_size is not None:
        downstream_edges = prune_geom_bins(downstream_edges, max_linear_size)
    downstream_edges = gene_end + downstream_edges
    downstream_df = pd.DataFrame({'chrom': chrom_name, 'start': downstream_edges[:-1], 'end': downstream_edges[1:]})

    if fillgene:
        gene_df = pd.DataFrame({'chrom': chrom_name, 'start': gene_start, 'end': gene_end}, index=[0])
        bins_df = pd.concat((upstream_df, gene_df, downstream_df), ignore_index=True)
    else:
        bins_df = pd.concat((upstream_df, downstream_df), ignore_index=True)
    return bins_df


make_cis_bins3 = make_cis_bins


def make_linear_bins(bin_size: int, chromsizes: Dict[str, int]):
    """Makes equal sized bins in linear scale for given chromosomes.

    Args:
        bin_size (int): linear bin size.
        chromsizes_df (pd.DataFrame): chromosome sizes dataframe
            with 2 columns: "chrom" (chromosome name) and "length"
            (chromosome length).

    Returns:
        pd.DataFrame: a dataframe with bins of 3 columns:
            "chrom", "start", "end".
    """
    chrom_dfs = list()
    bin_size = round(bin_size)
    for chrom, length in chromsizes.items():
        break_points = list(range(0, length, bin_size))
        starts = break_points[:-1]
        ends = break_points[1:]
        if len(break_points) <= 1:
            starts = [0]
            ends = [length]
        else:
            if ends[-1] < length:
                starts.append(ends[-1])
                ends.append(length)
        chrom_df = pd.DataFrame(((chrom, start, end)
                                 for start, end in zip(starts, ends)),
                                columns=['chrom', 'start', 'end'])
        chrom_dfs.append(chrom_df)
    total_bins = pd.concat(chrom_dfs, ignore_index=True)
    total_bins['chrom'] = total_bins['chrom'].astype('category')
    return total_bins


def make_trans_bins(bin_size: int, chromsizes: Dict[str, int], gene_chrom: str):
    """Makes trans bins of equal size.

    Uses `make_linear_bins` for binning.

    Args:
        bin_size (int): linear bin size.
        chromsizes_df (pd.DataFrame): chromosome sizes dataframe
            with 2 columns: "chrom" (chromosome name) and "length"
            (chromosome length).
        gene_chrom (str): a gene chromosome to exclude it
            from binning.

    Returns:
        pd.DataFrame: a dataframe with bins of 3 columns:
            "chrom", "start", "end".
    """
    _chromsizes = chromsizes.copy()
    del _chromsizes[gene_chrom]
    return make_linear_bins(bin_size, _chromsizes)


def calculate_bins_coverage(bins_df, contacts_df):
    """Calculates how many contacts overlap every bin.

    Args:
        contacts_df (pd.DataFrame): a contacts dataframe
            of 3 columns: "chrom", "start", "end".
        bins_df (pd.DataFrame): a bins dataframe of
            3 columns: "chrom", "start", "end".

    Returns:
        pd.DataFrame: a 4-column dataframe
            ("chrom", "start", "end", "count") with
            bins from bins_df and amount of contacts
            overlapping with each of them.
    """
    return bf.count_overlaps(bins_df, contacts_df)


def interval_centers(intervals_df):
    """Make annotation of centers of provided intervals.

    Args:
        intervals_df (pd.DataFrame): a 3-column dataframe with
            genomic intervals: "chrom", "start", "end".

    Returns:
        pd.DataFrame: a 3-column dataframe ("chrom", "start", "end")
            of interval centers.
    """
    intervals_centers = (intervals_df['end'] + intervals_df['start']) // 2
    intervals_centers_df = pd.DataFrame({'chrom': intervals_df['chrom'],
                                         'start': intervals_centers,
                                         'end': intervals_centers + 1})
    return intervals_centers_df


make_interval_centers = interval_centers


def make_rel_dist(point,
                  start,
                  end):
    """Finds relative distance of a point to gene with [start; end] coordinates.

    Returns:
        int: relative distance.
    """
    if point < start:
        return start - point
    elif point > end:
        return point - end
    else:
        return 0


def make_rel_dist_vector(points,
                         start,
                         end):
    """Finds relative distance of points to gene with [start; end] coordinates.

    Args:
        points (np.array): genomic points.
        start (int): gene start.
        end (int): gene end.

    Returns:
        np.array: relative distances.
    """
    return np.where(points < start, start - points, np.where(points > end, points - end, 0))


def calculate_rel_dist_from_centers(cis_bins_df,
                                    gene_start,
                                    gene_end):
    """Finds relative distance of bins centers to gene with [start; end] coordinates.

    Args:
        cis_bins_df (pd.DataFrame): genomic bins on the same chromosome as gene.
        start (int): gene start.
        end (int): gene end.

    Returns:
        np.array: relative distances.
    """
    bins_centers = (cis_bins_df['start'] + cis_bins_df['end']) // 2
    rel_dists = make_rel_dist_vector(bins_centers, gene_start, gene_end)
    return rel_dists


def make_subtrack(bins_df,
                  centers_df,
                  bg_track,
                  impute=False,
                  ifactor=0,
                  bg_count_name='score'):
    """Computes signal and bg track in given bins.

    Args:
        bins_df (pd.DataFrame): a bed3 df with bins.
        centers_df (pd.DataFrame): a contacts centers df.
        bg_track (pd.DataFrame): a bed4 df with binned bg counts
            ('chrom', 'start', 'end', 'count').
        impute (bool): whether to impute missing bg values or not
            (default: False).
        ifactor (float): an imputation value as a bg count per nt
            (default: None).
        bg_count_name (str): 4th column name in bg_track (default: 'count').

    Returns:
        pd.DataFrame: a df with 7 columns: 'chrom', 'start', 'end',
            'count', 'signal_prob', 'bg_count', 'bg_prob'.

    Raises:
        ValueError: in case `impute` is True and `ifactor` is None.
    """
    bins_coverage = calculate_bins_coverage(bins_df, centers_df)
    bins_coverage['signal_prob'] = bins_coverage['count'] / bins_coverage['count'].sum()
    overlap_with_bg = bf.overlap(bins_coverage,
                                 bg_track,
                                 how="left",
                                 return_input=True,
                                 keep_order=True,
                                 suffixes=('', '_bg'),
                                 return_index=True)\
                        .groupby('index')\
                        .agg({'start': 'min',
                              'end': 'max',
                              'start_bg': 'min',
                              'end_bg': 'max',
                              f'{bg_count_name}_bg': 'sum'})
    bin_sizes = bins_coverage['end'] - bins_coverage['start']
    bg_bin_sizes = overlap_with_bg['end_bg'].astype('int64') - overlap_with_bg['start_bg'].astype('int64')
    rescaled_bg_counts = overlap_with_bg[f'{bg_count_name}_bg'].astype('int64') / bg_bin_sizes * bin_sizes
    if impute:
        if ifactor is None:
            raise ValueError("ifactor must be a positive float, not None.")
        dropout = (bins_coverage['count'] > 0) & (rescaled_bg_counts == 0)
        imputed_rescaled_bg_counts = np.where(dropout,
                                              bin_sizes * ifactor,
                                              rescaled_bg_counts)
        bins_coverage['impute'] = dropout
        bins_coverage['bg_count'] = imputed_rescaled_bg_counts
    else:
        bins_coverage['bg_count'] = rescaled_bg_counts
        bins_coverage['impute'] = False
    bins_coverage['bg_prob'] = bins_coverage['bg_count'] / bins_coverage['bg_count'].sum()
    return bins_coverage


compute_track = make_subtrack


def make_track(dna_parts: pd.DataFrame,
               rna_name: str,
               chrom_dict: Dict[str, int],
               selection_dict: Dict,
               annot_dict: Dict,
               bg_track: pd.DataFrame,
               impute: bool = False,
               ifactor: float = 0) -> pd.DataFrame:
    rna_cis_factor = selection_dict[rna_name]['cis_factor']
    rna_trans_size = selection_dict[rna_name]['trans_bin_size']
    cis_start_size = 5000
    rna_chrom = annot_dict[rna_name]['chrom']
    rna_chrom_length = chrom_dict[rna_chrom]
    rna_start = annot_dict[rna_name]['start']
    rna_end = annot_dict[rna_name]['end']

    specific_dna_parts = dna_parts.query('name == @rna_name')
    specific_dna_centers = make_interval_centers(specific_dna_parts)
    trans_bins = make_trans_bins(rna_trans_size, chrom_dict, rna_chrom)
    trans_track = make_subtrack(trans_bins, specific_dna_centers, bg_track, impute=impute, ifactor=ifactor, bg_count_name='score')
    trans_track['raw_bg_prob'] = trans_track['bg_prob']
    trans_track['scaling_factor'] = 1
    cis_bins = make_cis_bins(rna_cis_factor, cis_start_size, rna_chrom, rna_chrom_length, rna_start, rna_end, rna_trans_size)
    cis_track = make_subtrack(cis_bins, specific_dna_centers, bg_track, impute=impute, ifactor=ifactor, bg_count_name='score')
    cis_track['raw_bg_prob'] = cis_track['bg_prob']
    cis_track['scaling_factor'] = 1
    total_track = bf.sort_bedframe(pd.concat((trans_track, cis_track), ignore_index=True))
    return total_track
