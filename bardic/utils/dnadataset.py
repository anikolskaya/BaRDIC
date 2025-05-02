from typing import Dict, Literal, Optional

import bioframe as bf
import pandas as pd

from ..api.formats import DnaDataset
from ..api.schemas import GeneCoord
from ..api.convert import annotation_to_dict

def validate_dna_frame(dna_frame: pd.DataFrame,
                       chromsizes_dict: Dict[str, int],
                       annotation: Optional[pd.DataFrame] = None,
                       type: Literal['annotation', 'dna_frame'] = 'dna_frame') -> bool:
        """
        Validates a dataframe with DNA parts of contacts.

        Validation checks:
        1. The dataframe is a valid BED dataframe (according to bioframe).
        2. For every RNA:
            2.1. RNA name is present in the annotation dictionary.
            2.2. Every chromosome name is present in the chromosome
                 sizes dictionary.
        3. For every DNA part:
            2.1. End does not exceed the chromosome size.
            2.2. Start is nonnegative.

        Returns
        -------
        bool
            True if all validation checks are successful.
            Otherwise, raises exceptions.

        Raises
        ------
        Exception
            If one of the validation rules is broken.
        """
        print(f'Validating {type}...')
        invalid_rna_names = []
        if not bf.is_bedframe(dna_frame):
            raise Exception
        for rna_name in dna_frame['name'].unique():
            if '/' in rna_name:
                invalid_rna_names.append(rna_name)
                print(f'[Warning]RNA name "{rna_name}" contains invalid characters — replacing with "_"')
                old_rna_name = rna_name
                rna_name = old_rna_name.replace('/', '_')
                dna_frame.loc[dna_frame['name'] == old_rna_name, 'name'] = rna_name

            if (type == 'dna_frame') and (rna_name not in annotation['name'].unique()):
                print(f'[Warning]RNA name "{rna_name}" not found in annotation dictionary — skipping.')
                dna_frame = dna_frame[dna_frame['name'] != rna_name]

        for chrom_name in dna_frame['chrom'].unique():
            if chrom_name not in chromsizes_dict:
                print(f'[Warning]Chromosome name "{chrom_name}" not found in chromosome sizes dictionary — skipping.')
                print(f'RNAs from this chromosome: {dna_frame[dna_frame["chrom"] == chrom_name]["name"].unique()} would be skipped.')
                dna_frame = dna_frame[dna_frame['chrom'] != chrom_name]

        dna_frame = clean_dna_frame(dna_frame, chromsizes_dict)
        return dna_frame


def clean_dna_frame(dna_frame: pd.DataFrame, chromsizes_dict: Dict[str, int]) -> pd.DataFrame:
    mask_end_invalid = dna_frame['end'] > dna_frame['chrom'].map(chromsizes_dict)
    if mask_end_invalid.any():
       bad_chroms = dna_frame.loc[mask_end_invalid, 'chrom'].unique()
       for chrom in bad_chroms:
           print(f'[Warning] End exceeds chromosome size for "{chrom}" — removing these contacts.')
       dna_frame = dna_frame[~mask_end_invalid]

    mask_start_invalid = dna_frame['start'] < 0
    if mask_start_invalid.any():
       bad_chroms = dna_frame.loc[mask_start_invalid, 'chrom'].unique()
       for chrom in bad_chroms:
           print(f'[Warning] Start is negative for "{chrom}" — removing these contacts.')
       dna_frame = dna_frame[~mask_start_invalid]

    return dna_frame

def check_duplicate_rna_names(dna_frame: pd.DataFrame) -> bool:
    if dna_frame['name'].duplicated().any():
        print('[Warning] Duplicate RNA names found in DNA frame — removing duplicates.')
        dna_frame = dna_frame[~dna_frame['name'].duplicated()]
    return dna_frame


def bed2h5(bed_fname: str,
           h5_fname: str,
           chromsizes: Dict[str, int],
           annotation_df: pd.DataFrame) -> DnaDataset:
    dna_frame = bf.read_table(bed_fname, schema='bed6')

    annotation_df = validate_dna_frame(annotation_df, chromsizes, type='annotation')
    annotation_df = check_duplicate_rna_names(annotation_df)

    refined_dna_frame = validate_dna_frame(dna_frame, chromsizes, annotation_df, type='dna_frame')

    annotation = annotation_to_dict(annotation_df)

    present_rnas = refined_dna_frame['name'].unique()
    refined_annotation = {rna_name: gene_coord
                          for rna_name, gene_coord in annotation.items()
                          if rna_name in present_rnas}
    dataset = DnaDataset(h5_fname, chromsizes, refined_annotation)
    dataset.write_dna_parts_batch(refined_dna_frame)
    return dataset
