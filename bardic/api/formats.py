"""
Describes two types of data used in BaRDIC and corresponding file types:
1. DnaDataset and corresponding .dnah5 file type.
2. Rdc and corresponding .rdc file type.
"""


from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import bioframe as bf
import h5py
import numpy as np
import pandas as pd

from .schemas import GeneCoord, RnaAttrs, RnaPixelRecord, SplineResult


class StatusProperty:
    """
    A descriptor class that represents a status property.

    This class allows accessing and setting a boolean status property
    stored in an HDF5 file. The property is lazily loaded from the file
    when accessed for the first time.

    Attributes
    ----------
    private_name : str
        The name of the private attribute where the property value is stored.
    public_name : str
        The name of the public property.

    Methods
    -------
    __set_name__(self, owner, name: str) -> None:
        Sets the private and public names of the property.
    __get__(self, obj, objtype=None) -> bool:
        Retrieves the value of the property.
    __set__(self, obj, value: bool) -> None:
        Sets the value of the property.

    Usage
    -----
    Define a class attribute using the `StatusProperty` descriptor
    to represent a status property.

    Example
    -------
    class MyClass:
        status : StatusProperty
            Represents the status of something.
        is_active : StatusProperty
            Represents the active status of something.
    """

    def __set_name__(self, owner, name: str) -> None:
        """
        Set the name of the descriptor.

        Parameters
        ----------
        owner : type
            The class that owns the descriptor.
        name : str
            The name of the descriptor.

        Returns
        -------
        None
            This method does not return anything.
        """
        self.private_name = '_' + name
        self.public_name = name

    def __get__(self, obj, objtype=None) -> bool:
        """
        Retrieve the value of the attribute from the object.

        If the attribute value is `None`, it retrieves the value from the HDF5 file associated with the object.
        The attribute value is expected to be a boolean, and if it is not, an exception is raised.

        Parameters
        ----------
        obj : object
            The object on which the attribute is accessed.
        objtype : type, optional
            The type of the object.

        Returns
        -------
        value : bool
            The value of the attribute.

        Raises
        ------
        ValueError
            If the attribute value is not a boolean.

        Examples
        --------
        >>> obj = MyClass()
        >>> obj.my_attribute
        True
        """
        value = getattr(obj, self.private_name)
        if value is None:
            with h5py.File(obj.fname, 'r') as f:
                value = f.attrs[self.public_name]
                if value in (True, False):
                    setattr(obj, self.private_name, value)
                else:
                    raise ValueError("Invalid value. The value must be a boolean.")
        return value

    def __set__(self, obj, value: bool) -> None:
        """
        Set the value of the attribute.

        Parameters
        ----------
        obj : object
            The object on which the attribute is being set.
        value : bool
            The value to be set.

        Raises
        ------
        ValueError
            If the value is not a boolean.

        Returns
        -------
        None
        """
        if value not in (True, False):
            raise ValueError("Invalid value. The value must be a boolean.")
        with h5py.File(obj.fname, 'a') as f:
            f.attrs[self.public_name] = value
        setattr(obj, self.private_name, value)


class VersionProperty:
    """
    A descriptor class that represents a version property.

    This class allows accessing and setting a version property of an object.
    The version property is stored in an HDF5 file and can only have values
    from a predefined list of supported versions.

    Parameters:
    -----------
    supported_versions : Tuple[str]
        A tuple of strings representing the supported versions.

    Attributes:
    -----------
    supported_versions : Tuple[str]
        A tuple of strings representing the supported versions.

    private_name : str
        The private name of the version property.

    public_name : str
        The public name of the version property.
    """

    def __init__(self, supported_versions: Tuple[str]) -> None:
        """
        Initialize the Formats class.

        Parameters
        ----------
        supported_versions : Tuple[str]
            A tuple of supported versions.

        Returns
        -------
        None
        """
        self.supported_versions = supported_versions

    def __set_name__(self, owner, name: str) -> None:
        """
        Set the name of the version property.

        This method is automatically called by the Python interpreter when
        the descriptor is assigned to a class attribute.

        Parameters:
        -----------
        owner : type
            The owner class of the version property.

        name : str
            The name of the version property.
        """
        self.private_name = '_' + name
        self.public_name = name

    def __get__(self, obj, objtype=None) -> str:
        """
        Get the value of the version property.

        This method is automatically called when the version property is accessed.

        Parameters:
        -----------
        obj : object
            The object instance that the version property belongs to.

        objtype : type, optional
            The type of the object instance.

        Returns:
        --------
        value : str
            The value of the version property.

        Raises:
        -------
        ValueError
            If the version property is not found in the HDF5 file or if the value
            is not in the list of supported versions.
        """
        value = getattr(obj, self.private_name)
        if value is None:
            with h5py.File(obj.fname, 'r') as f:
                value = f.attrs[self.public_name]
                if value in self.supported_versions:
                    setattr(obj, self.private_name, value)
                else:
                    raise ValueError(f"Unsupported version: {value}. Supported versions are {', '.join(self.supported_versions)}")
        return value

    def __set__(self, obj, value: str) -> None:
        """
        Set the value of the version property.

        This method is automatically called when the version property is set.

        Parameters:
        -----------
        obj : object
            The object instance that the version property belongs to.

        value : str
            The value to set for the version property.

        Raises:
        -------
        ValueError
            If the value is not in the list of supported versions.
        """
        if value not in self.supported_versions:
            raise ValueError(f"Unsupported version: {value}. Supported versions are {', '.join(self.supported_versions)}")
        with h5py.File(obj.fname, 'a') as f:
            f.attrs[self.public_name] = value
        setattr(obj, self.private_name, value)


class DnaDataset:
    """
    Handles DNA parts of contacts that are stored in an HDF5 .dnah5 file.

    Parameters
    ----------
    fname : str
        Filename of the corresponding .dnah5 file.
        If the file doesn't exist, it will be created.
    chromsizes : Optional[Dict[str, int]], default None
        Dictionary in the form of `{chromosome_name: chromosome_size}`,
        that contains chromosome sizes.
    annotation : Optional[Dict[str, GeneCoord]], default None
        Dictionary in the form of `{gene_name: gene_coordinates}`,
        that contains gene coordinates.

    Attributes
    ----------
    chrom_groupname : str
        Name of the chromosome sizes group in the .dnah5 file.
    dna_groupname : str
        Name of the dna parts group in the .dnah5 file.
    are_binsizes_selected : bool
        Indicates whether bin sizes were selected.
    version : str
        Version of .dnah5 schema.

    Raises
    ------
    ValueError
        If the version of the existing .dnah5 file is not `"1"`.

    """

    supported_versions = ("1", )
    chrom_groupname: str = "chrom_sizes"
    dna_groupname: str = "dna_parts"

    are_binsizes_selected = StatusProperty()
    version = VersionProperty(supported_versions)

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None,
                 annotation: Optional[Dict[str, GeneCoord]] = None) -> None:
        """
        Constructor

        Parameters
        ----------
        fname : str
            Filename of the corresponding .dnah5 file.
            If the file doesn't exist, it will be created.
        chromsizes : Optional[Dict[str, int]], default None
            Dictionary in the form of `{chromosome_name: chromosome_size}`,
            that contains chromosome sizes.
        annotation : Optional[Dict[str, GeneCoord]], default None
            Dictionary in the form of `{gene_name: gene_coordinates}`,
            that contains gene coordinates.

        Raises
        ------
        ValueError
            If the version of the existing .dnah5 file is not `"1"`.

        """
        self.fname: Path = Path(fname)
        if not self.fname.exists():
            with h5py.File(self.fname, 'w'):
                pass
            self.are_binsizes_selected = False
            self.version = "1"

        self._are_binsizes_selected: Optional[bool] = None
        self._version: Optional[str] = None
        if self.version not in self.supported_versions:
            raise ValueError(f"Unsupported version: {self.version}. Supported versions are {', '.join(self.supported_versions)}")

        self._chromsizes: Optional[Dict[str, int]] = chromsizes
        if chromsizes is not None:
            self._write_chromsizes()
        self._annotation: Optional[Dict[str, GeneCoord]] = annotation

    @property
    def chromsizes(self) -> Dict[str, int]:
        """
        Gets chromosome sizes from the .dnah5 file.

        Returns
        -------
        dict
            A dictionary of chromosome sizes in the form
            of `{chromosome_name: chromosome_size}`.

        Raises
        ------
        Exception
            If there are no chromosome sizes in the file.

        """
        if self._chromsizes is None:
            try:
                self._chromsizes = self._read_chromsizes()
            except FileNotFoundError:
                raise Exception("Failed to read chromosome sizes from the .dnah5 file: File not found.")
            except Exception as e:
                raise Exception("Failed to read chromosome sizes from the .dnah5 file: " + str(e))
        return self._chromsizes

    @chromsizes.setter
    def chromsizes(self, chromsizes_dict: Dict[str, int]) -> None:
        """
        Writes new chromosome sizes to the .dnah5 file.

        Parameters
        ----------
        chromsizes_dict : dict
            A dictionary of chromosome sizes in the form
            of `{chromosome_name: chromosome_size}`.

        """
        self._chromsizes = chromsizes_dict
        self._write_chromsizes()

    @property
    def annotation(self) -> Dict[str, GeneCoord]:
        """
        Gets gene annotation from the .dnah5 file.

        Returns
        -------
        dict
            Dictionary in the form of `{gene_name: gene_coordinates}`,
            that contains gene coordinates.

        Raises
        ------
        Exception
            If there is no gene annotation in the file.

        """
        if self._annotation is None:
            try:
                self._annotation = self._get_annotation()
            except Exception:
                raise Exception("Failed to get gene annotation.")
        return self._annotation

    @annotation.setter
    def annotation(self, annotation_dict: Dict[str, GeneCoord]) -> None:
        """
        Writes new gene annotation in the .dnah5 file.

        Parameters
        ----------
        annotation_dict : dict
            Dictionary in the form of `{gene_name: gene_coordinates}`,
            that contains gene coordinates.
        """
        self._annotation = annotation_dict
    

    def _read_chromsizes(self) -> Dict[str, int]:
        """
        Read chromosome sizes from the file.

        Returns
        -------
        Dict[str, int]
            A dictionary mapping chromosome names to their corresponding sizes.

        Raises
        ------
        Exception
            If the file does not exist, or if the required groups and datasets are not found.
        """
        if not self.fname.exists():
            raise Exception

        with h5py.File(self.fname, 'r') as f:
            if self.chrom_groupname not in f:
                raise Exception
            chromsizes_group = f[self.chrom_groupname]
            if 'chrom' not in chromsizes_group:
                raise Exception
            names = chromsizes_group['chrom'].asstr()[()]
            if 'size' not in chromsizes_group:
                raise Exception
            sizes = chromsizes_group['size'][()]

        return dict(zip(names, sizes))

    def _write_chromsizes(self) -> None:
        """
        Write chromosome sizes to an HDF5 file.

        This method creates a group in the HDF5 file and stores the chromosome names and sizes as datasets within that group.

        Returns
        -------
        None
        """
        chromsizes_dict = self.chromsizes
        with h5py.File(self.fname, 'a') as f:
            names = np.array(list(chromsizes_dict.keys()), dtype='O')
            sizes = np.array(list(chromsizes_dict.values()), dtype='int64')
            if self.chrom_groupname in f:
                del f[self.chrom_groupname]
            chromsizes_group = f.create_group(self.chrom_groupname)
            chromsizes_group.create_dataset('chrom', data=names)
            chromsizes_group.create_dataset('size', data=sizes)

    def _read_dna_parts(self, f: h5py.File, rna_name: str) -> pd.DataFrame:
        """
        Read DNA parts from an HDF5 file.

        Parameters
        ----------
        f : h5py.File
            The HDF5 file object.
        rna_name : str
            The name of the RNA.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the DNA parts.

        Raises
        ------
        Exception
            If the DNA group or the specified RNA name is not found in the file.
        """
        if self.dna_groupname not in f:
            raise Exception
        dna_parts_group = f[self.dna_groupname]
        if rna_name not in dna_parts_group:
            raise Exception
        rna_group = dna_parts_group[rna_name]
        dna_parts_list = list()
        for chrom_name, chrom_group in rna_group.items():
            starts = chrom_group['start'][()]
            ends = chrom_group['end'][()]
            dna_parts_list.append(pd.DataFrame({'chrom': chrom_name, 'start': starts, 'end': ends}))
        return pd.concat(dna_parts_list, ignore_index=True)

    def read_dna_parts(self, rna_name: str) -> pd.DataFrame:
        """
        Read DNA parts for a given RNA name.

        Parameters
        ----------
        rna_name : str
            The name of the RNA.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the DNA parts.
        """
        with h5py.File(self.fname, 'r') as f:
            dna_parts = self._read_dna_parts(f, rna_name)
        return dna_parts

    def _write_dna_parts(self,
                         f: h5py.File,
                         rna_name: str,
                         dna_parts: pd.DataFrame) -> None:
        """
        Write DNA parts of a single RNA to an HDF5 file.

        Parameters
        ----------
        f : h5py.File
            The HDF5 file object.
        rna_name : str
            The name of the RNA.
        dna_parts : pd.DataFrame
            The DNA parts data as a pandas DataFrame.

        Returns
        -------
        None
        """
        annotation_dict = self.annotation
        rna_annot = annotation_dict[rna_name]

        if self.dna_groupname not in f:
            dna_parts_group = f.create_group(self.dna_groupname)
        else:
            dna_parts_group = f[self.dna_groupname]

        if rna_name in dna_parts_group:
            del dna_parts_group[rna_name]
        rna_group = dna_parts_group.create_group(rna_name)

        for key, value in asdict(rna_annot).items():
            rna_group.attrs[key] = value
        rna_group.attrs['total_contacts'] = dna_parts.shape[0]
        for chrom_name, chrom_df in dna_parts.groupby('chrom'):
            rna_chrom_group = rna_group.create_group(chrom_name)
            starts = chrom_df['start'].values.astype('int64')
            ends = chrom_df['end'].values.astype('int64')
            rna_chrom_group.create_dataset('start', data=starts)
            rna_chrom_group.create_dataset('end', data=ends)

    def write_dna_parts(self,
                        rna_name: str,
                        dna_parts: pd.DataFrame) -> None:
        """
        Write DNA parts of a single RNA to the HDF5 file.

        Parameters
        ----------
        rna_name : str
            The name of the RNA.
        dna_parts : pd.DataFrame
            The DNA parts data.

        Returns
        -------
        None
        """
        with h5py.File(self.fname, 'a') as f:
            self._write_dna_parts_single(f, rna_name, dna_parts)

    def write_dna_parts_batch(self,
                              dna_frame: pd.DataFrame,
                              rna_col: str = 'name') -> None:
        """
        Write DNA parts of multiple RNAs to the HDF5 file.

        Parameters
        ----------
        dna_frame : pd.DataFrame
            DataFrame containing DNA parts.
        rna_col : str, optional
            Column name in `dna_frame` to group DNA parts by RNA name. Defaults to 'name'.

        Returns
        -------
        None
        """
        with h5py.File(self.fname, 'a') as f:
            for rna_name, dna_parts in dna_frame.groupby(rna_col):
                self._write_dna_parts(f, rna_name, dna_parts)

    def _get_coordinates(self, f: h5py.File, rna_name: str) -> GeneCoord:
        """
        Get the coordinates of a gene from an HDF5 file.

        Parameters
        ----------
        f : h5py.File
            The HDF5 file object.
        rna_name : str
            The name of the gene.

        Returns
        -------
        GeneCoord
            The coordinates of the gene.
        """
        if self.dna_groupname not in f:
            raise Exception
        dna_group = f[self.dna_groupname]
        if rna_name not in dna_group:
            raise Exception
        rna_group = dna_group[rna_name]
        return GeneCoord(chrom=rna_group.attrs['chrom'],
                         start=rna_group.attrs['start'],
                         end=rna_group.attrs['end'])

    def get_coordinates(self, rna_name: str) -> GeneCoord:
        """
        Retrieves the genomic coordinates of a given RNA.

        Parameters
        ----------
        rna_name : str
            The name of the RNA.

        Returns
        -------
        GeneCoord
            The genomic coordinates of the RNA.
        """
        with h5py.File(self.fname, 'r') as f:
            annot = self._get_coordinates(f, rna_name)
        return annot

    def _get_annotation(self) -> Dict[str, GeneCoord]:
        """
        Retrieve the annotation information from the specified file.

        Returns
        -------
        annotation_dict : Dict[str, GeneCoord]
            A dictionary containing the annotation information.
            The keys are RNA names and the values are instances of the GeneCoord class.
        """
        with h5py.File(self.fname, 'r') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            annotation_dict = {rna_name: GeneCoord(chrom=rna_group.attrs['chrom'],
                                                   start=rna_group.attrs['start'],
                                                   end=rna_group.attrs['end'])
                               for rna_name, rna_group in dna_group.items()}
        return annotation_dict

    def get_num_contacts(self) -> Dict[str, int]:
        """
        Get the number of contacts for each RNA molecule in the DNA group.

        Returns
        -------
        Dict[str, int]
            A dictionary where the keys are RNA molecule names and the values are the corresponding number of contacts.
        """
        with h5py.File(self.fname, 'r') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            sizes = {rna_name: rna_group.attrs['total_contacts'] for rna_name, rna_group in dna_group.items()}
        return sizes

    def read_rna_attribute_batch(self, attrname: str) -> Dict[str, Any]:
        """
        Read the specified RNA attribute for all RNAs in the file.

        Parameters
        ----------
        attrname : str
            The name of the RNA attribute to read.

        Returns
        -------
        Dict[str, Any]
            A dictionary mapping RNA names to their corresponding attribute values.
        """
        with h5py.File(self.fname, 'r') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            data = {rna_name: rna_group.attrs.get(attrname) for rna_name, rna_group in dna_group.items()}
        return data

    def write_rna_attribute_batch(self, attrname: str, data: Dict[str, Any]) -> None:
        """
        Write the specified RNA attribute for all RNAs to the HDF5 file.

        Parameters
        ----------
        attrname : str
            The name of the attribute to be written.
        data : Dict[str, Any]
            A dictionary containing RNA names as keys and attribute values as values.

        Returns
        -------
        None
        """
        with h5py.File(self.fname, 'a') as f:
            if self.dna_groupname not in f:
                raise Exception
            dna_group = f[self.dna_groupname]
            for rna_name, value in data.items():
                if rna_name in dna_group:
                    dna_group[rna_name].attrs[attrname] = value


class Rdc:
    """
    Represents an RDC (RNA-DNA contacts) file.

    The Rdc class provides methods to read and write data from/to an RDC file.
    It supports different versions of the RDC file format and provides convenient
    access to the data stored in the file.

    Parameters
    ----------
    fname : str
        The path to the RDC file.
    chromsizes : Optional[Dict[str, int]], optional
        A dictionary mapping chromosome names to their sizes. Defaults to None.

    Attributes
    ----------
    supported_versions : tuple
        A tuple of supported RDC file versions.
    pixels_cols_by_version : dict
        A dictionary mapping RDC file versions to the columns present in the pixels dataset.
    chrom_groupname : str
        The name of the group that stores chromosome sizes.
    bg_groupname : str
        The name of the group that stores background track data.
    pixels_groupname : str
        The name of the group that stores pixel data.
    is_scaling_fitted : StatusProperty
        A property indicating whether scaling has been fitted for the RDC file.
    are_peaks_estimated : StatusProperty
        A property indicating whether peaks have been estimated for the RDC file.
    version : VersionProperty
        The version of the RDC file.

    Raises
    ------
    Exception
        If the RDC file version is not supported.
    """
    supported_versions = ("1", "1.1")  # version 1 is deprecated and will be removed in the future.

    pixels_cols_by_version = {'1': {'start': 'int64',
                                    'end': 'int64',
                                    'signal_count': 'int64',
                                    'signal_prob': 'float',
                                    'impute': 'bool',
                                    'bg_count': 'float',
                                    'bg_prob': 'float',
                                    'raw_bg_prob': 'float',
                                    'scaling_factor': 'float',
                                    'fc': 'float',
                                    'pvalue': 'float',
                                    'qvalue': 'float'},
                              '1.1': {'start': 'int64',
                                      'end': 'int64',
                                      'signal_count': 'int64',
                                      'signal_prob': 'float',
                                      'impute': 'bool',
                                      'bg_count': 'float',
                                      'bg_prob': 'float',
                                      'raw_bg_prob': 'float',
                                      'scaling_factor': 'float',
                                      'fc': 'float',
                                      'pvalue': 'float',
                                      'qvalue_global': 'float',
                                      'qvalue_rna': 'float'}}

    chrom_groupname: str = "chrom_sizes"
    bg_groupname: str = "background"
    pixels_groupname: str = "pixels"

    is_scaling_fitted = StatusProperty()
    are_peaks_estimated = StatusProperty()
    version = VersionProperty(supported_versions)

    def __init__(self,
                 fname: str,
                 chromsizes: Optional[Dict[str, int]] = None) -> None:
        self.fname: Path = Path(fname)
        if not self.fname.exists():
            with h5py.File(self.fname, 'w'):
                pass
            self.is_scaling_fitted = False
            self.are_peaks_estimated = False
            self.version = "1.1"  # we write the newest version

        self._chromsizes: Optional[Dict[str, int]] = chromsizes
        if chromsizes is not None:
            self._write_chromsizes()

        self._annotation: Optional[Dict[str, GeneCoord]] = None
        self._is_scaling_fitted: Optional[bool] = None
        self._are_peaks_estimated: Optional[bool] = None

        self._version: Optional[str] = None
        if self.version not in self.supported_versions:
            raise Exception(f"The RDC file version {self.version} is not supported.")

        self.pixels_cols = self.pixels_cols_by_version[self.version]

    @property
    def chromsizes(self) -> Dict[str, int]:
        if self._chromsizes is None:
            try:
                self._chromsizes = self._read_chromsizes()
            except Exception:
                raise Exception
        return self._chromsizes

    @chromsizes.setter
    def chromsizes(self, chromdict: Dict[str, int]) -> None:
        self._chromsizes = chromdict
        self._write_chromsizes()

    def _read_chromsizes(self) -> Dict[str, int]:
        if not self.fname.exists():
            raise Exception

        with h5py.File(self.fname, 'r') as f:
            if self.chrom_groupname not in f:
                raise Exception
            chromsizes_group = f[self.chrom_groupname]
            if 'chrom' not in chromsizes_group:
                raise Exception
            names = chromsizes_group['chrom'].asstr()[()]
            if 'size' not in chromsizes_group:
                raise Exception
            sizes = chromsizes_group['size'][()]

        return dict(zip(names, sizes))

    def _write_chromsizes(self) -> None:
        chromsizes_dict = self.chromsizes
        with h5py.File(self.fname, 'a') as f:
            names = np.array(list(chromsizes_dict.keys()), dtype='O')
            sizes = np.array(list(chromsizes_dict.values()), dtype='int64')
            if self.chrom_groupname in f:
                del f[self.chrom_groupname]
            chromsizes_group = f.create_group(self.chrom_groupname)
            chromsizes_group.create_dataset('chrom', data=names)
            chromsizes_group.create_dataset('size', data=sizes)

    def read_bg_track(self) -> pd.DataFrame:
        with h5py.File(self.fname, 'r') as f:
            if self.bg_groupname not in f:
                raise Exception
            bg_group = f[self.bg_groupname]
            bg_dfs = list()
            for chrom_name, chrom_group in bg_group.items():
                chrom_df = pd.DataFrame({key: chrom_group[key][()]
                                         for key in ('start', 'end', 'count')})
                chrom_df['chrom'] = chrom_name
                chrom_df = chrom_df[['chrom', 'start', 'end', 'count']]
                bg_dfs.append(chrom_df)
            return pd.concat(bg_dfs, ignore_index=True)

    def write_bg_track(self, bg_track: pd.DataFrame) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.bg_groupname in f:
                del f[self.bg_groupname]
            bg_group = f.create_group(self.bg_groupname)
            for chrom_name, chrom_df in bg_track.groupby('chrom'):
                chrom_group = bg_group.create_group(chrom_name)
                for key in ('start', 'end', 'count'):
                    chrom_group.create_dataset(key, data=chrom_df[key].values)

    def _read_pixels(self,
                     f: h5py.File,
                     rna_name: str,
                     value_fields: Optional[List] = None,
                     chrom_type: Optional[str] = None) -> pd.DataFrame:
        if self.pixels_groupname not in f:
            raise Exception
        pixels_group = f[self.pixels_groupname]

        if rna_name not in pixels_group:
            raise Exception
        rna_group = pixels_group[rna_name]

        if chrom_type is None:
            valid_chroms = list(rna_group.keys())
        elif chrom_type == 'cis':
            valid_chroms = [rna_group.attrs['chrom']]
        elif chrom_type == 'trans':
            valid_chroms = [chrom for chrom in rna_group.keys() if chrom != rna_group.attrs['chrom']]
        else:
            raise ValueError

        mandatory_fields = ['start', 'end']
        all_fields = mandatory_fields
        chrom_dfs = list()

        for chrom_name, chrom_group in rna_group.items():
            if chrom_name not in valid_chroms:
                continue

            if value_fields is None:
                value_fields = [item for item in chrom_group.keys() if item not in mandatory_fields]
            all_fields = mandatory_fields + value_fields

            data = {field_name: field_data[()]
                    for field_name, field_data in chrom_group.items()
                    if field_name in all_fields}
            data['chrom'] = chrom_name
            chrom_dfs.append(pd.DataFrame(data)[['chrom'] + all_fields])

        result = pd.concat(chrom_dfs, ignore_index=True)
        return result

    def read_pixels(self,
                    rna_name: str,
                    value_fields: Optional[List] = None,
                    chrom_type: Optional[str] = None) -> pd.DataFrame:
        with h5py.File(self.fname, 'r') as f:
            return self._read_pixels(f, rna_name, value_fields=value_fields, chrom_type=chrom_type)

    def _write_pixels(self,
                      f: h5py.File,
                      rna_name: str,
                      pixels_df: pd.DataFrame,
                      rna_coords: GeneCoord,
                      rna_attrs: RnaAttrs) -> None:
        if self.pixels_groupname not in f:
            pixels_group = f.create_group(self.pixels_groupname)
        else:
            pixels_group = f[self.pixels_groupname]

        if rna_name in pixels_group:
            del pixels_group[rna_name]
        rna_group = pixels_group.create_group(rna_name)

        rna_dict = asdict(rna_attrs)
        rna_dict.update(asdict(rna_coords))
        for key, value in rna_dict.items():
            if value is not None:
                rna_group.attrs[key] = value

        for chrom_name, chrom_df in pixels_df.groupby('chrom'):
            chrom_group = rna_group.create_group(chrom_name)
            for col_name, col_dtype in self.pixels_cols.items():
                if col_name in chrom_df:
                    col_data = chrom_df[col_name].values.astype(col_dtype)
                    chrom_group.create_dataset(col_name, data=col_data)

    def write_pixels(self, rna_name, pixels_df, rna_coords, rna_attrs) -> None:
        with h5py.File(self.fname, 'a') as f:
            self._write_pixels(f, rna_name, pixels_df, rna_coords, rna_attrs)

    def write_pixels_batch(self, data: Dict[str, RnaPixelRecord]) -> None:
        with h5py.File(self.fname, 'a') as f:
            for rna_name, rna_record in data.items():
                pixels_df = rna_record.pixels
                rna_coord = rna_record.gene_coord
                rna_attrs = rna_record.rna_attrs
                self._write_pixels(f, rna_name, pixels_df, rna_coord, rna_attrs)

    def _write_pixels_column(self, f: h5py.File, rna_name: str, col_name: str, col_frame: pd.DataFrame) -> None:
        if self.pixels_groupname not in f:
            raise Exception
        pixels_group = f[self.pixels_groupname]
        if rna_name not in pixels_group:
            raise Exception
        rna_group = pixels_group[rna_name]

        for chrom_name, chrom_df in col_frame.groupby('chrom'):
            if chrom_name not in rna_group:
                raise Exception
            chrom_group = rna_group[chrom_name]
            if col_name not in self.pixels_cols:
                raise Exception
            arr_dtype = self.pixels_cols[col_name]
            arr_data = chrom_df[col_name].values.astype(arr_dtype)
            if col_name in chrom_group:
                chrom_group[col_name][...] = arr_data
            else:
                chrom_group.create_dataset(col_name, data=arr_data)

    def write_pixels_column(self, rna_name: str, col_name: str, col_frame: pd.DataFrame) -> None:
        with h5py.File(self.fname, 'a') as f:
            self._write_pixels_column(f, rna_name, col_name, col_frame)

    def write_pixels_column_batch(self, col_name: str, col_frames: Dict[str, pd.DataFrame]) -> None:
        with h5py.File(self.fname, 'a') as f:
            for rna_name, col_frame in col_frames.items():
                self._write_pixels_column(f, rna_name, col_name, col_frame)

    def read_rna_attribute_batch(self, attrname) -> Dict[str, Any]:
        with h5py.File(self.fname, 'r') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            data = {rna_name: rna_group.attrs.get(attrname)
                    for rna_name, rna_group in pixels_group.items()}
        return data

    def write_rna_attribute_batch(self, attrname: str, data: Dict[str, Any]) -> None:
        with h5py.File(self.fname, 'a') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            for rna_name, value in data.items():
                if rna_name in pixels_group:
                    pixels_group[rna_name].attrs[attrname] = value

    def read_scaling(self, rna_name: str) -> SplineResult:
        if not self.is_scaling_fitted:
            raise Exception
        with h5py.File(self.fname, 'r') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            if rna_name not in pixels_group:
                raise Exception
            rna_group = pixels_group[rna_name]
            spline_result = SplineResult(t=rna_group.attrs['scaling_spline_t'],
                                         c=rna_group.attrs['scaling_spline_c'],
                                         k=rna_group.attrs['scaling_spline_k'])
        return spline_result

    def read_scaling_batch(self) -> Dict[str, SplineResult]:
        rna_spline_t = self.read_rna_attribute_batch('scaling_spline_t')
        rna_spline_c = self.read_rna_attribute_batch('scaling_spline_c')
        rna_spline_k = self.read_rna_attribute_batch('scaling_spline_k')
        scaling_splines = {rna_name: SplineResult(t=rna_spline_t[rna_name],
                                                  c=rna_spline_c[rna_name],
                                                  k=rna_spline_k[rna_name])
                           for rna_name in rna_spline_t}
        return scaling_splines

    def write_scaling_batch(self, scaling_splines: Dict[str, SplineResult]) -> None:
        rna_spline_t = {rna_name: tck.t for rna_name, tck in scaling_splines.items()}
        rna_spline_c = {rna_name: tck.c for rna_name, tck in scaling_splines.items()}
        rna_spline_k = {rna_name: tck.k for rna_name, tck in scaling_splines.items()}
        self.write_rna_attribute_batch('scaling_spline_t', rna_spline_t)
        self.write_rna_attribute_batch('scaling_spline_c', rna_spline_c)
        self.write_rna_attribute_batch('scaling_spline_k', rna_spline_k)

    def _get_annotation(self) -> Dict[str, GeneCoord]:
        with h5py.File(self.fname, 'r') as f:
            if self.pixels_groupname not in f:
                raise Exception
            pixels_group = f[self.pixels_groupname]
            annotation_dict = {rna_name: GeneCoord(chrom=rna_group.attrs['chrom'],
                                                   start=rna_group.attrs['start'],
                                                   end=rna_group.attrs['end'])
                               for rna_name, rna_group in pixels_group.items()}
        return annotation_dict

    @property
    def annotation(self) -> Dict[str, GeneCoord]:
        """
        Returns a dictionary of gene coordinates for the current sequence.

        If the annotation has not been loaded yet, it will be loaded from the appropriate file.

        Returns
        -------
        Dict[str, GeneCoord]
            A dictionary of gene coordinates for the current sequence.
        """
        if self._annotation is None:
            self._annotation = self._get_annotation()
        return self._annotation
