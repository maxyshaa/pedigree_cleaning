"""This module contains functions for reading initial pedigree files and accompanied genotypes"""

import pandas as pd


def read_pedigree_sheets(name: str) -> pd.DataFrame:
    """"""

    sheet_to_df_map = pd.read_excel(name, 
                                sheet_name=None, 
                                dtype=str,
                                na_values=[" ", "", "None"])
    
    pedid_match = sheet_to_df_map["PedIDMatch"]
    ped_df = sheet_to_df_map["PedNew"]
    geno_id = sheet_to_df_map["GenotypeIDs"]

    return pedid_match, ped_df, geno_id


def get_pedigree_csv(name: str) -> pd.DataFrame:
    """"""
    ped_addit = pd.read_csv(name, 
                        dtype=str, na_values=[" ", "", "None"], 
                        usecols=range(12)).query("id.notna().values") # to read only those rows that has an actual records

    return ped_addit


def read_fam(name: str) -> set:
    """
    Read a .fam file and extract a set of IDs
    """
    fam_df = pd.read_csv(name, sep="\s+", header=None)
    ids = set(fam_df[1].unique())
    return ids

    