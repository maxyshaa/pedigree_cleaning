"""This module contains functions for clearing genotype ids"""

import pandas as pd


def filter_by_fam(geno_ids: pd.DataFrame, id_name: str, fam_ids: set) -> pd.DataFrame:
    """
    Filter pedid_match to keep only equinomeIDs that are present in the bed file
    """
    filtered = geno_ids[geno_ids[id_name].isin(fam_ids)]
    abcent = geno_ids[~geno_ids[id_name].isin(fam_ids)]["id"]
    print(f"The number of ID of pedidIDmach that are not presented in fam files {len(abcent)}")
    return filtered

def filter_by_chip(init_file: pd.DataFrame, snp_col:str, batch_col:str, equinomeid:str, bedid:str, 
                   removed_bedids: list) -> pd.DataFrame:
    """
    Filters genotype records based on SNPChip priority and batchID, and stores removed bedids.

    Parameters:
    -----------
    init_file : pd.DataFrame
        The input DataFrame containing genotype records.
    snp_col : str
        Column name indicating the SNPChip information.
    batch_col : str
        Column name for the batch information.
    equinomeid : str
        Column name for equinome ID to drop duplicates based on.
    bedid : str
        Column name for the bed ID.
    removed_bedids : list
        A list to store removed bed IDs from the input DataFrame.

    Returns:
    --------
    pd.DataFrame
        The filtered DataFrame after removing duplicates based on SNPChip priority and batchID.
    """
    init_file = init_file.copy()
    snp_priority = {"SNP670":3, "SNP70_V2":2, "SNP70_PVL":2, "SNP70":2, "SNP50":1}
    init_file.loc[:, "SNPChip_prority"] = init_file[snp_col].map(snp_priority)
    df_sorted = init_file.sort_values(["SNPChip_prority", batch_col], ascending=[False, False])
    # drop duplicates in 'equinomeID' column and keep the first occurrence (highest priority SNPChip and maximum batchID)
    df_filtered = df_sorted.drop_duplicates(subset=equinomeid, keep="first").drop(columns="SNPChip_prority")
    print(f"Duplicated genotypes by chip array were filtered, the shape of filtered dataframe is {df_filtered.shape[0]}")

    duplicates_removed = init_file[~init_file.index.isin(df_filtered.index)]
    exclusive_rows = duplicates_removed[duplicates_removed[bedid].notna()]
    removed_bedids.extend(exclusive_rows[bedid].tolist())
    print(f"Removed {exclusive_rows.shape[0]} rows")
    
    return df_filtered

def update_idmatch(pedid_match: pd.DataFrame, new_geno: pd.DataFrame) -> pd.DataFrame:
    """
    Takes updated lists of EquinomeID of given data and update pedidmatch info for first df.
    """
    if "Equinome ID" not in pedid_match.columns or "equinomeID" not in new_geno.columns:
        raise ValueError("Expected columns 'Equinome ID' in pedid_match and 'equinomeID' in new_geno")
        
    new_pedid_match = pedid_match[pedid_match["Equinome ID"].isin(new_geno["equinomeID"])]

    return new_pedid_match

# I need to extract columns related only to pedigree thing, \
# make columns the same, add info if it's full pedigree or not, if it has genotype or not,\
#  drop duplicates from the additional pedigree
# so the output would be cleaned, combined ped file
# after some other cleaning we can add an additional column that would say \
# what genotype we took for the analysis (like corresponding id_bed or equinome id)