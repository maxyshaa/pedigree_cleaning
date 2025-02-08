"""This module contains functions to clean the data and bringing to the same format 2 files"""

import pandas as pd
import numpy as np


def clear_colour(init_file: pd.DataFrame, col_name: str) -> pd.DataFrame:
    """Correct the colour variable in the ped_df.
    Args:
        init_file: dataframe to be cleared
        col_name: name of the colour column in the file
    Returns:
        df: returns dataframe with removed unclear records
    """
    incor_color = ~init_file[col_name].isna() & init_file[col_name].str.contains(r'\d', regex=True)
    # modified_rows = init_file.loc[incor_color]
    modified_ids = init_file.loc[incor_color, "id"]
    init_file.loc[incor_color, col_name] = np.nan
    
    if not modified_ids.empty:
        print(f"The following IDs were modified in the column '{col_name}':\n{modified_ids.tolist()}")
    else:
        print(f"No rows were modified in the column '{col_name}'.")
    
    return init_file 


def fix_logic(init_file: pd.DataFrame) -> pd.DataFrame:
    """Check the pedigree in terms of logical errors and removes incorrect information
    Args:
        init_file: dataframe to be cleared
    Returns:
        df: returns fixed dataframe
    """
    # if some inidividuals have multiple sires or dams
    multipar = init_file \
                .groupby(["id"], as_index=False) \
                .agg({"sire_id": pd.Series.nunique, "dam_id":pd.Series.nunique}) \
                .query("sire_id > 1 | dam_id > 1")
    
    if len(multipar) > 0:
        error_message = "DataFrame contains rows with multiple sire/dam of an inidividual"
        raise ValueError(error_message) # since we don't have this issue so I just raise the error here
    else:
        print("Multiple parents are not detected")
    
    # if sex code matches with partnership; 1 codes a male
    sex_1dam = init_file.query("(sex=='1') and (id in dam_id)")
    sex_2sire = init_file.query("(sex=='2') and (id in sire_id)")

    if len(sex_1dam) or len(sex_2sire) > 0:
        error_message = "DataFrame contains rows with wrong sex assigned"
        raise ValueError(error_message)
    else:
        print("Sex conflict is not detected")

    # if some individuals have parents younger than themselves
    slice_cols = ["id", "sire_id", "dam_id", "YOB"]
    wrongage_d = init_file[slice_cols] \
                    .merge(init_file[slice_cols], left_on="dam_id", right_on="id", suffixes=("_child", "_dam")) \
                    .query("YOB_child < YOB_dam")
    
    wrongage_s = init_file[slice_cols] \
                    .merge(init_file[slice_cols], left_on="sire_id", right_on="id", suffixes=("_child", "_sire")) \
                    .query("YOB_child < YOB_sire")
    
    if len(wrongage_d) > 0:
        wrong_ids = wrongage_d["id_child"].tolist()
        init_file.loc[init_file["id"].isin(wrong_ids), "dam_id"] = np.nan
        print(f"")
        print(f"dam_id for individual {init_file.loc[init_file['id'].isin(wrong_ids), 'id'].to_string(index=False)} is yonger than the daughter hence dam_id record for the individuals is removed")
    elif len(wrongage_s) > 0:
        wrong_ids = wrongage_s["id_child"].tolist()
        init_file.loc[init_file["id"].isin(wrong_ids), "sire_id"] = np.nan
        print(f"sire_id for individual {init_file.loc[init_file['id'].isin(wrong_ids),'id'].to_string(index=False)} is removed")
    else:
        print("Wrong parent-offspring relation based on age is not detected")
        pass
    
    return init_file

def dup_indices(i, ped_df):
    """
    Find indices of duplicates for a given row index `i`.
    """
    return ped_df.index[(ped_df['horse_name'].str.lower() == ped_df.loc[i, 'horse_name'].lower()) &
                        (abs(ped_df['YOB'] - ped_df.loc[i, 'YOB']) < 2)].tolist()

def remove_indices(i, dups, ped_df):
    """
    Identify which duplicate index to remove based on the least informative row.
    """
    dup_index = dups[i]
    dup_info = [ped_df.loc[idx].astype(str).str.len().sum() for idx in dup_index]
    to_remove = dup_index[np.argmin(dup_info)]
    print(f"{to_remove} were removed as a duplicate")
    return to_remove

def clean_by_name(ped_df):
    """
    Remove duplicate rows based on horse_name and YOB.
    """
    # Generate duplicate indices for each row
    dups = {i: dup_indices(i, ped_df) for i in ped_df.index}

    # Identify indices with duplicates to process
    indices_with_dups = [i for i, dup in dups.items() if len(dup) > 1]

    # Find indices to remove
    indices_remove = set(remove_indices(i, dups, ped_df) for i in indices_with_dups)

    # Drop rows from DataFrame
    ped_df = ped_df.drop(indices_remove)

    print(f"After dups by name and DOB: {ped_df.shape}")
    return ped_df
