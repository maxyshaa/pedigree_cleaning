"""This module contains functions for merging 2 datasets together and modifying columns"""

import pandas as pd
from preprocessing.utils import clear_string_val, change_sex

def merge_1stdataframes(ped_df: pd.DataFrame, pedid_match: pd.DataFrame, geno_id: pd.DataFrame) -> pd.DataFrame:
    """
    Merge filtered ped_df, pedid_match, and filtered geno_id into a single dataframe
    """

    geno_id = geno_id.copy()
    geno_id["sex"] = geno_id["sex"].map({"Female": 2, "Male": 1}).fillna(0).astype(int).astype(str)

    merged_df = (ped_df \
        .merge(pedid_match, left_on="id", right_on="horse_id", how="left") \
        .drop(columns=["horse_id"]) \
        .rename(columns={"Equinome ID": "equinomeID"}) \
        .merge(geno_id, on="equinomeID", how="left") \
        .rename(columns={"id_x": "horse_id", "id_y": "bed_id"})) \
        

    # mismatch_condition = (merged_df["sex_x"] != merged_df["sex_y"]) & \
    #                      ~merged_df["sex_x"].isna() & \
    #                      ~merged_df["sex_y"].isna() & \
    #                      ~merged_df["equinomeID"].isna()
    
    # if mismatch_condition.any():
    #     raise ValueError("Mismatch found between 'sex_x' and 'sex_y' columns for those who are not genotyped")
    #     # print(merged_df[mismatch_condition])

    # changing columns format
    merged_df = merged_df.drop(columns=["sex_y", "Horse Name", "Year of Birth"])
    merged_df["genotyped"] = ~merged_df["bed_id"].isna()
    merged_df = merged_df.rename(columns={"name":"horse_name", "sex_x":"sex", "equinomeID":"equinome_id",
                                          "batchID":"batch", "SNPChip": "snp_chip", 
                                          "Country Reported":"country_reported"})
    # print(f"The merged sheets are now as a 1st pedigree dataset with the size {merged_df.shape}")
    return merged_df

# there are 82 inidividuals who have sex mismatch but it's better to verify through genotype. For now I will use sex info from pedigree
# there are 127 individuals with years of birth not matching between pedigree and genotype data, will keeping yobs from pedigree


def clear_ped_additional(ped2nd: pd.DataFrame) -> pd.DataFrame:
    """This function changes the format of the additional dataset of pedigree to be able to concatenate with the base one"""
    ped2nd_cleaned = ped2nd.copy()
    ped2nd_cleaned = ped2nd_cleaned.drop(columns=["Sire", "Dam"])
    ped2nd_cleaned = ped2nd_cleaned.rename(columns={"id":"bed_id", "batchID":"batch",
                                    "equinomeID":"equinome_id", "SNPChip":"snp_chip",
                                    "Year of Birth":"YOB", "Country Reported":"country_reported",
                                    "Horse Name":"horse_name", "Month of Birth":"MOB",
                                    "Country of Birth":"COB"})
    
    ped2nd_cleaned = ped2nd_cleaned.assign(status="44444", sire_id=None, dam_id=None, genotyped=True)
    ped2nd_cleaned["sex"] = change_sex(ped2nd_cleaned["sex"])
    ped2nd_cleaned["horse_id"] = ped2nd_cleaned["bed_id"].copy()
    
    return ped2nd_cleaned

def modifying_countries(ped_1st: pd.DataFrame, ped2nd: pd.DataFrame, 
                        country_unific_dict: dict) -> pd.DataFrame:
    """"""
    cols = ["country_reported", "COB"]

    ped_1st[cols] = ped_1st[cols].apply(clear_string_val)
    ped2nd[cols] = ped2nd[cols].apply(clear_string_val)

    ped_1st[cols] = ped_1st[cols].replace(country_unific_dict)
    ped2nd[cols] = ped2nd[cols].replace(country_unific_dict)
    
    ped_1st[cols] = ped_1st[cols].apply(lambda x: x.str.upper())
    ped2nd[cols] = ped2nd[cols].apply(lambda x: x.str.upper())

    return ped_1st, ped2nd

import pandas as pd

def concat_peds(ped_1st: pd.DataFrame, ped2nd: pd.DataFrame, col_order) -> pd.DataFrame:
    """Combine all pedigree datasets into one"""
    
    alltogether_pedigree = pd.concat([ped_1st, ped2nd], ignore_index=True)    
    alltogether_pedigree = alltogether_pedigree.reindex(columns=col_order)
    alltogether_pedigree["YOB"] = alltogether_pedigree["YOB"].astype("float").astype("Int64")
    # print(alltogether_pedigree.shape)
    
    return alltogether_pedigree
