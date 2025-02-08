import os
import pandas as pd
from preprocessing.load_data import read_pedigree_sheets, get_pedigree_csv, read_fam
from preprocessing.filter_geno import filter_by_fam, filter_by_chip, update_idmatch
from preprocessing.clean_data import clear_colour, fix_logic, clean_by_name
from preprocessing.match_n_merge import merge_1stdataframes, clear_ped_additional, modifying_countries, concat_peds
from preprocessing.utils import save_file, clear_string_val


geno_path = "data"
results_path = "results"

country_unific_dict = {"australia": "AUS",
                       "aus": "AUS",
                       "argentina": "ARG",
                       "arg": "ARG",
                       "brazil": "BRZ",
                       "brz": "BRZ",
                       "canada": "CAN",
                       "can": "CAN",
                       "great britain": "GB",
                       "germany":"GER",
                       "united kingdown": "GB",
                       "gb": "GB",
                       "ireland": "IRL",
                       "IRELAND": "IRL",
                       "irl": "IRL",
                       "ire":"IRL",
                       "ity":"ITY",
                       "New Zealand": "NZ",
                       "new zealand": "NZ",
                       "nz": "NZ",
                       "South Africa": "SAF",
                       "south africa": "SAF",
                       "saf": "SAF",
                       "France": "FR",
                       "france": "FR",
                       "fr": "FR",
                       "Germany": "GER",
                       "ger": "GER",
                       "Japan": "JAP",
                       "japan": "JAP",
                       "jan": "JAP",
                       "Uruguay": "URU",
                       "uruguay": "URU",
                       "uru": "URU",
                       "United Arab Emirates": "UAE",
                       "united arab emirates": "UAE",
                       "uae": "UAE",
                       "usa": "USA"}

col_order = ["horse_id", "status", "horse_name", "sire_id", 
         "dam_id", "YOB", "MOB", "sex", "colour", "COB",
         "equinome_id", "bed_id", "batch", "snp_chip", 
         "country_reported", "genotyped"]

if os.path.exists(os.path.join(results_path, "bedids2exclude.txt")):
    os.remove(os.path.join(results_path, "bedids2exclude.txt"))
    print("bedid2exclude has been removed")

def main() -> None:

    # initializing a list to hold removed bed_ids
    removed_bedids = []

    ### Step 1: Load all data
    # 1st part of pedigree
    pedid_match, ped_df, geno_id = read_pedigree_sheets(
        os.path.join("data/raw_pedigree.xlsx")
        )
    # 2nd part of pedigree
    ped_addit = get_pedigree_csv(
        os.path.join("data/raw_notpedigree.csv")
        )
    # plink fam files
    fam_11k_ids = read_fam(os.path.join(geno_path, "TB_11K.fam"))
    fam_6k_ids = read_fam(os.path.join(geno_path, "TB_6K.fam"))
    bed_ids = fam_11k_ids.union(fam_6k_ids)


    ### Step 2: Clean 1st pedigree DataFrame
    ped_cleaned = fix_logic(clear_colour(init_file=ped_df, col_name="colour"))
    print("----------Cleaning of 1st pedigree has finished----------")


    ### Step 3: Clean genotype information and remove duplicates based on SNPchip
    geno_id_chip_filtered = filter_by_chip(
        filter_by_fam(geno_id, "id", bed_ids), 
        snp_col="SNPChip", batch_col="batchID", equinomeid="equinomeID", bedid="id", 
        removed_bedids=removed_bedids)
    print("----------Filtering on presence in fam and prioritizing on chip info of genotypes has finished; geno_id has updated for 1st pedigree----------")

    ped_addit_chip_filtered = filter_by_chip(
        filter_by_fam(ped_addit, "id", bed_ids), 
        snp_col="SNPChip", batch_col="batchID", equinomeid="equinomeID", bedid="id", 
        removed_bedids=removed_bedids)
    print("----------Filtering on presence in fam and prioritizing on chip info of genotypes has finished; geno_id has updated for 2nd pedigree----------")


    # Step 4: Update pedid_match with the cleaned geno_id
    pedid_match_chip_filtered = update_idmatch(pedid_match=pedid_match, 
                                               new_geno=geno_id_chip_filtered)
    print("----------pedid_match has updated: no extra SNPchips----------")


    ### Step 5: Handle duplicated horse_id entries with different EquinomeID [only for the first ped_df]

    # finding duplicates by horse_id and getting their equinome id
    duplicated_horse_ids = pedid_match_chip_filtered[pedid_match_chip_filtered["horse_id"].duplicated(keep=False)]
    dup_merged = duplicated_horse_ids.merge(geno_id_chip_filtered, left_on="Equinome ID", 
                                            right_on="equinomeID", how="left")
    
    # filtering those equinome id by chip priority
    no_eqdup = filter_by_chip(dup_merged, "SNPChip", "batchID", "horse_id", "id", removed_bedids)
    removed_ids = dup_merged[~dup_merged["Equinome ID"].isin(no_eqdup["Equinome ID"])]["Equinome ID"].tolist()
    print(f"removing these equinome ids {removed_ids}")
    
    # update new_geno
    geno_id_nodup = geno_id_chip_filtered[~geno_id_chip_filtered["equinomeID"].isin(removed_ids)]
    print("----------geno_id has updated: no extra equinomeID per horse----------")

    # update pedid_match
    pedid_match_nodup = update_idmatch(pedid_match=pedid_match_chip_filtered, new_geno=geno_id_nodup)
    print("----------pedid_match has updated: no extra equinomeID per horse----------")
    save_file(pedid_match_nodup, "results/", "updates_pedid.csv")

    # Step 6: Modifying columns and prepare pedigres
    ped_1stmerged = merge_1stdataframes(ped_df=ped_cleaned, 
                                         pedid_match=pedid_match_nodup, 
                                         geno_id=geno_id_nodup)
    ped_addit_prep = clear_ped_additional(ped_addit_chip_filtered)
    
    ped_1stmerged, ped_addit_prep = modifying_countries(ped_1stmerged, ped_addit_prep,
                                       country_unific_dict)
    
    # Step 7: Concatenating 2 datasets
    final_pedigree = concat_peds(ped_1stmerged, ped_addit_prep, col_order)
  
    # saving results
    save_file(final_pedigree, "results/", "combined_pedigree.csv") # intermediate datatable
    print(f"The merged sheets are now with the size {final_pedigree.shape}")

    # make only meaningfull entries for pedigree
    only_pedigree = final_pedigree[~final_pedigree["horse_name"].str.lower().str.contains("unnamed", na=False) & 
                                   final_pedigree["horse_name"].notna() &
                                   (final_pedigree["horse_name"] != "")
                                   ].copy()
    
    # remove potential duplicates by name
    only_pedigree["horse_name"] = only_pedigree["horse_name"].str.strip()    
    only_pedigree = clean_by_name(only_pedigree)
    save_file(only_pedigree, "results/", "clean_pedigree.csv")
    print(f"The only pedigree info has size: {only_pedigree.shape}")

    # Write all removed bedids to a file
    print(f"Removing less informative bed ids, total amount: {len(removed_bedids)}")
    with open("results/bedids2exclude.txt", 'w') as f:
        for bedid in removed_bedids:
            f.write(f"{bedid} {bedid}\n")

if __name__ == "__main__":
    main()