## **Pedigree Data Preprocessing and Analysis**

This Python module cleans and prepares pedigree datasets for downstream analysis, ensuring a single file is ready for use with R pedigree packages and BLUP models.

### Data overview

When customers request animal genotypes, they must provide the sex, birth month, and country. They may also provide the horse's name, sire, and dam, which form the pedigree.

Potential issues:

1) Due to differences between the northern and southern hemispheres, the birth month might be reversed (e.g., July instead of January), potentially causing year-of-birth discrepancies. Use Equiline.com to verify birth years.
2) Animals in the back pedigree born before 1960, marked as founders, may have multiple records with different IDs. Consider merging these entries.
3) Some genotypes don't have sex information - it's something we could infer using plink.
4) Consider using IBD or KING to perform a quality check to ensure, for instance, that all individuals with the same sire are genetically related.

### Input data

The raw data are located in the `data/ `folder and consist of two files: `raw_pedigree_part1.xlsx` and `raw_pedigree_part2.csv`, representing **Pedigree #1**  and **Pedigree #2**, respectively.

**Pedigree #1**

This is an `.xlsx` file containing individuals with back pedigree data, spread across three sheets:

* **pedIDmatch**

| Equinome ID | Horse Name         | horse_id |
| :---------- | :----------------- | -------: |
| 82          | A Goodlookin Broad |       62 |
| 82_2        | A Goodlookin Broad |       62 |
| CM012       | Sadlers Wells      |       69 |
| CM006       | Fairy King         |       98 |
| CM015       | Turtle Island      |      106 |

* **PedNew**

|     id | status | name          | sire_id | dam_id |  YOB | sex | colour | COB |
| -----: | -----: | :------------ | ------: | -----: | ---: | --: | :----- | :-- |
|      7 |      0 | battle joined |       3 |      6 | 1959 |   1 | b.     | USA |
| 177841 |      0 | armed venus   |       3 | 452020 | 1958 |   2 | b.     | USA |
| 448331 |      0 | gentility     |       3 | 448330 | 1961 |   2 | b.     | USA |
|     15 |      0 | ack ack       |       7 |     14 | 1966 |   1 | b.     | USA |
| 383725 |      0 | jungle war    |       7 | 383724 | 1964 |   2 | b.     | USA |

* **GenotypeIDs**

| id                  |  batchID | equinomeID | SNPChip | Year of Birth | Country Reported | sex    |
| :------------------ | -------: | :--------- | :------ | ------------: | :--------------- | :----- |
| 20160112_JSB001_409 | 20160112 | JSB001_409 | SNP670  |          1971 | Ireland          | Female |
| 20141029_CM108      | 20141029 | CM108      | SNP70   |          1973 | Ireland          | Female |
| 20100101_CM002      | 20100101 | CM002      | SNP50   |          1974 | Ireland          | Male   |
| 20100101_CM104      | 20100101 | CM104      | SNP50   |          1975 | Ireland          | Female |
| 20141029_CM110      | 20141029 | CM110      | SNP70   |          1976 | Ireland          | Female |

**Pedigree #2**

This is a `.csv` file containing genotyped animals not found in PedigreeQuery.

| id              |  batchID | equinomeID | SNPChip | Year of Birth | sex    | Country Reported | Horse Name      | Sire          | Dam           | Month of Birth | Country of Birth |
| :-------------- | -------: | :--------- | :------ | ------------: | :----- | :--------------- | :-------------- | :------------ | :------------ | :------------- | :--------------- |
| 20100101_103    | 20100101 | 103        | SNP50   |          2003 | Male   | New Zealand      | Condesire       | Traditionally | Fine Feather  | September      | New Zealand      |
| 20100101_GB5    | 20100101 | GB5        | SNP50   |          1999 | Female | Ireland          | Susan Glen      | Sacrament     | SanvacGlen    | May            | Ireland          |
| 20100101_GF01   | 20100101 | GF01       | SNP50   |          2000 | Male   | Ireland          | Rossmore Castle | WitnessBox    | LatinQuarter  | June           | Ireland          |
| 20100101_GF03   | 20100101 | GF03       | SNP50   |          2003 | Male   | Ireland          | Radiator Rooney | Elnadim       | QueenofTheMay | April          | Ireland          |
| 20100101_JSB049 | 20100101 | JSB049     | SNP50   |          2004 | Female | Ireland          | Ambereen        | Lils Boy      | Aeraiocht     | February       | Ireland          |

### Output data

The preprocessing produces two output files in the `results/` folder:

`bedids2exclude.txt` **—** a file containing bed_ids of genotypes that were excluded from the analysis, formatted for PLINK QC.

| 20141029_WKS001_177 20141029_WKS001_177 |
| 20180601_BTY001_002 20180601_BTY001_002 |
| 20180601_BTY001_003 20180601_BTY001_003 |
| 20150528_SHS001_262 20150528_SHS001_262 |
| 20180111_DRK001_063 20180111_DRK001_063 |
| 20180111_RSY001_047 20180111_RSY001_047 |

`cleaned_pedigree.csv` — the cleaned and combined pedigree file.

| horse_id | status | horse_name    | sire_id | dam_id |  YOB | MOB | sex | colour | COB | equinome_id | bed_id | batch | snp_chip | country_reported | genotyped |
| -------: | -----: | :------------ | ------: | -----: | ---: | --: | --: | :----- | :-- | ----------: | -----: | ----: | -------: | ---------------: | :-------- |
|        7 |      0 | battle joined |       3 |      6 | 1959 | nan |   1 | b.     | USA |         nan |    nan |   nan |      nan |              nan | False     |
|   177841 |      0 | armed venus   |       3 | 452020 | 1958 | nan |   2 | b.     | USA |         nan |    nan |   nan |      nan |              nan | False     |
|   448331 |      0 | gentility     |       3 | 448330 | 1961 | nan |   2 | b.     | USA |         nan |    nan |   nan |      nan |              nan | False     |
|       15 |      0 | ack ack       |       7 |     14 | 1966 | nan |   1 | b.     | USA |         nan |    nan |   nan |      nan |              nan | False     |
|   383725 |      0 | jungle war    |       7 | 383724 | 1964 | nan |   2 | b.     | USA |         nan |    nan |   nan |      nan |              nan | False     |

### **Requirements**

* Python => 3.10.10
* R => 4.3.2

## Synopsis

**Preprocessing steps:**

1) **Pedigree #1 undergoes QC** , which includes:

   * Removing invalid digit entries from the colour column.
   * Checking for individuals with multiple sires or dams.
   * Ensuring parents are always older than their offspring.
   * Verifying that the sex code (1 for sires, 2 for dams) matches the individual's role (sire or dam).
2) **Genotype records from both files** are filtered:

   * EquinomeIDs are checked against the provided .fam file.
   * Duplicates based on EquinomeID are sorted by SNPChip density and date of genotyping.
   * Higher-density genotypes and the most recent data are selected.
3) **Duplicates in pedIDmatch** are resolved by:

   * Identifying all EquinomeIDs associated with the same horse_id.
   * Prioritizing duplicates based on SNPChip density and genotyping date, as done previously.
4) **All three dataframes from Pedigree #1** are merged, and columns are standardized to allow concatenation with Pedigree #2.
5) **Pedigree #2** is adjusted by aligning its columns with Pedigree #1 for concatenation.
6) **Pedigree #1 and Pedigree #2** are concatenated into a single file, `cleaned_pedigree.csv`.
7) **The `bedids2exclude.txt` file** is generated, containing all bed_ids excluded from the analysis, formatted for PLINK QC.

## Running preprocessing

To run the preprocessing pipeline, follow the steps below:

```bash
pip install -r requirements.txt

python main.py
```
