# Correcting the Met Data at Duke
# ===================================
# by Bharat Sharma and Anthony Walker
# contact: sharmabd@ornl.gov
# -----------------------------------
# Running this file
# python FixingDupicateDuke.py -file_ci /Users/ud4/repos/GitHub/FATESFACE/Jupyter_Notebooks/DuplicateDukeDataCorrectIndexOnly.txt -path_data /Users/ud4/Documents/FACEMDS/MET_Data_Processing/Oren_2022_DUKE_Met/data/ -replace_file yes
# req: python3, packages: pandas, glob, argparse, numpy


import pandas as pd
import glob
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    "--file_correct_idx",
    "-file_ci",
    help="provide the path of correct index txt e.g. (./DuplicateDukeDataCorrectIndexOnly.txt)",
    type=str,
    default="./DuplicateDukeDataCorrectIndexOnly.txt",
)
parser.add_argument(
    "--path_met_data",
    "-path_data",
    help="provide the path to the parent dir of met data (~/data/)",
    type=str,
    default="~/data/",
)
parser.add_argument(
    "--replace_existing_file",
    "-replace_file",
    help="replace the existing files after correction ('yes'); else you can view results during code run",
    type=str,
    default="yes",
)

args = parser.parse_args()
filename_correct_idx = str(args.file_correct_idx)
path_data = str(args.path_met_data)
replace_existing = str(args.replace_existing_file)

df_correct_idx = pd.read_csv(f"{filename_correct_idx}", delim_whitespace=True)

selected_cols = ["Year", "JDT", "DOY", "Time"]
selected_cols_filter = ["Year", "JDT"]

all_files = glob.glob(f"{path_data}*/*_gf.csv")

for file in all_files:
    df_tmp = pd.read_csv(file, delimiter=",")
    # Check if there are common rows of duplicates
    df_common_rows = pd.merge(
        df_tmp[selected_cols],
        df_correct_idx[selected_cols_filter],
        on=["Year", "JDT"],
        how="inner",
    )
    if len(df_common_rows) > 0:
        # selecting the common dataframe in commom rows dataframe
        common_df_in_correct = df_correct_idx[
            df_correct_idx["JDT"].isin(df_common_rows["JDT"])
        ]

        # Checking if the common data frame is correct, if yes get out of the for loop
        if (
            np.sum(df_common_rows.values == common_df_in_correct[selected_cols].values)
            == 4
        ):
            continue

        # using the filter based on 'JDT' update the dataframe with correct values in 'selected_cols'
        # Find indices of common rows in df_tmp
        common_row_indices = df_tmp[df_tmp["JDT"].isin(df_common_rows["JDT"])].index
        # Updating the current file
        df_tmp.loc[common_row_indices, selected_cols] = np.asarray(
            common_df_in_correct[selected_cols].values, dtype=int
        )

        print(f"{file} has {len(df_common_rows)} common rows\n")
        print(f"Incorrect rows in the data:")
        print(df_common_rows)
        print(f"Correct rows in the data:")
        print(common_df_in_correct[selected_cols])
        if replace_existing == "yes":
            df_tmp.to_csv(f"{file}")
            print(f"{file} updated with correct values")
        else:
            print(
                f"{file} is not updated with correct values. If you want the {file} to update, add flag '-replace_file yes'."
            )
