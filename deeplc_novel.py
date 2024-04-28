"""
This script is used to evaluate the performance of DeepLC on the GCA dataset. The GRCh38 peptides are used to train DeepLC,
and the GCA peptides are used to evaluate the performance of DeepLC. For each PSM, the best PEP is selected.
The training is performed by SampleID and the evaluation is performed by SampleID. The final output is a csv file with the
name gca_peptides_for_deeplc_99thperc.csv and gca_peptides_for_deeplc_95thperc.csv. The 99th percentile file contains the
peptides with an error percentile lower than 99%, and the 95th percentile file contains the peptides with an error percentile.

@author: Yasset Perez-Riverol (https://github.com/ypriverol)
"""

import click
# specific packages
import pandas as pd
import tensorflow as tf
# import deeplc packages
from deeplc import DeepLC
from matplotlib import pyplot as plt

# Default

tf.get_logger().setLevel("INFO")

from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

from scipy.stats import percentileofscore

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


@click.command("filter_deeplc", help="Run the ms2pip filtering process to remove low-quality peptides.")
@click.option("-p", "--canonical_peptide_file", type=click.Path(exists=True),
              help="The path to the canonical peptide file.", required=True)
@click.option("-n", "--novel_peptide_file", type=click.Path(exists=True), help="The path to the novel peptide file.",
              required=True)
@click.option("-o", "--output_folder", type=click.Path(), help="The path to the output folder.", required=True)
@click.option("-f", "--output_file_95perc", type=click.Path(),
              help="The name of the output file for the 95th percentile.", required=True)
@click.option("-g", "--output_file_99perc", type=click.Path(),
              help="The name of the output file for the 99th percentile.", required=True)
@click.option("-c", "--num_cores", type=int, help="The number of cores to use.", default=5)
@click.option("-v", "--verbose", is_flag=True, help="Print verbose output.")
@click.option("-s", "--num_samples",type=int, help="The number of samples to use.", default=10000)
def filter_deeplc(canonical_peptide_file: str, novel_peptide_file: str, output_folder: str, output_file_95perc: str,
                  output_file_99perc: str, num_cores: int, verbose: bool = False, num_samples: int = 10000):
    all_gca = []
    max_inst_train = 10000

    # Reading both files
    if ".gz" in canonical_peptide_file:
        df = pd.read_csv(canonical_peptide_file, sep=",",compression="gzip")
    else:
        df = pd.read_csv(canonical_peptide_file, sep=",")

    if ".gz" in novel_peptide_file:
        df_gca = pd.read_csv(novel_peptide_file, sep=",",compression="gzip")
    else:
        df_gca = pd.read_csv(novel_peptide_file, sep=",")

    # In LFQ experiments, every PSM is associated with a SampleID but if is a TMT experiment, then each PSM is associated with more
    # than one SampleID. In this case, we use reference_file_name to identify the SampleID.
    if 'sample_id' not in df.columns:
        df['sample_id'] = df['reference_file_name']
    if 'sample_id' not in df_gca.columns:
        df_gca['sample_id'] = df_gca['reference_file_name']

    ## Only use peptides that passed spectrumAI filtering
    df_gca = df_gca[(df_gca['position'] == 'non-canonical') | (df_gca['flanking_ions_support'] == 'YES')]

    df = df[df["tr"] < 25000]

    # Create indexes for both dataframes
    df.fillna("", inplace=True)
    df.index = df['sample_id'] + "+" + df["seq"] + "+" + df["modifications"]

    df_gca.fillna("", inplace=True)
    df_gca.index = df_gca['sample_id'] + "+" + df_gca["seq"] + "+" + df_gca["modifications"]

    sample_count = 0 # Counter for the number of samples processed
    # count total number of samples
    total_samples = len(df["sample_id"].unique())
    print(f"Total number of samples: {total_samples} \n")

    for name, sub_df in df.groupby("sample_id"):
        sub_df_gca = df_gca[df_gca["sample_id"] == name]

        # Use the best score grch modified peptide to train
        sub_df.sort_values("posterior_error_probability", inplace=True)
        sub_df_unique = sub_df.drop_duplicates(["seq", "modifications"])

        if len(sub_df_unique.index) > max_inst_train:
            sub_df_unique = sub_df_unique.iloc[0:max_inst_train, :]

        if num_samples < len(sub_df_unique.index):
            sub_df_unique = sub_df_unique.iloc[0:num_samples, :]

        dlc = DeepLC(
            deeplc_retrain=True,
            n_epochs=20,
            n_jobs=num_cores,verbose=verbose,
        )

        # Perform calibration, make predictions and calculate metrics
        dlc.calibrate_preds(seq_df=sub_df_unique)

        sub_df_gca["preds_tr"] = dlc.make_preds(seq_df=sub_df_gca)
        sub_df_gca["error"] = sub_df_gca["tr"] - sub_df_gca["preds_tr"]
        sub_df_gca["abserror"] = abs(sub_df_gca["error"])
        sub_df_gca["error_percentile"] = sub_df_gca["abserror"].apply(
            lambda x: percentileofscore(sub_df_gca["abserror"], x)
        )
        all_gca.append(sub_df_gca)
        sample_count += 1
        print(f"\nSample {name} done." + " % samples processed = " + str(round(sample_count / total_samples * 100, 4)) + "%" + " Number of Canonical PSMs: " + str(len(sub_df.index)) + " Number of GCA PSMs: " + str(len(sub_df_gca.index)), end="\n", flush=True) # Count samples: " + str(sample_count))

    all_gca_df = pd.concat(all_gca)

    plt.scatter(
        all_gca_df["tr"],
        all_gca_df["preds_tr"],
        s=3,
        alpha=0.1,
        c=all_gca_df["error_percentile"],
    )
    plt.plot([500, 6500], [500, 6500], c="grey")
    plt.colorbar()
    plt.savefig(output_folder + "/all_predictions.png")
    plt.close()

    plt.hist(all_gca_df["error_percentile"], bins=100)
    plt.vlines(99, 0, 800, color="black", label="99th percentile")
    plt.vlines(95, 0, 800, color="grey", label="95th percentile")
    plt.legend()
    plt.savefig(output_folder + "/all_error_dist.png")
    plt.close()

    plt.scatter(all_gca_df["error_percentile"], all_gca_df["error"], s=1)
    plt.vlines(99, -4000, 4000, color="black", label="99th percentile")
    plt.vlines(95, -4000, 4000, color="grey", label="95th percentile")
    plt.savefig(output_folder + "/all_error_perc.png")
    plt.close()

    plt.scatter(
        all_gca_df["error_percentile"],
        all_gca_df["posterior_error_probability"],
        s=4,
        alpha=0.1,
    )
    plt.vlines(99, 0, 0.01, color="black", label="99th percentile")
    plt.vlines(95, 0, 0.01, color="grey", label="95th percentile")
    plt.legend()
    plt.savefig(output_folder + "/all_error_perc_pep.png")
    plt.close()

    all_gca_df[all_gca_df["error_percentile"] < 95].to_csv(output_file_95perc, index=False,compression="gzip")
    all_gca_df[all_gca_df["error_percentile"] < 99].to_csv(output_file_99perc, index=False,compression="gzip")


cli.add_command(filter_deeplc)
if __name__ == "__main__":
    cli()
