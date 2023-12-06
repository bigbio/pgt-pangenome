# import deeplc packages
from deeplc import DeepLC
from deeplcretrainer import deeplcretrainer

# Default
from collections import Counter
import os
import urllib.request

# specific packages
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import pearsonr, spearmanr
import numpy as np

import tensorflow as tf

tf.get_logger().setLevel("INFO")

from tensorflow.python.eager import context
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")

from scipy.stats import percentileofscore


def main():
    all_gca = []
    max_inst_train = 10000

    df = pd.read_csv("grch_peptides_for_deeplc.csv")
    df.fillna("", inplace=True)
    df.index = df["seq"] + "+" + df["modifications"]

    df.sort_values("Posterior error probability", inplace=True)
    df = df[df["tr"] < 25000]

    df_gca = pd.read_csv("gca_peptides_for_deeplc.csv")
    df_gca.fillna("", inplace=True)
    df_gca.index = df_gca["seq"] + "+" + df_gca["modifications"]

    for name, sub_df in tqdm(df.groupby("sample_id")):
        sub_df_gca = df_gca[df_gca["sample_id"] == name]

        df_pred = sub_df.drop_duplicates(["seq", "modifications"], keep="first")
        df_first = sub_df.sort_values("tr").drop_duplicates(["seq", "modifications"])
        df_first.sort_values("Posterior error probability", inplace=True)

        if len(df_first.index) > max_inst_train:
            df_first = df_first.iloc[0:max_inst_train, :]

        dlc = DeepLC(
            deeplc_retrain=True,
            n_epochs=20,
            n_jobs=5,
        )

        # Perform calibration, make predictions and calculate metrics
        dlc.calibrate_preds(seq_df=df_first)
        preds_calib = dlc.make_preds(seq_df=df_first)

        df_first["preds_tr"] = preds_calib
        df_first["error"] = df_first["tr"] - df_first["preds_tr"]
        df_first["abserror"] = abs(df_first["error"])

        sub_df_gca["preds_tr"] = dlc.make_preds(seq_df=sub_df_gca)
        sub_df_gca["error"] = sub_df_gca["tr"] - sub_df_gca["preds_tr"]
        sub_df_gca["abserror"] = abs(sub_df_gca["error"])
        sub_df_gca["error_percentile"] = sub_df_gca["abserror"].apply(
            lambda x: percentileofscore(df_first["abserror"], x)
        )
        all_gca.append(sub_df_gca)

    all_gca_df = pd.concat(all_gca)

    # In[ ]:

    plt.scatter(
        all_gca_df["tr"],
        all_gca_df["preds_tr"],
        s=3,
        alpha=0.1,
        c=all_gca_df["error_percentile"],
    )
    plt.plot([500, 6500], [500, 6500], c="grey")
    plt.colorbar()
    plt.savefig("./all_predictions.png")
    plt.close()

    plt.hist(all_gca_df["error_percentile"], bins=100)
    plt.vlines(99, 0, 800, color="black", label="99th percentile")
    plt.vlines(95, 0, 800, color="grey", label="95th percentile")
    plt.legend()
    plt.savefig("./all_error_dist.png")
    plt.close()

    plt.scatter(all_gca_df["error_percentile"], all_gca_df["error"], s=1)
    plt.vlines(99, -4000, 4000, color="black", label="99th percentile")
    plt.vlines(95, -4000, 4000, color="grey", label="95th percentile")
    plt.savefig("./all_error_perc.png")
    plt.close()

    plt.scatter(
        all_gca_df["error_percentile"],
        all_gca_df["Posterior error probability"],
        s=4,
        alpha=0.1,
    )
    plt.vlines(99, 0, 0.01, color="black", label="99th percentile")
    plt.vlines(95, 0, 0.01, color="grey", label="95th percentile")
    plt.legend()
    plt.savefig("./all_error_perc_pep.png")
    plt.close()

    all_gca_df[all_gca_df["error_percentile"] < 95].to_csv(
        "gca_peptides_for_deeplc_95thperc.csv"
    )
    all_gca_df[all_gca_df["error_percentile"] < 99].to_csv(
        "gca_peptides_for_deeplc_99thperc.csv"
    )


if __name__ == "__main__":
    main()
