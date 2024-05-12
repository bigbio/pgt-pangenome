import os

import pandas as pd
from matplotlib import pyplot as plt
from pyopenms import MSExperiment, MzMLFile
from tqdm import tqdm

mzmls_path = "mzmls/"

## Read the SDRF files
sdrf_010154 = pd.read_table("PXD010154/sdrf_parquet/PXD010154.sdrf.tsv")

# Read the SDRF files
sdrf_016999_first = pd.read_table("PXD016999/PXD016999-first-instrument.sdrf.tsv")
sdrf_016999_second = pd.read_table("PXD016999/PXD016999-second-instrument.sdrf.tsv")
sdrf_016999 = pd.concat([sdrf_016999_first, sdrf_016999_second], axis=0, ignore_index=True)


def get_sdrf_info(sdrf):
    sdrf["reference_file_name"] = sdrf['comment[data file]'].str.split(".", expand=True)[0]
    file_tissue = sdrf[['characteristics[organism part]', 'reference_file_name']].drop_duplicates()
    file_tissue = file_tissue.groupby('reference_file_name')['characteristics[organism part]'].apply(list).reset_index()

    file_sample = sdrf[['source name', 'reference_file_name']].drop_duplicates()
    file_sample = file_sample.groupby('reference_file_name')['source name'].apply(list).reset_index()

    # file_tissue_map: file -> tissue
    file_tissue_map = file_tissue.set_index('reference_file_name').to_dict()['characteristics[organism part]']
    # file_sample_map: file -> sample
    file_sample_map = file_sample.set_index('reference_file_name').to_dict()['source name']

    # tissue_fileNumbers: tissue -> files_number
    tissue_fileNumbers = sdrf[['characteristics[organism part]', 'reference_file_name']].drop_duplicates()[
        'characteristics[organism part]'].value_counts()
    # tissue_sampleNumbers: tissue -> samples_number
    tissue_sampleNumbers = sdrf[['characteristics[organism part]', 'source name']].drop_duplicates()[
        'characteristics[organism part]'].value_counts()

    tissue_file = sdrf[['characteristics[organism part]', 'reference_file_name']].drop_duplicates()
    tissue_file = tissue_file.groupby('characteristics[organism part]')['reference_file_name'].apply(list).reset_index()
    # tissue_file_map: tissue -> file
    tissue_file_map = tissue_file.set_index('characteristics[organism part]').to_dict()['reference_file_name']

    return [file_tissue_map, file_sample_map, tissue_fileNumbers, tissue_sampleNumbers, tissue_file_map]

PXD010154_sdrf_info = get_sdrf_info(sdrf_010154)
PXD016999_sdrf_info = get_sdrf_info(sdrf_016999)


def plot_Spectra(RAWs_PATH, sdrf_info, path):
    tissue_MS = dict()
    for tissue, mzmls in tqdm(sdrf_info[4].items()):
        msms_num = 0
        for mzml in tqdm(mzmls, leave=False):
            file_path = RAWs_PATH + mzml.split(".")[0] + '.mzML'
            exp = MSExperiment()
            if os.path.exists(file_path):
                try:
                    MzMLFile().load(file_path, exp)
                    for scan in exp:
                        if scan.getMSLevel() == 2:
                            msms_num += 1
                except RuntimeError:
                    print(file_path + "ERROR！")
            else:
                # print(file_path + "not found")
                continue

        tissue_MS[tissue] = msms_num

    sorted_tissue_MS = sorted(tissue_MS.items(), key=lambda x: x[1], reverse=True)
    sorted_tissues = [tissue[0] for tissue in sorted_tissue_MS]
    sorted_values = [tissue[1] for tissue in sorted_tissue_MS]

    plt.figure(figsize=(10, 9))
    bar_width = 0.8
    bars = plt.bar(sorted_tissues, sorted_values, width=bar_width, align='center')
    plt.xlabel('Tissue', fontsize=14)
    plt.ylabel('Number of Spectra', fontsize=14)
    plt.xticks(rotation='vertical', fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()

    nums_list = [sdrf_info[2].get(tissue, "") for tissue in sorted_tissues]

    for i, bar in zip(nums_list, bars):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height, '%d' % int(i), ha='center', va='bottom', fontsize=8)

    plt.savefig(path, format='png')

plot_Spectra(mzmls_path, PXD010154_sdrf_info, "count-plots/PXD010154_spectra_MSMS_distribution.png")
plot_Spectra(mzmls_path, PXD016999_sdrf_info, "count-plots/PXD016999_spectra_MSMS_distribution.png")