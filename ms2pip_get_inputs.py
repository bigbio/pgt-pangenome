from ftplib import FTP
from typing import List

import pandas as pd
from ms2pip.ms2pipC import MS2PIP
from pyopenms import *
import os

from tqdm import tqdm
from matplotlib import pyplot as plt
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

use_ftp = False

def download_file_with_progress(ftp, remote_filename, local_filename, chunk_size=8192):
    with open(local_filename, 'wb') as local_file, tqdm(
            unit='B', unit_scale=True, unit_divisor=1024, miniters=1,
            desc=f'Downloading {os.path.basename(local_filename)}', leave=True) as progress:
        def callback(data):
            local_file.write(data)
            progress.update(len(data))

        ftp.retrbinary(f'RETR {remote_filename}', callback, chunk_size)


def ftp_list_files(ftp, path='.'):
    files = []
    ftp.cwd(path)
    ftp.retrlines('LIST', files.append)
    return files


original_df = pd.read_csv("gca_peptides_for_deeplc_95thperc_observations.tsv", sep="\t")
df_peprec = original_df[["usi", "seq", "modifications", "charge", 'scan_number', 'reference_file_name']]

df_peprec = df_peprec.rename(columns={'usi': 'spec_id', 'seq': 'peptide'})

params = {
    "ms2pip": {
        "ptm": [
            "Oxidation,15.994915,opt,M",
            "Carbamidomethyl,57.021464,opt,C",
            "Acetyl,42.010565,opt,N-term",
            "Deamidated,0.984016,opt,N",
            "Deamidated,0.984016,opt,Q",
        ],
        "model": "HCD2019",
        "frag_method": "HCDch2",
        "frag_error": 0.04,
        "out": "csv",
        "sptm": [],
        "gptm": [],
    }
}

## Get mzml file lists from the filesystem in the folder raw files
# Remote path on the FTP server
ftp_host = 'ftp.pride.ebi.ac.uk'
remote_path = '/pub/databases/pride/resources/proteomes/proteogenomics/noncanonical-tissues-2023/mzmls/'

## Get mzml file lists from the filesystem in the folder raw files
filesystem_path = '/Users/yperez/work/pgt-pangenome/'

# Local path to save downloaded files
cache_local_path = './'

# Connect to the FTP server with an anonymous account
file_list = []

if use_ftp:
    with FTP(ftp_host) as ftp:
        ftp.login()

        # List files in the remote path
        file_list = ftp_list_files(ftp, remote_path)
        print(f"Files in {remote_path}:")
else:
    # read files with mzML extesnion from filesystem
    file_list = [f for f in os.listdir(filesystem_path) if f.endswith(".mzML")]
    print(f"Files in {filesystem_path}:")

# convert df_peprec to a dictionary group by reference_file_name
df_peprec_dict = df_peprec.groupby('reference_file_name').apply(lambda x: x.to_dict(orient='records')).to_dict()
# sort dictionary by the number of records in the value list, in descending order
df_peprec_dict = dict(sorted(df_peprec_dict.items(), key=lambda x: len(x[1]), reverse=True))


# Get spectra from the mzML in the dd_peprec_dict
spectra = []

def read_spectra_from_mzml(local_file_path: str, peptides: List, spectra: List) -> List:
    exp = MSExperiment()
    MzMLFile().load(local_file_path, exp)
    for peptide in peptides:
        for spec in exp:
            scan_id = "scan={}".format(peptide['scan_number'])
            try:
                if scan_id in str(spec.getNativeID()):
                    spectra_mgf = {
                        "TITLE": peptide['spec_id'],
                        "PEPMASS": spec.getPrecursors()[0].getMZ(),
                        "CHARGE": str(peptide['charge']) + '+',
                        "INTENSITIES": spec.get_peaks()[1].tolist(),
                        "MZS": spec.get_peaks()[0].tolist(),
                    }
                    print("Reading spectrum: {}".format(peptide['spec_id']))
                    spectra.append(spectra_mgf)
            except Exception as e:
                ms_level = spec.getMSLevel()
                print("Error reading spectrum: {} {} which is MS level {}".format(peptide['spec_id'], e, ms_level))


    print("Number of spectra: {} until the file {}".format(len(spectra), local_file_path))
    return spectra

count_files = 0 # count the number of files
if use_ftp:
    with FTP(ftp_host) as ftp:
        ftp.login()
        for ref_file, peptides in df_peprec_dict.items():
            print("Reading File: {}".format(ref_file))
            for mzml_file in file_list:
                if ref_file in mzml_file:
                   name = mzml_file.split()[8]
                   remote_file_path = os.path.join(remote_path, name)
                   local_file_path = os.path.join(cache_local_path, name)
                   # Download the file
                   if not os.path.exists(local_file_path):
                        download_file_with_progress(ftp, remote_file_path, local_file_path)
                   print(mzml_file)
                   spectra = read_spectra_from_mzml(local_file_path, peptides, spectra)
                   break
                print("Number of spectra: {}".format(len(spectra)))
            count_files += 1
            if count_files == 20:
               break
else:
    for ref_file, peptides in df_peprec_dict.items():
        print(ref_file)
        for mzml_file in file_list:
            if ref_file in mzml_file:
                local_file_path = os.path.join(filesystem_path, mzml_file)
                print(mzml_file)
                spectra = read_spectra_from_mzml(local_file_path, peptides, spectra)
                break
        print("Number of spectra: {}".format(len(spectra)))
        count_files += 1
        if count_files == 20:
            break

## Write spectra to a mgf file
mgf_file = "gca_peptides_for_deeplc_95thperc_observations.mgf"
with open(mgf_file, 'w') as f:
    for spec in spectra:
        f.write("BEGIN IONS\n")
        f.write("TITLE={}\n".format(spec['TITLE']))
        f.write("PEPMASS={}\n".format(spec['PEPMASS']))
        f.write("CHARGE={}\n".format(spec['CHARGE']))
        for i in range(len(spec['MZS'])):
            f.write("{} {}\n".format(spec['MZS'][i], spec['INTENSITIES'][i]))
        f.write("END IONS\n")

## Get the ms2pipC predictions
ms2pip = MS2PIP(pep_file=df_peprec, spec_file=mgf_file, params=params, return_results=True, num_cpu=4,
                compute_correlations=True)
predictions = ms2pip.run()

predictions.replace(np.nan, 0.0001, inplace=True)

# convert pearsonr column to a new two columns depending on the type of ion column.
ion_columns = predictions.pivot(index='spec_id', columns='ion', values='pearsonr')
ion_columns.columns = [f'pearsonr_{ion}' for ion in ion_columns.columns]
predictions = pd.merge(predictions, ion_columns, left_on='spec_id', right_index=True)
#remove ion column and duplicate records
predictions.drop(columns='ion', inplace=True)
predictions.drop(columns='pearsonr', inplace=True)
predictions.drop_duplicates(inplace=True)
predictions['dot_product'] = abs(predictions['pearsonr_B']) * abs(predictions['pearsonr_Y'])

#merge predictions with original df usi column and spec_id are the same.
original_df = original_df.rename(columns={'usi': 'spec_id'})
original_df = pd.merge(original_df, predictions, on='spec_id')
original_df.rename(columns={'spec_id': 'usi'}, inplace=True)
original_df.to_csv("gca_peptides_for_deeplc_95thperc_observations_ms2pip.tsv", sep="\t", index=False)

plt.scatter(
        predictions["pearsonr_B"],
        predictions["pearsonr_Y"],
        s=3,
        alpha=0.1,
    )
plt.savefig("./ms2pip_all_predictions.png")
plt.close()
