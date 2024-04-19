"""
This function is used to create the mgf file and run the ms2pipC predictions. In addition, it computes additional metrics for
each spectrum including the number of peaks, signal-to-noise, difference between the highest and lowest peaks, etc.
It contains two methods:
- create_mgf: This method creates the mgf file from the pep file and the folder with the mzMLs.
- run_ms2pip: This method runs the ms2pipC predictions for a given mgf file.
"""

from ftplib import FTP
from typing import List

import pandas as pd
from ms2pip.ms2pipC import MS2PIP
from pyopenms import *
from pyteomics import mgf
import os
import click

from tqdm import tqdm
from matplotlib import pyplot as plt
import ssl
ssl._create_default_https_context = ssl._create_unverified_context


def download_file_with_progress(ftp, remote_filename, local_filename, chunk_size=8192):
    """
    Download a file from an FTP server with a progress bar displayed in the terminal.
    :param ftp: ftp server
    :param remote_filename: remove file name
    :param local_filename: local file name
    :param chunk_size: chuck to download
    :return:
    """
    with open(local_filename, 'wb') as local_file, tqdm(
            unit='B', unit_scale=True, unit_divisor=1024, miniters=1,
            desc=f'Downloading {os.path.basename(local_filename)}', leave=True) as progress:
        def callback(data):
            local_file.write(data)
            progress.update(len(data))

        ftp.retrbinary(f'RETR {remote_filename}', callback, chunk_size)


def ftp_list_files(ftp, path='.'):
    """
    List files in a given path on an FTP server.
    :param ftp: ftp server
    :param path: path to list files
    :return: list of files
    """
    files = []
    ftp.cwd(path)
    ftp.retrlines('LIST', files.append)
    return files


def read_spectra_from_mzml(local_file_path: str, peptides: List, spectra: List) -> List:
    exp = MSExperiment()
    MzMLFile().load(local_file_path, exp)
    look = SpectrumLookup()
    look.readSpectra(exp, "((?<SCAN>)\d+$)")
    for peptide in peptides:
        index = look.findByScanNumber(peptide['scan_number'])
        spec = exp.getSpectrum(index)
        try:
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


params = {
    "ms2pip": {
        "ptm": [
            "Oxidation,15.994915,opt,M",
            "Carbamidomethyl,57.021464,opt,C",
            "Acetyl,42.010565,opt,N-term",
            "Deamidated,0.984016,opt,N",
            "Deamidated,0.984016,opt,Q",
            "TMT6plex,229.162932,opt,T",
            "TMT6plex,229.162932,opt,S",
            "TMT6plex,229.162932,opt,H",
            "TMT6plex,229.162932,opt,N-term",
            "TMT6plex,229.162932,opt,K",
        ],
        "model": "TMT",
        "frag_method": "HCD",
        "frag_error": 0.5,
        "out": "csv",
        "sptm": [],
        "gptm": [],
    }
}


def run_ms2pip_on_mgf(peptide_file: pd.DataFrame, mgf_file: str, num_cpus: int = 4):
    ## Get the ms2pipC predictions
    ms2pip = MS2PIP(pep_file=peptide_file, spec_file=mgf_file, params=params, return_results=True, num_cpu=num_cpus,
                    compute_correlations=True)
    predictions = ms2pip.run()

    predictions.replace(np.nan, 0.01, inplace=True)

    # convert pearsonr column to a new two columns depending on the type of ion column.
    ion_columns = predictions.pivot(index='spec_id', columns='ion', values='pearsonr')
    ion_columns.columns = [f'pearsonr_{ion}' for ion in ion_columns.columns]
    predictions = pd.merge(predictions, ion_columns, left_on='spec_id', right_index=True)
    # remove ion column and duplicate records
    predictions.drop(columns='ion', inplace=True)
    predictions.drop(columns='pearsonr', inplace=True)
    predictions.drop_duplicates(inplace=True)
    predictions['dot_product'] = abs(predictions['pearsonr_B']) * abs(predictions['pearsonr_Y'])

    ms2pip = MS2PIP(pep_file=peptide_file, spec_file=mgf_file, params=params, return_results=True, num_cpu=num_cpus)
    predictions_ions = ms2pip.run()
    grouped = predictions_ions.groupby(['spec_id', 'ion']).size().reset_index(name='count')

    # Pivot the table to have 'A' and 'B' as columns
    pivot_table = grouped.pivot(index='spec_id', columns='ion', values='count').reset_index()

    # Fill NaN values with 0
    pivot_table = pivot_table.fillna(0)

    # Rename the columns
    pivot_table.columns.name = None  # Remove the 'column2' header
    pivot_table.columns = ['spec_id', 'count_B', 'count_Y']

    predictions = pd.merge(predictions, pivot_table, on='spec_id')

    corrected_pearsonr_B = 0.3 * (predictions['pearsonr_B'] * np.log(predictions['count_B']))
    corrected_pearsonr_Y = 0.7 * (predictions['pearsonr_Y'] * np.log(predictions['count_Y']))
    corrected_dot_product = corrected_pearsonr_B + corrected_pearsonr_Y

    predictions['corrected_pearsonr_B'] = corrected_pearsonr_B
    predictions['corrected_pearsonr_Y'] = corrected_pearsonr_Y
    predictions['corrected_dot_product'] = corrected_dot_product
    predictions['total_ions'] = predictions['count_B'] + predictions['count_Y']

    return predictions


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    This is the main tool that gives access to all commands to convert SDRF files into pipeline-specific configuration files
    """
    pass


def compute_number_misscleavages(original_df: pd.DataFrame) -> pd.DataFrame:
    """
    This function look in the original df for the column number_misscleavages, and if it doesn't exist, it computes the
    sequence and the number of misscleavages
    :param original_df:
    :return:
    """
    if 'number_misscleavages' not in original_df.columns:
        original_df['number_misscleavages'] = original_df.apply(lambda x: x['seq'].count('K') + x['seq'].count('R'),
                                                                axis=1)
    return original_df


@click.command("create-mgf", help="Create mgf and perform the ms2pipC predictions")
@click.option("--peptide_file", help="Peptide file with observations in csv", required=True)
@click.option("--mgf_file", help="The mgf to be created", required=True)
@click.option("--mzml_path", help="the path with all the mzMLs", required=False)
@click.option("--ftp_server", help="The FTP server to download the files", required=False)
@click.option("--ftp_path", help="The FTP path to download the files", required=False)
@click.option("--local_cache_path", help="The local path to save the files", required=False)
@click.pass_context
def create_mgf(cxt, peptide_file: str, mgf_file: str, mzml_path: str = None, ftp_server: str = None,
               ftp_path: str = None, local_cache_path: str = './'):
    """
    peptide_file: Peptide file with observations in csv
    mgf_file: The mgf to be created
    mzml_path: the path with all the mzMLs
    ftp_server: The FTP server to download the files
    ftp_path: The FTP path to download the files
    local_cache_path: The local path to save the files
    """
    if mzml_path is None and (ftp_server is None or ftp_path is None):
        raise ValueError("One of the following parameters must be provided: mzml_path, ftp_server, ftp_path, "
                         "local_cache_path")

    original_df = pd.read_csv(peptide_file, sep=",")

    original_df = compute_number_misscleavages(original_df)  # add number_misscleavages column
    df_peprec = original_df[["usi", "seq", "modifications", "charge", 'scan_number', 'reference_file_name']]
    df_peprec = df_peprec.rename(columns={'usi': 'spec_id', 'seq': 'peptide'})

    # Connect to the FTP server with an anonymous account
    file_list = []

    if ftp_server is not None and ftp_path is not None:
        use_ftp = True
    else:
        use_ftp = False

    if use_ftp:
        with FTP(ftp_server) as ftp:
            ftp.login()

            # List files in the remote path
            file_list = ftp_list_files(ftp, ftp_path)
            print(f"Files in {ftp_path}:")
    else:
        # read files with mzML extension from filesystem
        file_list = [f for f in os.listdir(mzml_path) if f.endswith(".mzML")]
        print(f"Files in {mzml_path}:")

    # convert df_peprec to a dictionary group by reference_file_name
    df_peprec_dict = df_peprec.groupby('reference_file_name').apply(lambda x: x.to_dict(orient='records')).to_dict()

    # sort dictionary by the number of records in the value list, in descending order
    df_peprec_dict = dict(sorted(df_peprec_dict.items(), key=lambda x: len(x[1]), reverse=True))

    spectra = []

    count_files = 0  # count the number of files
    if use_ftp:
        with FTP(ftp_server) as ftp:
            ftp.login()
            for ref_file, peptides in df_peprec_dict.items():
                print("Reading File: {}".format(ref_file))
                for mzml_file in file_list:
                    if ref_file in mzml_file:
                        name = mzml_file.split()[8]
                        remote_file_path = os.path.join(ftp_path, name)
                        local_file_path = os.path.join(local_cache_path, name)
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
                    local_file_path = os.path.join(mzml_path, mzml_file)
                    print(mzml_file)
                    spectra = read_spectra_from_mzml(local_file_path, peptides, spectra)
                    break
            print("Number of spectra: {}".format(len(spectra)))
            count_files += 1
            if count_files == 20:
                break

    with open(mgf_file, 'w') as f:
        for spec in spectra:
            f.write("BEGIN IONS\n")
            f.write("TITLE={}\n".format(spec['TITLE']))
            f.write("PEPMASS={}\n".format(spec['PEPMASS']))
            f.write("CHARGE={}\n".format(spec['CHARGE']))
            for i in range(len(spec['MZS'])):
                f.write("{} {}\n".format(spec['MZS'][i], spec['INTENSITIES'][i]))
            f.write("END IONS\n")


def compute_signal_to_noise(intensities):
    """
    Compute the signal-to-noise ratio for a given spectrum
    :param mz_array: mz values
    :param intesity_array: intesity values
    :return:
    """
    rmsd = np.sqrt(np.mean(np.square(intensities)))

    # Calculate SNR
    snr = np.max(intensities) / rmsd

    return snr


def difference_between_highest_lowest_peaks(intensity_array):
    """
    Compute the difference in intensity between the highest and lowest peaks
    :param intensity_array:
    :return: 
    """
    sqrt_intesities = [np.sqrt(x) for x in intensity_array]
    # compute the difference in intensity between the highest and lowest peaks
    diff_highest_lowest = max(sqrt_intesities) - min(sqrt_intesities)
    return diff_highest_lowest


def get_mgf_spectrum_properties(predictions, mgf_file):
    """
    Read mgf in the title we have the usi of each spectrum. Compute all properties around the spectrum.
    :param predictions: predictions
    :param mgf_file: mgf file to be read.
    :return:
    """
    # read mgf file and get all the properties
    spectra_properties = []
    with mgf.read(mgf_file) as spectra:
        spectrum = next(spectra)
        while spectrum:
            number_peaks = len(spectrum['m/z array'])
            usi = spectrum['params']['title']
            signal_to_noise = compute_signal_to_noise(spectrum['intensity array'])
            diff_highest_lowest = difference_between_highest_lowest_peaks(spectrum['intensity array'])
            spectra_properties.append({'spec_id': usi, 'number_peaks': number_peaks,
                                       'signal_to_noise': signal_to_noise,
                                       'diff_highest_lowest': diff_highest_lowest})

            try:
                spectrum = next(spectra)
            except StopIteration:
                break

    spectra_properties_df = pd.DataFrame(spectra_properties)
    predictions = pd.merge(predictions, spectra_properties_df, on='spec_id')
    return predictions


@click.command("run-ms2pip", help="Run the ms2pip for a given peptide and mgf")
@click.option("--peptide_file", help="Peptide file with observations in csv", required=True)
@click.option("--mgf_file", help="The mgf to be created", required=True)
@click.option("--output_file", help="The output file with the predictions", required=True)
@click.option("--ms2pip_cpus", help="Number of CPUS to run ms2pip", required=False, default=4)
@click.option("--filter_aa", help="Filter peptides with less than filter_aa amino acids", required=False, default=7)
@click.pass_context
def run_ms2pip(cxt, peptide_file: str, mgf_file: str, output_file: str, ms2pip_cpus: int = 4, filter_aa: int = 7):
    original_df = pd.read_csv(peptide_file, sep=",")
    if filter_aa > 0:
        original_df = original_df[original_df['seq'].str.len() >= filter_aa]
    original_df = compute_number_misscleavages(original_df)  # add number_misscleavages column
    df_peprec = original_df[["usi", "seq", "modifications", "charge", 'scan_number', 'reference_file_name']]
    df_peprec = df_peprec.rename(columns={'usi': 'spec_id', 'seq': 'peptide'})

    predictions = run_ms2pip_on_mgf(df_peprec, mgf_file, num_cpus=ms2pip_cpus)
    predictions = get_mgf_spectrum_properties(predictions, mgf_file)

    # merge predictions with original df usi column and spec_id are the same.
    original_df = original_df.rename(columns={'usi': 'spec_id'})
    original_df = pd.merge(original_df, predictions, on='spec_id')
    original_df.rename(columns={'spec_id': 'usi'}, inplace=True)
    original_df.to_csv(output_file, sep=",", index=False)

    plt.scatter(predictions["pearsonr_B"], predictions["pearsonr_Y"], s=3, alpha=0.1)
    plt.savefig("./ms2pip_all_predictions.png")
    plt.close()

@click.command("filter-ms2pip", help="Run the ms2pip filtering process to remove low-quality peptides.")
@click.option( "-p", "--peptide_file", type=str, required=True, help="Peptide sequence to be used for the ms2pip")
@click.option("-f", "--folder_output", type=str, required=True, help="Output file after filtering ms2pip threshold")
@click.option("-o", "--output_file", type=str, required=True, help="Output file after filtering ms2pip threshold")
@click.option("--number_aa", type=int, default=8, help="Minimum number of amino acids in the peptide sequence")
def filter_ms2pip(peptide_file: str, folder_output: str, output_file: str, number_aa: int = 8):

    data = pd.read_csv(peptide_file, sep=',')
    print("Number from Previous Step and before filtering:", len(data))

    # Remove rows that the sequence_x length is less than 8 aa
    data = data[data['seq'].str.len() >= number_aa]

    # Remove ratio total_ions/number_peaks > 1 the number of total ions should be less than the number of peaks.
    data = data[data['total_ions'] / data['number_peaks'] < 1]

    # Dynamically set thresholds based on percentiles
    lower_percentile = 5
    upper_percentile = 95

    lower_threshold = data['signal_to_noise'].quantile(lower_percentile / 100)
    upper_threshold = data['signal_to_noise'].quantile(upper_percentile / 100)
    lower_threshold = round(lower_threshold, 2)
    upper_threshold = round(upper_threshold, 2)

    print(f"Lower Threshold for SNR: {lower_threshold}")
    print(f"Upper Threshold for SNR: {upper_threshold}")

    plt.hist(data['signal_to_noise'], bins=500)
    plt.title('Signal to Noise (SNR) Distribution')
    plt.vlines(lower_threshold, 0, 100, color='red', label=f'Lower 5th Percentile threshold {lower_threshold}')
    plt.vlines(upper_threshold, 0, 100, color='red', label=f'Upper 95th Percentile threshold {upper_threshold}')
    plt.xlabel('SNR')
    plt.ylabel('Frequency')
    plt.legend()
    plt.savefig(folder_output + "/signal_to_noise.png")

    # Remove rows with y-ions and b-ions wrongly correlated
    data = data[data['pearsonr_Y'] > 0.1]
    data = data[data['signal_to_noise'] > lower_threshold]
    data = data[data['signal_to_noise'] < upper_threshold]

    print("Number of rows after filtering SNR:", len(data))

    # Plot the distribution of probability scores
    plt.hist(data['corrected_pearsonr_B'], bins=50, alpha=0.6, label='b')
    plt.hist(data['corrected_pearsonr_Y'], bins=50, alpha=0.6, label='y')
    plt.xlabel('B and Y ion Pearson Correlations')
    plt.ylabel('Frequency')
    plt.title('B and Y ion Pearson Correlations')
    plt.show()

    # Choose a percentile threshold (e.g., 10th percentile)
    threshold_percentile = 5

    data_filtered_b = data
    data_filtered_y = data
    data_filtered_dp = data
    threshold_value_b = np.percentile(data_filtered_b['corrected_pearsonr_B'], threshold_percentile)
    threshold_value_y = np.percentile(data_filtered_b['corrected_pearsonr_Y'], threshold_percentile)
    threshold_value_dp = np.percentile(data_filtered_dp['corrected_dot_product'], threshold_percentile)

    # round to 2 decimal places
    threshold_value_b = round(threshold_value_b, 2)
    threshold_value_y = round(threshold_value_y, 2)
    threshold_value_dp = round(threshold_value_dp, 2)

    # Display the threshold value
    print(f"B Threshold Value (at {threshold_percentile}th percentile): {threshold_value_b}")
    print(f"Y Threshold Value (at {threshold_percentile}th percentile): {threshold_value_y}")
    print(f"Dot Product Threshold Value (at {threshold_percentile}th percentile): {threshold_value_dp}")

    # Plot the distribution of probability scores
    plt.hist(data['corrected_pearsonr_B'], bins=50, alpha=0.6, label='b')
    plt.hist(data['corrected_pearsonr_Y'], bins=50, alpha=0.6, label='y')
    plt.hist(data['corrected_dot_product'], bins=50, alpha=0.6, label='MS2PIP Adjusted Score')
    plt.xlabel('B and Y ion Pearson Distribution, and MS2PIP Adjusted Score')
    plt.vlines(threshold_value_dp, 0, 1600, color='red',
               label=f'{threshold_percentile}th Percentile  MS2PIP Adjusted Score {threshold_value_dp}')
    plt.ylabel('Frequency')
    plt.title('B and Y ion Pearson Correlations')
    plt.legend()
    plt.savefig(folder_output + "/ms2pip_filtered_score.png")

    data = data[data['corrected_dot_product'] > threshold_value_dp]
    data.to_csv(folder_output + '/' + output_file, sep=',', index=False)

cli.add_command(create_mgf)
cli.add_command(run_ms2pip)
cli.add_command(filter_ms2pip)

if __name__ == "__main__":
    cli()
