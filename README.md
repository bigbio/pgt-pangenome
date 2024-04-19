# Pangenome Proteogenomics

Protegenomics analysis based on Pangenome references

The aim of this project is to search normal tissue proteomics datasets to identify novel proteins using the latest genome assemblies published via the [PanGenome project](https://www.nature.com/articles/s41586-023-05896-x).

## Project Aims
- Develop a workflow based on [quantms](https://github.com/bigbio/quantms) to reanalyze public proteomics datasets with custom proteogenomics databases.  
- Develop a workflow that enables systematic validation of novel (non-canonical peptides) using multiple existing tools.

- Performing a comprehensive analysis of multiple `normal` tissue datasets from public domain using databases generated from the latest Pangenome assemblies.
  - Compare scores and FDR for known canonical and novel canonical peptides, check distributions, etc.
  - Revisit FDR calculations and significant measures for non-canonical peptides.
  - Analyze the novel canonical, locations, gene types, other evidence for expression, etc.

- Provide a fasta database with all the novel proteins observed. 
- Draft manuscript layout and sections.

### Proteogenomics workflow

![alt text](pangenome-workflow.svg)

Workflow components: 
- **Database generation**: The proteogenomics database is created with [pypgatk](https://github.com/bigbio/py-pgatk). 
- **quanmts peptide identification**: The proteomics data is searched against the database using [quantms](https://github.com/bigbio/quantms). The workflow uses three search engines including COMET, SAGE and MSGF+ to perform the peptide identification. Percolator is then used to boost the number of peptide identifications and proteomicsLFQ or proteinQuantifier tools are used to perform the quantification and statistical filter of peptides based on the TDA (Target-Decoy approach). 
- **Post-processing**: The identified peptides are then post-processed to identify novel peptides and perform a comprehensive analysis of the results.
  - **Peptide Alignment**: The identified peptides are aligned to a Global canonical protein sequences which includes (ENSEMBL, Uniprot TrEMBL) to identify novel peptides.
  - **Spectrum Validation**: Spectrum identification validation is based on [**MS2PIP**](https://github.com/compomics/ms2pip) and Signal-to-Noise ratio (**SNR**).
  - **Variant annotation**: The identified peptides that contain Single Aminoacid variants (SAAVs) are validated using [spectrumAI tool](https://github.com/bigbio/py-pgatk)
  - **Retention time prediction**: The retention time of the identified peptides is predicted using [DeepLC](https://github.com/compomics/DeepLC)
- **Manual inspection of results using USIs** and [PRIDE USI Viewer](https://www.ebi.ac.uk/pride/archive/usi/)

### Dataset information
- Information about the protemics dataset [PXD010154](https://www.ebi.ac.uk/pride/archive/projects/PXD010154) by [Wang et al.](https://www.embopress.org/doi/full/10.15252/msb.20188503) that is used as a reference to search the databases can be obtained from this associated [SDRF file](https://github.com/bigbio/pgt-pangenome/blob/main/PXD010154.sdrf.tsv).

### Database information
- [`py-pgatk`](https://github.com/bigbio/py-pgatk) package is used to generate protein sequences for [samples](https://ftp.ensembl.org//pub/rapid-release/species/Homo_sapiens/) provided by [Liao et al.](https://www.nature.com/articles/s41586-023-05896-x) in the [Human Pangenome Reference Consortium](https://humanpangenome.org/). The steps and detailed infromation about the database generation can be accessed through the [db_generation notebook](https://github.com/bigbio/pgt-pangenome/blob/main/db_generation.ipynb).
- Additionally, [`cdna`](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz) and [`ncrna`](https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz) files from ENSEMBL release-110 were incuded in the database.

### Search tools and parameters
- [quantms](https://github.com/nf-core/quantms/tree/dev) is used to search the data and identify peptides.
- any further analyses are conducted using scripts and notebooks under the [pgt-manuscript](https://github.com/bigbio/pgt-pangenome/tree/main) repository.

### Canonical search:
- Identified peptides file: [PXD010154-pangenome-canonical-0.01.parquet](http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/proteogenomics/noncanonical-tissues-2023/PXD010154-pangenome-canonical-0.01.parquet)
- Output files from `quantms` can be accessed [here](http://ftp.pride.ebi.ac.uk/pub/databases/pride/resources/proteomes/proteogenomics/noncanonical-tissues-2023/PXD010154-pangenome-canonical-0.01/) (need update)

### Manuscript
- [Document](https://docs.google.com/document/d/1UmF2FLh54rYBSoGRUGBg4dR7QFPcCBNSjdc2mSOEcgM/edit?usp=sharing) outlining main findings and a draft version of the manuscript
- [Data files](https://docs.google.com/document/d/1fnSzAbQ0kPVrwfwWE4MHfFsM_6GG-RLsWmeno1Tb1Jg/edit) folder for supplementary data files 
- [Figure files](https://emblebi-my.sharepoint.com/:p:/g/personal/yperez_ebi_ac_uk/EQa47o2xGFFDtq0LYGQ7oHoBtcl-uRVc4H_EaoFmHAYnYA?rtime=6Lyql8EV3Eg) folder for supplementary figure files
- [Excel file](https://docs.google.com/spreadsheets/d/1KbDpwPlrJugCX2NG5XBsamLPoxiIf-pSa7vW7xkj7D8/edit?usp=sharing) for any supplementary tables


### Detailed explanation of each file
#### plots
- [plots](https://github.com/bigbio/pgt-pangenome/tree/main/plots) Figures from the paper can be accessed here
#### tables
- [tables](https://github.com/bigbio/pgt-pangenome/blob/main/tables) The tables in the paper can be accessed here
#### notebooks and scripts
- [blast_with_canonical.ipynb](https://github.com/bigbio/pgt-pangenome/blob/main/blast_with_canonical.ipynb) Obtaining variant information with canonical protein by blast
- [Count.ipynb](https://github.com/bigbio/pgt-pangenome/blob/main/Count.ipynb) Statistics on novel peptide information
- [db_generation.ipynb](https://github.com/bigbio/pgt-pangenome/blob/main/db_generation.ipynb) Detailed operations to generate the project database
- [deeplc_gca.py](https://github.com/bigbio/pgt-pangenome/blob/main/deeplc_gca.py) For each sample, the Grch38 peptide of that sample was used to calibrate the model to predict the GCA peptide
- [gca_canonical_validation.ipynb](https://github.com/bigbio/pgt-pangenome/blob/main/gca_canonical_validation.ipynb) The GCA peptide in the results was compared with canonical protein to prevent misjudgment and modified to the input format supported by deeplc
- [gca_peptides_for_deeplc.ipynb](https://github.com/bigbio/pgt-pangenome/blob/main/gca_peptides_for_deeplc.ipynb) Use deeplc to predict GCA peptides and see the effect of deeplc
- [get_observations.ipynb](https://github.com/bigbio/pgt-pangenome/blob/main/get_observations.ipynb) Obtain observation numbers for peptides
- [get_tissue_distribution.ipynb](https://github.com/bigbio/pgt-pangenome/blob/main/get_tissue_distribution.ipynb) Obtain the distribution of initial and final result peptides in each tissue
- [ms2pip_gca.py](https://github.com/bigbio/pgt-pangenome/blob/main/ms2pip_gca.py) The peptides filtered by deeplc were predicted using ms2pip
- [ms2_spectra_filter.py](https://github.com/bigbio/pgt-pangenome/blob/main/ms2_spectra_filter.py) Filter with ms2pip prediction results
#### SDRF
- [PXD010154.sdrf.tsv](https://github.com/bigbio/pgt-pangenome/blob/main/PXD010154.sdrf.tsv) SDRF corresponding to PXD010154



