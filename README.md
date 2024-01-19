 # Code Repository

This repository contains the scripts accompanying the publication Friedrich VD & Pennitz P et al. "Neural Network-Assisted Humanization of COVID-19 Hamster scRNAseq Data Reveals Matching Severity States in Human Disease". 

## Folder Structure

- **Figure1_and_FigureS2**: Scripts for generating Figure 1 and and Figure S2.
- **Figure_2**: Scripts for generating Figure 2.
- **Figure_3**: Scripts for generating Figure 3.
     - **cross_species_Roborovski_hamster**: Scripts for generating Figure 3A.
     - **cross_species_Syrian_hamster**: Scripts for generating Figure 3B.
     - **single_species_Syrian_hamster**: Scripts for generating Figure 3C.
     - **single_species_Roborovski_hamster**: Scripts for generating Figure 3D.
- **Figure_4**: Scripts for generating Figure 4.
- **Figure5_and_FigureS5A,B**: Scripts for generating Figure 5 and and Figure S5A,B.
- **Figure6_and_FigureS5C**: Scripts for generating Figure 6 and and Figure S5C.
- **Figure7_FigureS6_and_FigureS7**: Scripts for generating Figure 7, Figure S6 and Figure S7.
- **Supplemental_Figure1**: Scripts for generating Figure S1.
- **supplement_S3_intermediate**:  Scripts for generating preliminary outputs for Figure S3.
- **Supplemental_Figure3**: Scripts for generating Figure S3.
- **supplement_S4**: Scripts for generating Figure S4.
- **helper_VAE**: Helper functions for VAE disease state matching.
- **VAE**: Scripts for model training and processing steps within VAE disease state matching.
   - **Roborovski_hamster**: Scripts for training VAE model for Roborovski hamster and human data at cell type-level.
   - **Syrian_hamster**: Scripts for training VAE model for Syrian hamster and human data at cell type-level.
   - **input_processing**: Scripts for input processing in VAE pipeline.
   - **model_selection**: Scripts for model selection.

## Usage

Navigate to the respective folder to find scripts for generating the figures or for preliminary data processing. Input and output files of the scripts can be found at zenodo XXXXXXX. An overview of the files at zenodo is provided in 'documentation_zenodo_files.txt'. Raw and processed data of Phodopus roborovskii can be found at GEO XXXXXXX.  

## Issues

If you encounter any issues or have questions, please open an issue in the repository.

## Data

Publicly available datasets that were used in this manuscript can be found at GEO:“GSE162208” (Syrian hamsters) and EGA:“EGAS00001004571” (humans) as indicated in the original publications (doi: 10.1038/s41467-021-25030-7; doi.org/10.1016/j.cell.2020.08.001).
Processed Human whole Blood ScRNAseq Data published by Jonas Schulte-Schrepping et al. 2020 (https://doi.org/10.1016/j.cell.2020.08.001) was downloaded from the FASTGenomics platform:
https://beta.fastgenomics.org/datasets/detail-dataset-1ad2967be372494a9fdba621610ad3f3#Files.
Data for Roborovski hamsters and processed data are available via GEO:XXXXXXX and Zenodo:YYYYYYYYY. 

