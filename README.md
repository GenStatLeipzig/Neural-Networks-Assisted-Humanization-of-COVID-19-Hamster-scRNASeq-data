\# Code Repository

This repository contains the scripts accompanying the publication Friedrich VD & Pennitz P et al. "Neural Network-Assisted Humanization of COVID-19 Hamster scRNAseq Data Reveals Matching Severity States in Human Disease".

## Folder Structure

-   **Figure1_and_FigureS2**: Scripts for generating Figure 1 and and Figure S2.
-   **Figure2**: Scripts for generating Figure 2.
-   **Figure3**: Scripts for generating Figure 3.
    -   **cross_species_Roborovski_hamster**: Scripts for generating Figure 3A.
    -   **cross_species_Syrian_hamster**: Scripts for generating Figure 3B.
    -   **single_species_Syrian_hamster**: Scripts for generating Figure 3C.
    -   **single_species_Roborovski_hamster**: Scripts for generating Figure 3D.
-   **Figure4**: Scripts for generating Figure 4.
-   **Figure5_and_FigureS5A,B**: Scripts for generating Figure 5 and and Figure S5A,B.
-   **Figure6_and_FigureS5C**: Scripts for generating Figure 6 and and Figure S5C.
-   **Figure7_FigureS6_and_FigureS7**: Scripts for generating Figure 7, Figure S6 and Figure S7.
-   **FigureS1**: Scripts for generating Figure S1.
-   **FigureS3**: Scripts for generating Figure S3.
    -   **Roborovski_hamster**: Scripts for generating preliminary output for Roborovski hamster and human data for Figure S3.
    -   **Syrian_hamster**: Scripts for generating preliminary output for Syrian hamster and human data for Figure S3.
-   **FigureS4**: Scripts for generating Figure S4.
-   **VAE**: Scripts for model training and processing steps within VAE disease state matching.
    -   **Roborovski_hamster**: Scripts for training VAE model for Roborovski hamster and human data at cell type-level.
    -   **Syrian_hamster**: Scripts for training VAE model for Syrian hamster and human data at cell type-level.
    -   **input_processing**: Scripts for input processing in VAE pipeline.
    -   **model_selection**: Scripts for model selection.
    -   **helper_VAE**: Helper functions for VAE disease state matching.
    -   **tool_versions**: Tool versions for the VAE pipeline. VAE pipeline was mainly carried out with tool versions described in *tool_versions_VAE.txt*. Minor steps were carried out with tool versions shown in *tool_versions_local_environment1.txt* and *tool_versions_local_environment2.txt*. In case problems arise when setting up the respective environemnts, feel free to contact the contributors or open an issue in the repository.
-   **hamster_scRNAseq**: Scripts for initial provision of Syrian hamster and Roborovski hamster Seurat input files.
-   **Supplementary_Tables:** Scripts documenting content creation of Supplementary tables

## Usage

Navigate to the respective folder to find scripts for generating the figures . Input data for these scripts is provided via zenodo (see section *Data*). Scripts used for creation of input data are also provided in respective folders, see file *documentation_zenodo_files.txt* for details.

## Data

Input and output files of the scripts can be found at **Zenodo: 10.5281/zenodo.10580154**. An overview of all files at Zenodo is available in *documentation_zenodo_files.txt*. Raw and processed data of Phodopus roborovskii (Roborovski hamster) can be found at **GEO GSE253845**. Publicly available datasets that were used in this manuscript can be found at **GEO: GSE162208** (Syrian hamsters) and **EGA: EGAS00001004571** (humans) as indicated in the original publications (doi: 10.1038/s41467-021-25030-7; doi.org/10.1016/j.cell.2020.08.001). Processed human whole blood scRNAseq data published by Schulte-Schrepping et al. 2020 (<https://doi.org/10.1016/j.cell.2020.08.001>) was downloaded from the FASTGenomics platform: <https://beta.fastgenomics.org/datasets/detail-dataset-1ad2967be372494a9fdba621610ad3f3#Files>.

## Issues

If you encounter any issues or have questions, please open an issue in the repository.

## To Dos

- version numbers packages
- add installation guide, readme  



## Disclaimer

THE SOFTWARE AND THE DATA IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
