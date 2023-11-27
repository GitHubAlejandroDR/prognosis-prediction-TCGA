# Colorectal Cancer Prognosis Prediction - TCGA Data

## Table of Contents

- [Project Description](#project-description)
- [Technologies](#technologies)
  - [Key Features](#key-features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [How to use it](#how-to-use-it)
- [Project Structure](#project-structure)




## Project Description

Machine learning (ML) techniques are transforming many aspects of our society nowadays. In the medical field and more specifically in oncology ML is being applied from laboratories to clinical practice. In this project, 4 ML models will be applied in the prognosis of 5-year survival using colorectal cancer omics data. At the same time, the bias in the predictions of the models will be evaluated considering as study variable 'RACE' and as sensitive attribute 'Black and African American'. The dataset used corresponds to a real cohort of 594 patients and was downloaded from the public repository The Cancer Genome Atlas TCGA. The integration of 3 types of omics data (RNAseq, degree of methylation and abundance of microorganisms in microbiome) associated with each sample was applied in the development of the ML models. The most significant variables from the omics datasets were selected by 2 sequentially applied feature selection methods. The 30 most significant feature from each omics dataset were applied in the development of ML models. The ML models were implemented using the lgbm libraries of Python and Caret of R. The performance of the models was evaluated by 5x2 Cross Validation. A selection of metrics were applied both in evaluating the performance of the ML models and in detecting bias in their predictions. The results showed significant differences in performance between the families of ML models applied. The analysis of the bias in the predictions was influenced by the unbalanced proportion of variables in the 'Race' variable, causing a decrease in the robustness of the results in the applied metrics.

## Technologies

| **Step** | **Technology** |
|---------|---------------|
| Configuration Modularity | R |
| Documentation | Roxygen|
| Data Analysis | R |
| Differential Analysis | limma |
| Model training | Caret, lgbm (Python) |
| Bias Analyisis | Fairness |
| Report | RMarkdown |

### Key Features

- **Configuration Modularity:** File modules for the main project components (configuration variables, functions, libraries, execution), ensuring reproducibility and modularity.

- **Project Steps Modularity:** RMD files for each one of the ML project steps (load, data analysis, preprocessing, feature selection, data integration, training, bias analysis).

- **Data Analysis:** Data analysis in R and its broad versatile libraries and functionalities.

- **Double Feature Selection:** Two-step feature selection, using the limma package and information gain filter.

- **Omics Data Integration:** Early data integration of the genomic, methylation, and microbiome datasets.

- **Model Training Versatility:** Model training flexibility using the Caret R library and lgbm Python library.

- **Bias Analysis:** Bias model prediction analysis, using the fairness R package.

- **Report:** HTML representation of the entire project workflow, through RMD modules files.


## Getting Started

<!--
![Template](docs/media/clinical-cancer-template_page-0001.jpg)

![Template](docs/mediaprueba_animated.gif)

<img src="docs/mediaprueba_animated.gif" width="300" alt="GitHub Logo">
 
-->

### Prerequisites

Before you begin working with this project, make sure you have the following prerequisites installed on your system:

- [R](https://cran.r-project.org/): This project is built using R, and you'll need R version 4.3.1 or higher. You can download R from the official [R website](https://cran.r-project.org/).

- [Python](https://www.python.org/downloads/): LGBM model training is developed in Python, and you'll need Python 3.6 or higher. You can download Python from the official [Python website](https://www.python.org/downloads/).

- Download the dataset Colorectal Adenocarcinoma (TCGA, PanCancer Atlas) from [cBioPortal](https://www.cbioportal.org/study/clinicalData?id=coadread_tcga_pan_can_atlas_2018).

- Adapt the directory paths to your working directory in `/code/configure.R`.

Once you've installed R and Python, downloaded the data, and adapted the directories, you can proceed to work with this project.

### How to use it

To get project up, follow these steps:

1. **Clone the project repository**:

   ```shell
   git clone https://github.com/GitHubAlejandroDR/prognosis-prediction-TCGA.git

2. Change Directory:

   ```shell
   cd prognosis-prediction-TCGA

3. Execution:

Ensure you have R installed. If not, install it following the official instructions.

```shell
   Rscript execute.R configure.R
```

## Project Structure

Below is an overview of the key directories and their respective roles within the project:

- **`/RMDs`**: RMD files for each one of the ML project steps

- **`/code`**: R files for project configuration, setting, and execution.

- **`/plots`**: Plots and summaries outputs

- **`/results`**: Date-organized project execution results. Results are formed by a set of model results .CSV tables and the final .html workflow report.
  
  - **`/project_Name_results_DATE`**: Result of one execution

- **`/data`**: Project data

  - **`/raw`**: Raw data from cBioPortal
  - **`/processed`**: Data processed and the final dataset resulted from omic data integration.

