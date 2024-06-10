# CRAU Review

This repository contains the CRAU Review project, which includes data analysis, results, and code for evaluating cloud mask and ET quicklook methods. 

## Repository Structure

Contains various datasets used for the analysis:
- **data/**
  - `2016_all_models_HUC8_WB.csv`
  - `et_2016.csv`
  - `usgs_water_watch_runoff_mm_WY.csv`
  - `2018_final_combined_HUC8_WB.csv`
  - `2009_final_combined_HUC8_WB.csv`
  - `2013_final_combined_HUC8_WB.csv`
  - `et_2011.csv`
  - `et_2013.csv`
  - `et_2009.csv`
  - `2011_final_combined_HUC8_WB.csv`
  - `et_2018.csv`

Contains the results of the analysis:
- **results/**
  - `crau/`

 Contains source code for data processing and analysis:
- **src/**
  - `review_tools.py`
  - `review_tools_new.py`
    
Contains Jupyter notebooks for testing and analysis:
- **testing/**
  - `landsat_c2_sr_cloud_mask.ipynb`
  - `ETquicklook.ipynb`

## Getting Started

### Prerequisites

Ensure you have the following software installed:
- Python 3.x
- Jupyter Notebook
- Necessary Python packages (listed in `requirements.txt` if available)

### Installation

1. Clone the repository:
   
```bash
git clone https://github.com/Open-ET/CRAU-Review.git
cd CRAU-Review
```
   
2. Install the required Python packages:
   
```bash
pip install -r requirements.txt
```

## Usage
To explore the datasets, navigate to the data/ directory and open the relevant CSV files.

To run the Jupyter notebooks for analysis:

```bash
jupyter notebook testing/landsat_c2_sr_cloud_mask.ipynb
jupyter notebook testing/ETquicklook.ipynb
```

Source code for data processing and analysis can be found in the src/ directory. You can run the Python scripts directly or import the functions in your own scripts.

## Contributing
If you would like to contribute to this project, please follow these steps:

### Fork the repository.
Create a new branch (git checkout -b feature-branch).
Make your changes and commit them (git commit -am 'Add new feature').
Push to the branch (git push origin feature-branch).
Create a new Pull Request.

## Acknowledgements
AJ Purdy
