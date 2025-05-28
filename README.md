# FluxesErrorHierarchy

## Overview
This project evaluates the performance of various configurations of the Canadian Regional Climate Model (CRCM6-GEM5) in simulating surface turbulent fluxes using data from AmeriFlux stations across North America.

## Prerequisites

- Python 3.8
- All dependencies are listed in `requirements.txt`

Install with:
```bash
python3.8 -m pip install -r requirements.txt
```
## Usage

### Data
The script expects the data from doi: [https://doi.org/10.5683/SP3/VWMCY0]

### Configuration
**Note:**
Edit the path variables in `config.py` and set them to the actual locations of your data and output folders before running the code.

### Preprocessing
Raw simulation and observational data are put in `preprocess.py` to organize and filter them by station, regime, and time. It computes summary statistics and saves the processed data. This step ensures that the main script can run quickly and efficiently.

> **⚠️ Warning:**  
> Running `preprocess.py` will generate approximately **3GB** of data in the output directory.  
> Make sure you have enough disk space before running this step.
```bash
python preprocess.py
```

### Main code
Execute the main script to start the evaluation:
```bash
python main.py
```

## Acknowledgements

This project is adapted from [WindErrorHierarchy](https://github.com/timwhittaker/WindErrorHierarchy).
