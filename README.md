# WindErrorHierarchy

## Overview
This project evaluates the performance of various configurations of the Canadian Regional Climate Model (CRCM6/GEM5) in simulating 10-meter wind speeds using data from AmeriFlux stations across North America. The study employs a hierarchy of error metrics and highlights the challenges in simulating near-surface wind speeds.

## Prerequisites

- Python 3.8
- Required Python packages: `numpy`, `matplotlib`, `tqdm`, `datetime`, `xarray`

Install the necessary packages using:
```bash
python3.8 -m pip install -r requirements.txt
```
## Usage

### Input data

The script expects the data from doi: [https://doi.org/10.5683/SP3/VWMCY0]

### Running the Script
First adjust the paths in 
```bash
config.py
```
Then preprocess the data with:
```bash
python preprocess.py
```
Then execute the main script to start the evaluation:
```bash
python main.py
```
