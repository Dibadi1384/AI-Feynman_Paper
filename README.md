# AI Feynman Data Generation & Execution

## Overview
This project generates example datasets, runs **AI Feynman** on them, and collects results with runtime statistics.

## Usage

### 1. Generate Data
1. Create a folder named `example_data`:
   ```bash
   mkdir example_data
   ```
2. Run the data generation script:
   ```bash
   python Datagen.py
   ```
3. Move the generated data files into the `example_data` folder:
   ```bash
   mv <generated_files> example_data/
   ```

### 2. Run AI Feynman
Run:
```bash
./aifeynman_job.sh
```
This will:
- Execute AI Feynman for **all** data files listed in `feynex.py`
- Print runtime statistics at the end

### 3. Modify Datasets to Run
To change which datasets are processed, edit:
```bash
feynex.py
```

### 4. Results
All results are automatically saved in the `results/` folder.
