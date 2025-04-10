# NGS Data Processing Package - User Guide

This guide will walk you through how to use the NGS data processing package step by step, even if you have no experience with Python.

## Prerequisites

- Windows computer
- R installed (version 4.3.2 or higher)
- Python installed (version 3.6 or higher)

## Setup

1. **Download the package:** Download and extract the package to a folder on your computer (e.g., `C:\NGS_Analysis`)

2. **Prepare your input files:** You need two zip files:
   - `input.zip` - Contains your input/control data
   - `output.zip` - Contains your output/test data

3. **Place your files:** Create a folder for your experiment (e.g., `C:\NGS_Analysis\Experiment1`) and put both zip files there

## Running the Analysis

### Step 1: Open Command Prompt

1. Press `Win + R` on your keyboard
2. Type `cmd` and press Enter
3. Navigate to the package folder:
   ```
   cd C:\NGS_Analysis
   ```

### Step 2: Generate Configuration File

Run the following command to create a configuration file for your experiment:

```
python generate_config.py --experiment YOUR_EXPERIMENT_NAME --input-dir C:\path\to\your\zip\files --output-dir C:\path\to\desired\output\folder
```

Replace:
- `YOUR_EXPERIMENT_NAME` with your experiment name (e.g., R3)
- `C:\path\to\your\zip\files` with the folder containing your input.zip and output.zip files
- `C:\path\to\desired\output\folder` with where you want the results to be saved

Example:
```
python generate_config.py --experiment EXP1 --input-dir C:\NGS_Analysis\Experiment1 --output-dir C:\NGS_Analysis\Results\EXP1
```

### Step 3: Run the Data Processing Script

Now run the main processing script:

```
python data_cleaning.py --config config.json --run-r
```

If you want to use a specific experiment name:
```
python data_cleaning.py --config config.json --run-r --experiment YOUR_EXPERIMENT_NAME
```

### Step 4: Check the Results

After processing completes, check the output folder for:
- Processed CSV files
- Differential expression results
- Visualization plot (PNG file)

## Troubleshooting

### Error: File not found

Make sure your zip files are in the correct location and named correctly (`input.zip` and `output.zip`).

### Error: Library not found

If you get a `KeyError: 'library'` error, it means your configuration file is missing the library parameter. Regenerate your config file using the generate_config.py script.

### R Script Issues

If the R script fails to run:
1. Make sure R is installed correctly
2. The path to Rscript.exe might need to be updated in data_cleaning.py
3. Try running the R script manually as suggested in the error message

## Advanced Options

### Skipping Steps

If you want to run only specific steps:
```
python data_cleaning.py --config config.json --skip-to-step NUMBER
```
Where NUMBER is:
1. Process input zip
2. Process output zip
3. Clean data
4. TMM normalization
5. Post-differential expression processing
6. Generate plots

### Changing Parameters

To modify settings like pattern matching or motifs, edit the config.json file after it's generated or update the generate_config.py script.

## Example Workflow

1. Create experiment folder: `C:\NGS_Analysis\R3_Test`
2. Copy input.zip and output.zip to that folder
3. Generate config:
   ```
   python generate_config.py --experiment R3 --input-dir C:\NGS_Analysis\R3_Test --output-dir C:\NGS_Analysis\Results\R3
   ```
4. Run the analysis:
   ```
   python data_cleaning.py --config config.json --run-r
   ```
5. View the results in the output directory

## Need Help?

If you encounter any issues not covered in this guide, please contact technical support or refer to the documentation provided with the package.
