#!/usr/bin/env python3
"""
Generate configuration files for NGS data processing pipeline
with customizable input and output filenames and file paths
"""

import argparse
import json
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate configuration for NGS data processing')
    parser.add_argument('--experiment', type=str, required=True,
                      help='Experiment name (e.g., R3)')
    parser.add_argument('--input-dir', type=str, required=True,
                      help='Directory containing input and output zip files')
    parser.add_argument('--output-dir', type=str, required=True,
                      help='Directory for final output files')
    parser.add_argument('--input-zip-name', type=str, default='input.zip',
                      help='Name of the input zip file (default: input.zip)')
    parser.add_argument('--output-zip-name', type=str, default='output.zip',
                      help='Name of the output zip file (default: output.zip)')
    parser.add_argument('--file-path', type=str, default=None,
                      help='Custom path for the R output file (default: constructed from output-dir)')
    parser.add_argument('--pattern', type=str, default='^S...C.........C',
                      help='Regex pattern for sequence filtering')
    parser.add_argument('--output-file', type=str, default='config.json',
                      help='Output config file path')
    parser.add_argument('--motif', type=str, default='DMT',
                      help='Motif for plot generation')
    return parser.parse_args()

def normalize_path(path):
    """
    Normalize path to use forward slashes for consistency across platforms
    """
    return path.replace('\\', '/')

def main():
    args = parse_arguments()
    
    # Normalize input and output directories
    input_dir = normalize_path(args.input_dir)
    output_dir = normalize_path(args.output_dir)
    
    # Default file path or custom file path
    if args.file_path:
        file_path = normalize_path(args.file_path)
    else:
        file_path = normalize_path(os.path.join(output_dir, f"{args.experiment}DE.txt"))
    
    # Create configuration with normalized paths and customizable filenames
    config = {
        "input_zip": normalize_path(os.path.join(input_dir, args.input_zip_name)),
        "output_zip": normalize_path(os.path.join(input_dir, args.output_zip_name)),
        "input_csv": normalize_path(os.path.join(input_dir, f"{args.experiment}input.csv")),
        "output_csv": normalize_path(os.path.join(input_dir, f"{args.experiment}output.csv")),
        "output_file": normalize_path(os.path.join(input_dir, f"{args.experiment}forDE.csv")),
        "pattern": args.pattern,
        "r_script_path": "./tmm_normalization.R",
        "file_path": file_path,
        "output_path": normalize_path(os.path.join(output_dir, "reads_cpm_FC.csv")),
        "df_path": normalize_path(os.path.join(output_dir, "reads_cpm_FC.csv")),
        "plot_output": normalize_path(os.path.join(output_dir, f"{args.experiment}plot.png")),
        'motif': args.motif,
        "top_n": 100,
        "min_val": 1,
        "max_val": 1000000
    }
    
    # Write configuration to file
    with open(args.output_file, 'w') as f:
        json.dump(config, f, indent=4)
    
    print(f"Configuration created at {args.output_file}")
    print(f"\nInput zip path: {config['input_zip']}")
    print(f"Output zip path: {config['output_zip']}")
    print(f"R output file path: {config['file_path']}")
    print("\nExample usage:")
    print(f"python data_cleaning.py --config {args.output_file} --run-r")
    print(f"# To skip to step 5 (after R processing):")
    print(f"python data_cleaning.py --config {args.output_file} --skip-to-step 5")

if __name__ == "__main__":
    main()
