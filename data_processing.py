import pandas as pd

def preprocess_data(df, pattern=r'^S...C.........C', library="WL0.SDB21", baseline="WL0.SDB17"):
    """
    Preprocess the input DataFrame:
    - Filter library peptides using SDB and drop the ones that don't fit the pattern 
    - Filter baseline peptides using SDB and drop the ones that look like the library
    - Remove rows containing '*' or 'blank'
    - Trim peptides to the first 15 amino acids
    """

    df = df[
        ((df["SDB"] == library) & df["peptide"].apply(pattern)) |  # Keep matching peptides for SDB21
        ((df["SDB"] == baseline) & ~df["peptide"].apply(pattern))   # Remove matching peptides for SDB17
    ]

    df_filtered 
    df = df[df['peptide'].str.contains(pattern)]
    df = df[~df['peptide'].str.contains('\\*')]
    df = df[~df['peptide'].str.contains('blank')]
    df['peptide'] = df['peptide'].str[:15]
    return df

def aggregate_reads(df):
    """
    Aggregate reads columns by summing them up for each unique peptide.
    Dynamically handles cases where one or more of the columns are missing.
    """
    # Define the columns you want to aggregate
    columns_to_aggregate = ['reads1', 'reads2', 'reads3']

    # Determine which columns are actually present in the DataFrame
    available_columns = [col for col in columns_to_aggregate if col in df.columns]

    # If no specified columns are present, raise an error or return the original DataFrame
    if not available_columns:
        raise ValueError("None of the specified columns ('reads1', 'reads2', 'reads3') are present in the DataFrame.")

    # Perform the aggregation only on the available columns
    aggregation_dict = {col: 'sum' for col in available_columns}
    
    # Group by 'peptide' and aggregate
    return df.groupby('peptide').agg(aggregation_dict).reset_index()

def merge_data(df1, df2):
    """
    Merge two DataFrames on the 'peptide' column and fill missing values with 0.
    """
    df_merged = pd.merge(df1, df2, on='peptide', how='outer', suffixes=('_input', '_output'))
    df_merged.fillna(0, inplace=True)
    return df_merged

# Assign groups based on patterns
def assign_group(sequence):
    if pd.notnull(sequence):  # Ensure the sequence is not NaN
        if pd.Series(sequence).str.contains(library, regex=True).any():
            return 1  # Group 1: library
        else:
            return 2  # Group 2: baseline
    return None  # Handle NaN values


def process_afterDE(file_path, output_path, library=r'^S...C.........C', inputtiter = '10^10', outputtiter = '10^7'):
    """
    Process the input file to calculate averages, CPM, and fold changes for input and output reads.
    
    Parameters:
    - file_path (str): Path to the input tab-separated .txt file.
    - output_path (str): Path to save the processed .csv file.
    
    The function dynamically detects input and output columns by their names.
    """
    # Read the tab-separated file
    df = pd.read_csv(file_path, sep='\t')

    # Identify input and output columns dynamically
    input_columns = [col for col in df.columns if 'input' in col]
    output_columns = [col for col in df.columns if 'output' in col]

    if not input_columns or not output_columns:
        raise ValueError("Input or output columns not found in the file.")

    # Calculate average input and output reads
    df['input_average'] = df[input_columns].mean(axis=1)
    df['output_average'] = df[output_columns].mean(axis=1)

    # Calculate total sums of averages
    sum_average_reads_input = df['input_average'].sum()
    sum_average_reads_output = df['output_average'].sum()

    # Calculate CPM values
    df['inputcpm'] = (df['input_average'] / sum_average_reads_input) * 1_000_000
    df['outputcpm'] = (df['output_average'] / sum_average_reads_output) * 1_000_000

    # Adjust CPM values for plotting
    df['inputcpmplot'] = df['inputcpm'] + 1
    df['outputcpmplot'] = df['outputcpm'] + 1

    # Calculate fold change (FC)
    df['FC'] = df['outputcpmplot'] / df['inputcpmplot']
    
    # Group the peptides for plotting
    df['group'] = df['sequence'].apply(assign_group)
    
    df['inputtiter'] = (df['input_average'] / sum_average_reads_input) * inputtiter
    df['outputtiter'] = (df['output_average'] / sum_average_reads_oput) * outputtiter

    # Compute FractionalRecovery
    df = df[df['inputtiter'] != 0]
    df = df.dropna(subset=['inputtiter'])
    df['FractionalRecovery'] = np.where(
        df['outputtiter'] == 0,
        0,
        df['outputtiter'] / df['inputtiter']
    )

    # Drop intermediate columns
    columns_to_drop = ['input_average', 'output_average']
    df = df.drop(columns=columns_to_drop)
    
    # Save the processed DataFrame to a CSV file
    df.to_csv(output_path, index=False)

    # Display the first few rows of the processed DataFrame
    print(df.head())
