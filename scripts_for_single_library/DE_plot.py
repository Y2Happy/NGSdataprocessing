
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_input_vs_output(df, top_n=100, motif='DMT', min_val=1, max_val=10**6, output_file='R3plot.png'):
    """
    Plots a log scale scatter plot of Input vs Output ppm with highlighted top N peptides and motif matches.
    
    Parameters:
    - df: pandas DataFrame containing the data with columns 'inputppmplot', 'outputppmplot', and 'sequence'.
    - top_n: int, the number of top entries to highlight (default is 100).
    - motif: str, the peptide sequence motif to highlight (default is 'DMT').
    - min_val: int, the minimum value for the diagonal line (default is 1).
    - max_val: int, the maximum value for the diagonal line (default is 10**6).
    - output_file: str, the file name for saving the plot (default is 'R2plot.png').
    """
    
    # Filter rows that contain the motif
    motif_rows = df[df['sequence'].str.contains(motif, case=False, na=False)]
    
    # Get top N rows based on the 'outputcpm' column
    top_n_rows = df.nlargest(top_n, 'outputcpm')

    # Find overlap between top N rows and motif matches
    overlap = pd.merge(top_n_rows, motif_rows, how='inner')
    
    # Create diagonal line points
    line_points = np.linspace(min_val, max_val, 100)
    
    # Create the plot
    plt.figure(figsize=(8, 8))
    
    # Scatter plot of all data
    sns.scatterplot(x='inputcpmplot', y='outputcpmplot', data=df, label='Peptides', linewidth=0.2, s=50)
    
    # Highlight the top N peptides in red
    sns.scatterplot(x='inputcpmplot', y='outputcpmplot', data=top_n_rows, color='red', label=f'Top {top_n}')
    
    # Highlight the motif-containing rows in gold
    sns.scatterplot(x='inputcpmplot', y='outputcpmplot', data=motif_rows, color='gold', edgecolor='black', linewidth=0.1, label=f'{motif} Motif')
    
    # Plot the diagonal reference line
    plt.plot(line_points, line_points, color='red', linestyle='--')
    
    # Set log scale
    plt.xscale('log')
    plt.yscale('log')
    
    # Ensure the aspect ratio is equal
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Set plot labels and title
    plt.title('Input vs Output ppm (Log Scale)')
    plt.xlabel('Input (log scale)')
    plt.ylabel('Output (log scale)')
    
    # Show legend
    plt.legend()
    
    # Save the plot to a file
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    # Display the plot
    plt.show()
