Step1.Download your NGS file (.zip) from 48hour discovery website
Step 2. If you only have 1 library in panning, go to /scripts_for_single_library, check cleaning.ipynb as an example for deleting empty rows, stuffer, combine identical peptides etc. If you have baseline, go to /scripts_for_2_SDB_libraries. Tatget ibrary and baselines are separated through SDB code, to filter out the duplicated peptides, stuffer and blanks
Step 3. After cleaning, Take the R3forDE.csv file and run the DE_analysis.R on R for TMM
Step 4. (Optional) Take the afterDE.txt, continue with example in cleaning.ipynb to run process_afterDE for calculate the mean cpm value for further use
Step 5. (Optional) Scatter plot can be obtained by running plot_input_vs_output function in DE_plot.py 
