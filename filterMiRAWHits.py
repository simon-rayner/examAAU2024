import pandas as pd

# set the local path to wherever you downloaded the data.
miraw_result_file='/home/simon/programming/python/MF915__exam/miRAW_target_data.tsv'

# load the data
df_miraw_results = pd.read_csv(miraw_result_file, delimiter='\t', comment="#")

## Task 1
# How many predictions in the raw data file?
print(str(len(df_miraw_results)))
print("s")
#ans 35231

# get column names
df_miraw_results.columns
#Index(['Unnamed: 0', 'GeneName', 'miRNA', 'SiteStart', 'SiteEnd', 'Prediction',
#       'PairStartinSite', 'SeedStart', 'SeedEnd', 'Pairs', 'WC', 'Wob', 'MFE',
#       'Canonical', 'SiteTranscript(5'to3')', 'MatureMiRNATranscript',
#       'BracketNotation', 'Filtered', 'FilteringReason',
#       'AdditionalProperties', 'Unnamed: 19'],
#      dtype='object')

# How many unique genes in the list?
df_miraw_results['GeneName'].nunique()
# ans 40

# How many unique genes in the list?
df_miraw_results['miRNA'].nunique()
# ans 1627

# Task 2
#Filter the data using the target_prediction_parsing.py to only keep predictions with an
#energy > 10 and probability > 0.999.
#Question 2: How many predictions are left after filtering?
df_miraw_results_filtered = df_miraw_results[(df_miraw_results['MFE']< -10) & (df_miraw_results['Prediction']> 0.999)]
len(df_miraw_results_filtered)
#9282

# Task 3 Targeting Types
#There are two types of targeting events that occur in miRNA targeting. These are called
#canonical and non-canonical targeting events. Canonical binding occurs when the seed
#region of the miRNA is involved when the miRNA binds to the target. Column N reports
#whether a target event is canonical or non-canonical.
#Question 3.1: What percentage of the reported targets are non-canonical?
len(df_miraw_results_filtered[df_miraw_results_filtered['Canonical']=='PITA: canonnical'])/len(df_miraw_results_filtered)*100
#54.3%
# 5041
#Question 3.2: Is this consistent with what is reported in the literature?


# Task 4: 3’UTR Lengths
# Question 4.1: What are the lengths of all the 3’UTRs? (show as a histogram)
# First of all, we need to get the position information for each 3'UTR. This is stored in the GeneName column
# For example
#   ENSG00000165699___TSC1_HUMAN___9-132891349-132896234n

ax = df_miraw_results_filtered.hist(column='session_duration_seconds', bins=25, grid=False, figsize=(12,8), color='#86bf91', zorder=2, rwidth=0.9)


# Task 5: miRNA binding along the 3’UTR
# Question 5.1: Plot the distribution of locations of the start positions for all 3’UTRs.
# Question 5.2: Are there any biases in the location of where the miRNAs bind on the 3’UTR?
# (for example, are there more towards the 5’ or 3’ end of the sequence?)
# Question 5.3: Are there any differences between where the canonical and non-canonical
# targeting events bind?
# Question 5.4: Display these two types of binding on the plot?


# Task 6: multiple miRNA binding events
# Question 6.1: Generate another histogram to show how many times one miRNA binds to the
# same 3’UTR for all miRNA and all 3’UTRs.

# Task 7: generate a miRNA--3’UTR targeting network
# Question 7.1: Similar to the plot above, generate a Connectivity Plot.
# Question 7.2: Which miRNAs target the most 3’UTRs?
# Question 7.3: What is the highest number of 3’UTRs targeted by a single miRNA?
# Question 7.4: What are the roles in disease of the three miRNAs that have the most 3’UTR
# targets?