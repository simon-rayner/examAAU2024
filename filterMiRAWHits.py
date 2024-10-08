import os.path

import pandas as pd
import matplotlib.pyplot as plt
import os

# set the local path to wherever you downloaded the data.
miraw_result_file='/home/simon/PycharmProjects/examAAU2024/miRAW_target_data.tsv'
miraw_result_folder = os.path.dirname(miraw_result_file)
# load the data
df_miraw_results = pd.read_csv(miraw_result_file, delimiter='\t', comment="#")

## Task 1
# How many predictions in the raw data file?
print("There are a total of <" + str(len(df_miraw_results)) + "> predictions in the raw data")

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
print("There are <" + str(df_miraw_results['GeneName'].nunique()) + "> unique 3'UTRs in the raw data")
# ans 40

# How many unique genes in the list?
df_miraw_results['miRNA'].nunique()
print("and <" + str(df_miraw_results['miRNA'].nunique()) + "> unique miRNAs ")
# ans 1627

# Task 2
#Filter the data using the target_prediction_parsing.py to only keep predictions with an
#energy > 10 and probability > 0.999.
#Question 2: How many predictions are left after filtering?
df_miraw_results_filtered = df_miraw_results[(df_miraw_results['MFE']< -10) & (df_miraw_results['Prediction']> 0.999)]
print("After filtering for MFE < -10 and PredictionProbability > 0.999,\nthere are a total of <" + str(len(df_miraw_results_filtered)) + "> predictions in the raw data")
#9282

# Task 3 Targeting Types
# There are two types of targeting events that occur in miRNA targeting. These are called
# canonical and non-canonical targeting events. Canonical binding occurs when the seed
# region of the miRNA is involved when the miRNA binds to the target. Column N reports
# whether a target event is canonical or non-canonical.

# Question 3.1: What percentage of the reported targets are non-canonical?
print("the percentage of canonical targets in the filtered data is "
      + f"{len(df_miraw_results_filtered[df_miraw_results_filtered['Canonical']=='PITA: canonnical'])/len(df_miraw_results_filtered)*100:.2f}" + "%" )
#54.3%
# 5041
#Question 3.2: Is this consistent with what is reported in the literature?


# Task 4: 3’UTR Lengths
# Question 4.1: What are the lengths of all the 3’UTRs? (show as a histogram)
# First of all, we need to get the position information for each 3'UTR. This is stored in the GeneName column
# For example
#   ENSG00000165699___TSC1_HUMAN___9-132891349-132896234n

# NOTE:
# Because of the way I am working with a subset of the original dataframe i will get warnings about
# 'A value is trying to be set on a copy of a slice from a DataFrame'
# So, this type of approach is not always a good idea. But i go ahead and using this approach because
# I am not making any changes to existing values in the dataframe. Nevertheless, you need to be aware
# about the risks with using this approach
#
df_miraw_results_filtered['start'] = df_miraw_results_filtered['GeneName'].str.split('-').str[1]
# for the stop position there is an 'n' at the end that needs to be replaced
df_miraw_results_filtered['stop'] = df_miraw_results_filtered['GeneName'].str.replace("n","").str.split('-').str[2]
df_miraw_results_filtered['length'] = df_miraw_results_filtered['stop'].astype(int) - df_miraw_results_filtered['start'].astype(int) + 1

# Now, we only want to count each gene once, so
df_miraw_results_filtered.reset_index(drop=True, inplace=True)
out = df_miraw_results_filtered.reset_index().groupby(['GeneName'])['index'].min().to_list()
df_miraw_results_unique=df_miraw_results_filtered.iloc[out]
#df_miraw_results_unique['start'] = df_miraw_results_unique['GeneName'].str.split('-').str[1]
# for the stop position there is an 'n' at the end that needs to be replaced
#df_miraw_results_unique['stop'] = df_miraw_results_unique['GeneName'].str.replace("n","").str.split('-').str[2]
#df_miraw_results_unique['length'] = df_miraw_results_unique['stop'].astype(int) - df_miraw_results_unique['start'].astype(int) + 1

# reset the index so we can subset

histLengths = df_miraw_results_unique['length'].hist(bins=8,
                                                        edgecolor='navy',  # Color of the border
                                                        color='blueviolet',
                                                        rwidth=0.8)
plt.title("Distribution of 3\'UTR lengths")
plt.xlabel('Length')
plt.ylabel('Count')
plotfileTask4File = os.path.basename(miraw_result_file).split(".")[0]+ "_task4__histogram_of_3pUTR_lengths.png"
plotfileTask4FilePath =os.path.join(miraw_result_folder, plotfileTask4File)
plt.savefig(plotfileTask4FilePath)
plt.close()
#

# Task 5: miRNA binding along the 3’UTR
# Question 5.1: Plot the distribution of locations of the start positions for all 3’UTRs.
# This is similar to Task 4, but this time we want to include all predictions
histStart = df_miraw_results_filtered['SiteStart'].hist(bins=20,
                                                        edgecolor='darkblue',  # Color of the border
                                                        color='dodgerblue',
                                                        rwidth=0.8)
plt.title('Distribution of miRNA binding starts in 3\'UTRs')
plt.xlabel('StartSite in 3\'UTR')
plt.ylabel('Frequency')

plotfileTask5File = os.path.basename(miraw_result_file).split(".")[0]+ "_task5__histogram_of_miRNA_start_sites.png"
plotfileTask5FilePath =os.path.join(miraw_result_folder, plotfileTask5File)
plt.savefig(plotfileTask5FilePath)
plt.close()

# Question 5.2: Are there any biases in the location of where the miRNAs bind on the 3’UTR?
# (for example, are there more towards the 5’ or 3’ end of the sequence?)
# Question 5.3: Are there any differences between where the canonical and non-canonical
# targeting events bind?
lCanonical = list(df_miraw_results_filtered[df_miraw_results_filtered['Canonical']=='PITA: canonnical']['SiteStart'])
lNonCanonical = list(df_miraw_results_filtered[df_miraw_results_filtered['Canonical']!='PITA: canonnical']['SiteStart'])
# Question 5.4: Display these two types of binding on the plot?
colors=['mediumslateblue', 'cornflowerblue']
names=['Canonical', 'Non-Canonical']
plt.hist([lCanonical, lNonCanonical], color=colors, label=names, density=True, stacked = True, rwidth=0.8)
plt.title('Canonical versus Non-Canonical targets by start site')
plt.xlabel('StartSite in 3\'UTR')
plt.ylabel('Frequency')
plt.legend(loc='upper right')
plotfileTask5bFile = os.path.basename(miraw_result_file).split(".")[0]+ "_task5b__histogram_of_canonical_vs_noncanonical_miRNA_start_sites.png"
plotfileTask5bFilePath =os.path.join(miraw_result_folder, plotfileTask5bFile)
plt.savefig(plotfileTask5bFilePath)
plt.close()
#plt.hist(df_miraw_results_filtered[df_miraw_results_filtered['Canonical']=='PITA: canonnical']['start'],
#         alpha=0.5, # the transaparency parameter
#         label='petal_length')
#df_miraw_results_filtered[df_miraw_results_filtered['Canonical']=='PITA: canonnical']['start'].hist(bins=20)
#df_miraw_results_filtered[df_miraw_results_filtered['Canonical']!='PITA: canonnical']['start'].hist(bins=20)
# Task 6: multiple miRNA binding events
# Question 6.1: Generate another histogram to show how many times one miRNA binds to the
# same 3’UTR for all miRNA and all 3’UTRs.
df_miraw_results_filtered['miRNA_gene']=df_miraw_results_filtered['GeneName'].str.replace("n","%")  + df_miraw_results_filtered['miRNA']
df_miraw_results_filtered['miRNA_gene'].value_counts().hist(bins=20,
                        edgecolor='darkblue', # Color of the border
                        color='violet',
                        rwidth=0.8)
plt.title('number of target/miRNAs within the same 3\'UTR')
plt.xlabel('no of targets')
plt.ylabel('Frequency')
plotfileTask6File = os.path.basename(miraw_result_file).split(".")[0]+ "_task6__histogram_of_miRNA_binding_sites_on_same_3pUTR.png"
plotfileTask6FilePath =os.path.join(miraw_result_folder, plotfileTask6File)
plt.savefig(plotfileTask6FilePath)
plt.close()
df_miraw_results_filtered['miRNA_gene'].value_counts()
df_miraw_results_filtered['miRNA_gene'].unique()

# Task 7: generate a miRNA--3’UTR targeting network
# Question 7.1: Similar to the plot above, generate a Connectivity Plot.
# this is the number of unique miRNA-gene pairs, counted by miRNA
dfU= pd.DataFrame(df_miraw_results_filtered['miRNA_gene'].unique())
dfU.columns=["genemir"]
# the first split will give me
#   hsa-miR-4422_MIMAT0018935_homo_sapiens_miR-4422
# so, if i split a second time on "_" and take the first element, this will give me
#   hsa-miR-4422
dfU["mir"] =dfU["genemir"].str.split("%").str[1].str.split("_").str[0]
dfU["mir"].value_counts().hist(bins=10,
                               edgecolor='indigo',  # Color of the border
                               color='lightskyblue',
                               rwidth=0.8)
plt.title('no of unique 3\'UTRs targeted by a single miRNA')
plt.xlabel('no of miRNAs')
plt.ylabel('no of 3\'UTR targets')
plotfileTask7File = os.path.basename(miraw_result_file).split(".")[0]+ "_task7__histogram_of_number_of_3pUTRs_per_3pUTR.png"
plotfileTask7FilePath =os.path.join(miraw_result_folder, plotfileTask7File)
plt.savefig(plotfileTask7FilePath)
plt.close()

# now i can just count up each of the occurences of each miRNA
dfU["mir"].value_counts()
# mir
# hsa-miR-874-5p     17
# hsa-miR-640        16
# hsa-miR-4531       16
# hsa-miR-6895-5p    16
# hsa-miR-1289       16

# Question 7.2: Which miRNAs target the most 3’UTRs?
dfU["mir"].value_counts()
# mir
# hsa-miR-874-5p     17
# hsa-miR-640        16
# hsa-miR-4531       16
# hsa-miR-6895-5p    16
# hsa-miR-1289       16


# Question 7.3: What is the highest number of 3’UTRs targeted by a single miRNA?
# 17 , hsa-miR-874-5p
# Question 7.4: What are the roles in disease of the three miRNAs that target the most 3’UTRs
# targets?

