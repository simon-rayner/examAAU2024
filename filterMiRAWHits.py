import pandas as pd

miraw_result_file='/home/simon/Downloads/pos.csv'
df_miraw_results = pd.read_csv(miraw_result_file, delimiter='\t', comment="#")

print(str(len(df_miraw_results)))
print("s")