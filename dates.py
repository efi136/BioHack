import pandas as pd

df = pd.read_csv("sequence_metadata.csv")
date_dict = {}
for row in df.index:
    date_dict[df.at[row, 'Accession']] = df.at[row, 'Collection_Date']

print(date_dict)