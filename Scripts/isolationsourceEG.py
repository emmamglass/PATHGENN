import pandas as pd
import numpy as np

df = pd.read_csv('Essential_Genes_Isolation_Source.csv')


unique_isosource = df.IsolationSource.unique()

DataFrameDict = {elem: pd.DataFrame() for elem in unique_isosource}

for key in DataFrameDict.keys():
	DataFrameDict[key] = df[:][df.IsolationSource == key]

df1 = pd.DataFrame()
for key, value in DataFrameDict.items():
	DataFrameDict[key] = DataFrameDict[key].loc[:, (DataFrameDict[key] != 0).any(axis=0)]
	DataFrameDict[key] = DataFrameDict[key].reset_index()
	for column in DataFrameDict[key]:
		df = DataFrameDict[key]
		total = df[column].sum()
		length = len(df[column])
		if column != 'IsolationSource':
			if total/length >= .5:
				print(key, column)
		'''if (df[column] == df[column][0]).all():
			#print(df)
			print(key, column)'''

