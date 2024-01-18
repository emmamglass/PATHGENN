import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cobra
from cobra.io import read_sbml_model
import glob, os
import natsort
import ast
import json

#Determine universal essential genes across models
#create empty dictionary, keys are PATRIC IDs'

if __name__ == '__main__':
	GenomeIds = pd.read_csv('GenomeIds.csv')

	IDs = GenomeIds['Genome Ids']

	new_IDs = []

	with open('genomeids.txt') as file:
		for line in file:
			line = line.strip()
			new_IDs.append(line)

	df = pd.DataFrame()

	i = 0
	for ID in new_IDs:

		file = str(ID)+'.sbml'
		print(file)
		model = read_sbml_model(file)
		essential_genes = cobra.flux_analysis.variability.find_essential_genes(model)

		if len(list(essential_genes)) > len(df.index):
			s = pd.Series(list(essential_genes))
			df = df.reindex(s.index)
			df[str(ID)] = list(essential_genes)

		else:
			df[str(ID)] = pd.Series(list(essential_genes))

		
		i+=1
		print(df)

	df.to_csv('EssentialGenes.csv')

	
