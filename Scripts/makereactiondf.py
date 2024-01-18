import numpy as np 
import pandas as pd
import argparse
import cobra
import glob
import csv
from sklearn.neighbors import DistanceMetric
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import pairwise_distances
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import matplotlib as mpl
import scipy.cluster.hierarchy as shc
import plotly.figure_factory as ff
from sklearn.manifold import TSNE

#first we will read in all sbmls, grab rxns present
def make_rxn_list():
	file_list = glob.glob("*.sbml")
	rxn_list = []
	strain = 0
	for file in file_list:
		model = cobra.io.read_sbml_model(file, low_memory=False)
		rxns = model.reactions
		for rxn in rxns:
			reaction = str(rxn).split()
			reaction = str(reaction[0:1])
			reaction = str(reaction).replace("['", '')
			reaction = str(reaction).replace("']",'')
			reaction = str(reaction).replace("'", '')
			reaction = str(reaction).replace(":", '')
			print(reaction)
			rxn_list.append(str(reaction))
		strain += 1
		print(strain)

	all_rxn_list = [i for n, i in enumerate(rxn_list) if i not in rxn_list[:n]]

	print(all_rxn_list)
	return all_rxn_list

def make_rxn_df():
	file_list = glob.glob("*.sbml")
	rxn_list= []
	strain = 0

	file = open('all_rxn_list.txt', 'r')
	reactome1 = file.read().split(',')
	reactome1[-1] = reactome1[-1].strip()
	print(reactome1)
	reactome = []
	for reaction in reactome1:
		reaction = str(reaction).replace("[", '')
		reaction = str(reaction).replace("'", '')
		reaction = str(reaction).replace("]", '')
		reaction = str(reaction).replace(" ", '')
		reactome.append(str(reaction))
	print(reactome)

	rows = file_list
	dfrows = []
	for row in rows:
		row = row.strip('.sbml')
		dfrows.append(row)

	df = pd.DataFrame(index = dfrows, columns = reactome)

	row_num = 0
	for file in file_list:
		model = cobra.io.read_sbml_model(file, low_memory = False)
		reactions = model.reactions
		rxn_presence = []
		for i in range(len(reactome)):
			rxn = df.columns[i]
			if rxn in reactions:
				rxn_presence.append(1)
			else:
				rxn_presence.append(0)
		#print(rxn_presence)
		df.loc[dfrows[row_num]] = rxn_presence
		row_num+=1
		print(row_num)

	df.to_csv('reactionpresence.csv')
	return df


##Main
make_rxn_df()
