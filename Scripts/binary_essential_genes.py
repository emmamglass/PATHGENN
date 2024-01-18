from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
import pandas as pd
import numpy as np
import math
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.neighbors import DistanceMetric

df = pd.read_csv('KO_essentialgenes.csv')

unique_essential_genes = set([])

for column in df.columns:
	for KO in df[column]:
		unique_essential_genes.add(KO)

unique_essential_genes = {x for x in unique_essential_genes if x==x}
unique_essential_genes.remove('No KEGG Ortholog Found')
unique_essential_genes = list(unique_essential_genes)

binary_essential_genes = pd.DataFrame(index = df.columns, columns = unique_essential_genes)

i = 0
for KO in binary_essential_genes.columns:
	for species in df.columns:
		for essential_gene in df[species]:
			if KO == essential_gene:
				binary_essential_genes.loc[species][KO] = 1

binary_essential_genes = binary_essential_genes.fillna(0)
		
binary_essential_genes.to_csv('Binary_gene_essentiality.csv')









