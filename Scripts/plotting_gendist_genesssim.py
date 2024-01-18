import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cobra
from cobra.io import read_sbml_model
import glob, os
import natsort
import ast
import json
from scipy.stats import gaussian_kde
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn.svm import SVR
import seaborn as sns
from scipy.signal import convolve
from scipy.optimize import curve_fit

gene_ess = pd.read_csv('pairwise_gene_essentiality_similarity.csv')
gene_ess = gene_ess.iloc[:,1:]

gen_dist = pd.read_csv('GeneticSimilarity245.csv')
gen_dist = gen_dist.iloc[:,1:]


gene_ess_vec = []

i = 0
for column in gene_ess.columns:
	i+=1
	for row in gene_ess.index:
		if row+i < len(gene_ess.index):
			gene_ess_vec.append(float(gene_ess.iloc[row+i][column]))
gene_ess_vec = np.asarray(gene_ess_vec)
print(gene_ess_vec)

gen_dist_vec = []
i = 0
for column in gen_dist.columns:
	i+=1
	for row in gen_dist.index:
		if row+i < len(gen_dist.index):
			gen_dist_vec.append(float(gen_dist.iloc[row+i][column]))
gen_dist_vec = np.asarray(gen_dist_vec)
print(gen_dist_vec)

gen_dist_vec = gen_dist_vec.tolist()
gene_ess_vec = gene_ess_vec.tolist()
for i in range(len(gen_dist_vec)):
	print(len(gen_dist_vec))
	print(i)
	if i < len(gen_dist_vec):
		if gen_dist_vec[i] == 1:
			del gen_dist_vec[i]
			del gene_ess_vec[i]

#gen_dist_vec = gen_dist_vec*100
#gene_ess_vec = gene_ess_vec*100
#caluclate point density
xy = np.vstack([gen_dist_vec,gene_ess_vec])
z = gaussian_kde(xy)(xy)

#log curvefit
# ignore any "invalid value in log" warnings internal to curve_fit() routine
import warnings
warnings.filterwarnings("ignore")

# alias data to match previous example
xData = np.array(gen_dist_vec, dtype=float)
yData = np.array(gene_ess_vec, dtype=float)

def func(x, m, t, b, c): # x-shifted log
    return m*np.exp(t + (c*x))+b

# these are the same as the scipy defaults
initialParameters = np.array([0, 0, 0, 0])

# curve fit the test data
fittedParameters, pcov = curve_fit(func, xData, yData, initialParameters)

modelPredictions = func(xData, *fittedParameters) 

absError = modelPredictions - yData

SE = np.square(absError) # squared errors
MSE = np.mean(SE) # mean squared errors
RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
Rsquared = 1.0 - (np.var(absError) / np.var(yData))

print('Parameters:', fittedParameters)
print('RMSE:', RMSE)
print('R-squared:', Rsquared)

print()

'''# linear regression
gen_dist_vec1 = gen_dist_vec.reshape((gen_dist_vec.shape[0],1))
data_expanded = gen_dist_vec1**(1/3)

linear_regression = LinearRegression()
linear_regression.fit(data_expanded, gene_ess_vec)
target_predicted = linear_regression.predict(data_expanded)
mse = mean_squared_error(gene_ess_vec, target_predicted)
r2_score = r2_score(gene_ess_vec, target_predicted)

'''
#p = np.polyfit(gen_dist_vec, np.log(gene_ess_vec), 1, w=np.sqrt(gene_ess_vec))
#p = np.polyfit(gen_dist_vec, np.log(gene_ess_vec), 1)
#a = np.exp(p[1])
#b = p[0]

'''from scipy.optimize import curve_fit

# Fit the function a * np.exp(b * t) + c to x and y
popt, pcov = curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, gen_dist_vec, gene_ess_vec)
a = popt[0]
b = popt[1]
c = popt[2]

x_fitted = np.linspace(np.min(gen_dist_vec), np.max(gen_dist_vec), 100)
y_fitted = a *np.exp(b *x_fitted) +c 



'''
fig, ax = plt.subplots(figsize = (15,6))

print("scattering")
plt.scatter(x = gen_dist_vec, y = gene_ess_vec, c = z, cmap = 'YlGnBu_r', s = 100)

#ax.plot(x_fitted, y_fitted, 'k')

xModel = np.linspace(min(xData), max(xData))
yModel = func(xModel, *fittedParameters)
plt.plot(xModel, yModel, c = 'k', linewidth = 2)
#plt.plot(xData, yData, c = '#b07a51', linewidth = 2)
ax.set_xlabel('Genetic Similarity \n(% Similarity)')
ax.set_ylabel('Pairwise Essential Gene Similarity \n(% Similarity)')
#ax.axhline(y = 0.18, color  = '#626363', linestyle = '--')
#ax.axhline(y = 0.34, color  = '#626363', linestyle = '--')
#ax.set_xlim(0,0.61)

#plt.yscale("log")
#ax.set_yticks(ticks = [0, 0.025, 0.050, 0.075, 0.1, 0.125, 0.150, 0.175, 0.20], labels = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0])
ax.set_xticks(ticks = [0.20, 0.30, 0.40, 0.50, 0.60, .70, .80, .90, 1], labels = [20, 30, 40, 50, 60, 70, 80, 90, 100])
#ax.set_xticks(ticks= [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6], labels = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
ax.set_yticks(ticks= [0.80, .85, 0.90, .95, 1.00], labels = [80, 85, 90, 95,100])
ax.set_ylim(0.78,1.01)

cbar = plt.colorbar()
cbar.set_label('Density of pathobiont pairs')
plt.show()

