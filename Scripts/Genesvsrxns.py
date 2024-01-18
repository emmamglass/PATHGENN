import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.optimize import curve_fit

data = pd.read_csv('taxonomyinfo.csv')
print(data)

genes = data['Genes']
genes = genes.iloc[0:913]
print(genes)
reactions = data['Reactions']
reactions = reactions.iloc[0:913]
print(reactions)
class_data = data['Class']
class_data = class_data.iloc[0:913]
#log curvefit
# ignore any "invalid value in log" warnings internal to curve_fit() routine
import warnings
warnings.filterwarnings("ignore")

# alias data to match previous example
xData = np.array(genes, dtype=float)
yData = np.array(reactions, dtype=float)

def func(x, a, b, c): # x-shifted log
    return a*np.log(x + b)+c

# these are the same as the scipy defaults
initialParameters = np.array([1.0, 1.0, 1.0])

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

colors = sns.color_palette(palette = ['#d76777', '#60b19f', '#5f9cbf', '#5fa7b9', '#5fb19f', '#977fc0', '#fc9756', '#a777bc', '#e7795f', '#bac462', '#67b26b', '#7eb56a', '#fdc778', '#ffffff', '#ce6fa1', '#9cb970', '#fdad58', '#dcc969'], n_colors = 18)
sns.set_palette(colors)

fig, axes = plt.subplots(figsize = (3,2))
#axes.plot(y = 1407, color = '#B8B8B8', linestyle = ':', zorder = 0)
xModel = np.linspace(min(xData), max(xData))
yModel = func(xModel, *fittedParameters)
sns.scatterplot(ax = axes, x='Genes', y='Reactions', data = data, hue = 'Class', s= 20)
plt.plot(xModel, yModel, c = 'k', linewidth = 1, alpha = 0.7)
plt.legend([],[], frameon=False)
#axes.set_xticklabels(['','','','','','','','',''])
axes.set_xlabel('Number of Genes')
axes.set_ylabel('Number of Reactions')
#axes.tick_params(axis = 'x', which = 'both', bottom = False)

plt.show()