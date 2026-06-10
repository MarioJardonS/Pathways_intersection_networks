# importing required libraries
# Tratamiento de datos
# ==============================================================================
import numpy as np
import pandas as pd
import statsmodels.api as sm

# Gráficos
# ==============================================================================
import matplotlib.pyplot as plt
import matplotlib.font_manager
from matplotlib import style
style.use('ggplot') or plt.style.use('ggplot')

# Preprocesado y modelado
# ==============================================================================
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import scale

# Configuración warnings
# ==============================================================================
import warnings
warnings.filterwarnings('ignore')


dataset = pd.read_table("../../distancia.tsv" )

dataset = dataset.iloc[ : , 1: ] 

scaler = StandardScaler()

X = scaler.fit_transform(dataset)

PCA = PCA(n_components=2)
components = PCA.fit_transform(X)
PCA.components_

metadata = pd.read_csv("../Data/metadata.csv" )
file = open("../Data/argumentos.txt")

Y = []

def cl(diag):
	if diag == "Healthy":
		return 1 
	if diag == "UC":
		return 2
	if diag == "CD":
		return 3
	if diag == "Obese":
		return 4

for line in file:
	muestra = line[0:(len(line)-1)]
	for i in range(metadata.shape[0]):
		if metadata.loc[ i , "SampleID" ] == muestra:
			c = cl(metadata.loc[i , "Diagnosis"]) 
			Y.append(c)


Z = pd.Series(Y)


componentsDf = pd.DataFrame(data = components, columns = ['PC1', 'PC2'])
pcaDf = pd.concat([componentsDf, Z], axis=1)

plt.figure(figsize=(12, 6))
plt.scatter(data=pcaDf, x="PC1", y="PC2" , c = Y)
plt.savefig("pca_edges.png")
#sns.scatterplot(data=pcaDf, x="PC1", y="PC2")

#salida.to_tsv("compo.tsv" , sep = "\t")
