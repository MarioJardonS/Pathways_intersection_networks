import pathways_intersection_networks
import pandas as pd

muestras = pathways_intersection_networks.samples
funciones = pathways_intersection_networks.rutas


matriz = []


for i in range(len(muestras)):
	s_i = []
	for j in range(len(funciones)):
		c = pathways_intersection_networks.comunidad(funciones[j] , muestras[i])
		s_i.append(pathways_intersection_networks.shan_index(c))
		
	matriz.append(s_i)
	
	
data = pd.DataFrame(matriz)

data.columns = funciones
#data.rows = muestras 

data.to_csv("shannon_index.tsv" , sep = "\t") 
