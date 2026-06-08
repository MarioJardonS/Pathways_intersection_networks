import pandas as pd

data = pd.read_table("../Data/pathways_species2.tsv" )

#conjunto de rutas
pathways = data["Pathway"]
rutas = pathways.unique()
rutas = list(rutas)

pathways = pathways.tolist()

#funcion que da taxonomia (nivel especie) de ruta en una sola muestra
def taxa(ruta , muestra):
	r_m = []
	for i in range(len(pathways)):
		if pathways[i] == ruta and data.loc[i,muestra] > 0:
			r_m.append(i)
	 	
	duta_r_m = data.loc[r_m , "OTU"].unique()
	duta_r_m = set(duta_r_m)
	return duta_r_m
 	
 	
# Función para calcular índice de Jaccard
def jaccard(conjunto1, conjunto2):
    union = conjunto1.union(conjunto2)

    if len(union) == 0:
        return 0

    interseccion = conjunto1.intersection(conjunto2)

    return len(interseccion) / len(union)

 	
#funcion que crea red de rutas para una muestra

def red(muestra):
    edges = []

    for i in range(len(rutas)-1):
        for j in range(i+1, len(rutas)):

            t_i = taxa(rutas[i], muestra)
            t_j = taxa(rutas[j], muestra)

            inter = len(t_i.intersection(t_j))
            jac = jaccard(t_i, t_j)

            if inter > 0:
                edges.append([rutas[i], rutas[j], inter, jac])

    return edges
	
	
	
	

	

		
	
