import sys
import pandas as pd

redes = sys.argv[2:]
output = sys.argv[1]


def edges(red):
	red = pd.read_table(f"../Networks/{red}.tsv" )
	red = red.iloc[ :  , 0:2  ]
	aristas = set()
	for i in range(red.shape[0]):
		r = red.iloc[i , 0:2].tolist()
		r = tuple(r)
		aristas.add(r)

	return aristas


all_edges = edges(redes[0])
conteo = {}

for e in all_edges:
	conteo[e] = 1



for j in range(1,len(redes)):
	aristas = edges(redes[j])
	aristas_nuevas = aristas.difference(all_edges)
	aristas = aristas.intersection(all_edges)
	for e in aristas_nuevas:
		conteo[e] = 1
	for e in aristas:
		conteo[e] = conteo[e] + 1
	all_edges = all_edges.union(aristas_nuevas)

#	if len(core_network.intersection(aristas_nuevas)) > 0:
#		core_network = core_network.intersection(aristas_nuevas)
#	else:
#		print(j)
#		break
	

file = open(f"../Networks/{output}.tsv" , "w")
file2 = open(f"../Networks/edges.tsv" , "w")

for edge in conteo:
	file.write(f"{edge[0]}\t{edge[1]}\t{conteo[edge]}\n")
	file2.write(f"{edge[0]}\t{edge[1]}\n")

file.close()
file2.close()



