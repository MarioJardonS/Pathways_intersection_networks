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


core_network = edges(redes[0])

for j in range(1,len(redes)):
	aristas_nuevas = edges(redes[j])
	if len(core_network.intersection(aristas_nuevas)) > 0:
		core_network = core_network.intersection(aristas_nuevas)
	else:
		print(j)
		break 


file = open(f"../Networks/{output}.tsv" , "w")

for edge in core_network:
	file.write(f"{edge[0]}\t {edge[1]} \n")

file.close()




