import pathways_intersection_networks
import sys 


edges = pathways_intersection_networks.red(f"{sys.argv[1]}")

red = open(f"../Networks/{sys.argv[1]}.tsv" , "w")

for edge in edges:
	red.write(f"{edge[0]}\t{edge[1]}\t{edge[2]}\t{edge[3]}\t{edge[4]}\n")
	
red.close()
	
