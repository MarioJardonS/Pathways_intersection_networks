import pathways_intersection_networks
import sys 


edges = pathways_intersection_networks.red(f"{sys.argv[1]}")

red = open(f"{sys.argv[1]}.csv" , "w")

for edge in edges:
	red.write(f"{edge[0]},{edge[1]},{edge[2]}\n")
	
red.close()
	
