import sys

file  = open("argumentos.sh" , "w")
file.write("#!/bin/bash\n")

#argv[1] = grupo
argumentos = open("../Data/{sys.argv[1]}.txt" , "r")


#argv[2] = medida
i = 0 
for line in argumentos: 
	i = i + 1
	file.write(f"python print_network.py  ")
	file.write(line[0:(len(line)-1)])
	file.write("&\n")
	if  i == 40:
		file.write("wait\n")
		i = 0

file.close()
argumentos.close()	
