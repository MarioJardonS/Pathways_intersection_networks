#!/bin/bash

touch argumentos.sh

echo "#!/bin/bash" >> argumentos.sh
echo -n  "python core_network.py  network " >> argumentos.sh

for file in ../Networks/*.tsv ; 
do 
	echo -n "${file:12:-4} " >> argumentos.sh  ;
done

bash argumentos.sh
rm argumentos.sh
