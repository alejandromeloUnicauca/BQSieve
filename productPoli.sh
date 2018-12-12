#!/bin/bash	
rm salidap.txt
prod=1
filas=$(wc -l vec.txt | cut -d' ' -f1);
echo "$filas"
i=1
while [ "$i" -le 64 ];
do
	rm salida.txt pos.txt
	#echo "$i"
	while read f1 <&3 && read f2 <&2; do 
		f2=$(echo "$f2" | cut -c"$i");
		#echo "$f1 $f2";
		if [[ $f2 -eq 1 ]];then
			#prod=$(echo "$prod"*"$f1" | bc)
			echo $(cut -d';' -f2 <<< "$f1") >> salida.txt;
			echo $(cut -d';' -f1 <<< "$f1") >> pos.txt;
		fi
	done 3<polinomioB.txt 2<vec.txt
	./mulPoli "$1" >> salidap.txt;
	wait;
	i=$((i+1))
done
