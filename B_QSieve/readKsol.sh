#!/bin/bash	
while read -r line; do
	line=$(echo "$line" | tr '[:lower:]' '[:upper:]');
	#echo "$line"
	bin=$(echo "obase=2; ibase=16; $line" | bc | tr -d '\n');
	len=${#bin};
	dif=$((64-$len));#ceros que faltan
	for value in $(seq -w 1 $dif);#loop para completar los ceros
	do
		printf "0";
	done
	echo "$bin";
done < "$1"
