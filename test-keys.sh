#!/bin/bash
#script para generar archivos pem que contienen llaves desde 32 hasta 2048 bits
cont=1
rm results/result.log
for bits in {72..128..8}
do
    echo "Prueba N°$cont de $bits bits ">>results/result.log
	modulus=$(openssl rsa -noout -modulus -in keys/key$bits.pem) #se obtiene el modulo de la llave
	modulus=$(echo $modulus | cut -d "=" -f 2) #se corta el string para obtiene el hex del modulo
	key=$(echo "obase=10; ibase=16; $modulus" | bc)
	echo "N=$key">>results/result.log
	echo "Digitos de N:${#key}">>results/result.log
    echo "Prueba con bloques">>results/result.log
	/usr/bin/time -p -o results/result.log -a ./facB.sh $key 20 > results/output.txt 
	grep Bloques results/output.txt >> results/result.log
	echo "" >> results/result.log
	wait
	#echo "Prueba con divisiones triviales">>results/result.log
    #/usr/bin/time -p -o results/result.log -a ./facT.sh $key > results/output.txt 
    #wait
    ((cont++))
    echo "" >> results/result.log
done
