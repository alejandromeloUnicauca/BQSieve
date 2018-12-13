#!/bin/bash
#script para generar archivos pem que contienen llaves desde 32 hasta 2048 bits
for bits in {32..2048}
do
	modulus=$(openssl rsa -noout -modulus -in keys/key$bits.pem) #se obtiene el modulo de la llave
	modulus=$(echo $modulus | cut -d "=" -f 2) #se corta el string para obtiene el hex del modulo
	key=$(echo "obase=10; $modulus" | bc)
	./prsaBlocks $key 10 >> results/result$bits.txt
	wait
done
