#!/bin/bash
#script para generar archivos pem que contienen llaves desde 32 hasta 2048 bits
for bits in {32..2048}
do
	key=$(openssl genrsa -out key$bits.pem $bits &> /dev/null) #se crea una llave
	#modulus=$(openssl rsa -noout -modulus -in key.pem) #se obtiene el modulo de la llave
	#modulus=$(echo $modulus | cut -d "=" -f 2) #se corta el string para obtiene el hex del modulo
	#echo $modulus >> keys.txt
done
