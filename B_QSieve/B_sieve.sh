#!/bin/bash	
if [[ "$1" == "-h" ]];then
	num=$(echo "ibase=16;obase=10; $2" | bc);
else
	num=$2;
fi
./B_QSieve "$1" "$num" "$3" "$4"
wait
exit_status=$?;
#echo $exit_status
wait;
if [[ "$exit_status" == "1" ]];then
	exit 1;
fi
/home/debian/Descargas/cado-nfs/build/debian12/linalg/bwc/mf_scan --ascii-in mfile=matrix.txt --binary-out ofile=matrix.bin --freq
wait
mkdir /tmp/bqsieve
cp matrix.bin matrix.rw.bin matrix.cw.bin /tmp/bqsieve/.
wait	
echo "Solucionando Matriz...";
/home/debian/Descargas/cado-nfs/build/debian12/linalg/bwc/bwc.pl :complete thr=2 m=64 n=64 nullspace=left interval=100 matrix=/tmp/bqsieve/matrix.bin wdir=/tmp/bqsieve/ interleaving=0 > outputbwc.txt
wait
cp /tmp/bqsieve/K.sols0-64.0.txt /home/debian/prsa/B_QSieve/.
wait
./readKsol.sh K.sols0-64.0.txt > vec.txt
wait
./productPoli.sh "$num"
wait
cat salidap.txt
wait
rm matrix.bin matrix.rw.bin matrix.cw.bin matrix.txt residuos.txt outputbwc.txt K.sols0-64.0.txt salidap.txt vec.txt pos.txt salida.txt polinomio.txt
