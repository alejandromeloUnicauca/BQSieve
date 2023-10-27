#!/bin/bash	
./prsa "$1"
wait
/home/ingesis/cado-nfs-2.3.0/cado-nfs-2.3.0/build/debian9/linalg/bwc/mf_scan --ascii-in mfile=matrix.txt --binary-out ofile=matrix.bin --freq
wait
cp matrix.bin matrix.rw.bin matrix.cw.bin /tmp/.
wait
rm -r /tmp/c60/
wait
/home/ingesis/cado-nfs-2.3.0/cado-nfs-2.3.0/build/debian9/linalg/bwc/bwc.pl :complete thr=2 m=64 n=64 nullspace=left interval=100 matrix=/tmp/matrix.bin wdir=/tmp/c60/c60.bwc interleaving=0
wait
cp /tmp/c60/c60.bwc/K.sols0-64.0.txt /home/ingesis/prsa/.
wait
./readKsol.sh K.sols0-64.0.txt > vec.txt
wait
./productPoli.sh "$1"




