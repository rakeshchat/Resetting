#!/bin/bash
#PBS -e /shokef/rakeshc/kcm1/aerr.txt
#PBS -o /shokef/rakeshc/kcm1/aout.txt

#=============================
#cd /shokef/
#cd ~/mycodes/
#icc kcm1.c
#icc kcm1.c -o job1
#qsub -q shokef or nano2 or shokefr or shokefmem rcjob.sh
#qstat -u rakeshc
#qdel <job number> e.g. 7317948
#=============================


cd /shokef/rakeshc/kcm1/

./job1

