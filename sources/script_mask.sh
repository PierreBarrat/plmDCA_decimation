#!/bin/bash

M=$2  # Length of sequence
NN=$3 # Half of M


    # cut to the first NN
    for i in ` seq 1 $M ` ; do
      awk -v i=$i 'BEGIN{FS=","}{if((NF==3)&&($1==i||$2==i))print $1" " $2" "$3}' $1 | sort -k 3 -g -r > p
      head -n ${NN} p >> J_L
      head -n 20 p >> J_20
    done

awk -v M=$M '{a[$1,$2]=1;a[$2,$1]=1;}END{for(i=1;i<=M;i++)for(j=i+1;j<=M;j++)if(a[i,j]!=1)print i,j}' J_L > mask_tLh.txt
awk -v M=$M '{a[$1,$2]=1;a[$2,$1]=1;}END{for(i=1;i<=M;i++)for(j=i+1;j<=M;j++)if(a[i,j]!=1)print i,j}' J_20 > mask_t20.txt

rm J_L J_20