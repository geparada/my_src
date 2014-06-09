#!/bin/bash

while [ 1 ];
do
if ps -A | grep -e bowtie
then
       echo "Bowtie is Running"
else
       echo "Terminó el alineamiento de bowtie de los SJ con MapSplice"  | mail -s "Termino Bowtie" geparada88@gmail.com

       exit 0
fi
sleep 60
done
