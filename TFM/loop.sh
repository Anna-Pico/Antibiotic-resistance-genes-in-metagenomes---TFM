#!/bin/bash

#recorre tots els fastqs de la carpeta i executa el script 
for arxiu in $(ls [0-9]*.fastq);do
	./executa_1_fitxer.sh $arxiu
done
