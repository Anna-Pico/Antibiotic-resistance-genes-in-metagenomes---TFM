#!/bin/bash

#recorre tots els fastqs de la carpeta
#per fer analisi 
for arxiu in $(ls [0-9]*.fastq);do
	./executa_1_fitxer.sh $arxiu
done
