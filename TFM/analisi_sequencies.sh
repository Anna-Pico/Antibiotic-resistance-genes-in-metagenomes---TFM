#!/bin/bash

##### Resistome analysis
#Carreguem els paquets que utilitzarem
echo "RESISTOME ANALYSIS"
echo "Loading tools (spack & module)..."
spack load diamond
spack load fastx-toolkit
module load metaxa/2.2.3

# Càrrega de la base de dades, només cal fer-ho una vegada per a obtenir l'arxiu .dmnd
# echo "Setting ARGminer-v1.1.1.A.fasta db..."
# diamond makedb --in ../ARGminer-v1.1.1.A_RA.fasta -d referenceARG

#Agafem un arxiu .La variable arxiu té l'extensió fastq
arxiu=$1;

echo "carregant arxiu $arxiu...";

#arxiu fastq amb "qual" davant
arxiu_qual="qual_$arxiu"

#Filtre de qualitat: fastq_quality_filter -Q33 -q20 -p80
echo "Quality filter..."
fastq_quality_filter -Q33 -q20 -p80 -i "$arxiu" -o "$arxiu_qual"

#arxiu fastq amb "qual" davant
arxiu_fasta="$arxiu_qual.fasta"

# Fastq a fasta: fastq_to_fasta -Q33 -i
echo ".fastq to .fasta..."
fastq_to_fasta -Q33 -i "$arxiu_qual" -o "$arxiu_fasta"

# DIAMOND 

# arxiu daa de sortida
arxiu_daa="$arxiu_fasta.daa"

#blastx: diamond blastx -d database -q input -a output --query-cover 80 --id 90 -k1 -c1 -p16
echo "DIAMOND blastx in progress..."
diamond blastx -d referenceARG.dmnd -q "$arxiu_fasta" -a "$arxiu_daa" --query-cover 80 --id 90 -k1 -c1 -p16

#arxiu m8
arxiu_m8="$arxiu_daa.m8"

#Visualització dels resultats: diamond view -a input -o output
echo ".daa to .m8..."
diamond view -a "$arxiu_daa" -o "$arxiu_m8"

# METAXA (16s analysis)

#crear una carpeta
arxiu_sense_extensio=$(echo $arxiu | sed "s/\.fastq//")
nom_carpeta="metaxa-$arxiu_sense_extensio-output"

#elimina la carpeta si ja existeix
rm -rf "$nom_carpeta" > /dev/null

echo "Creating $nom_carpeta folder..."

mkdir "$nom_carpeta"
mv "$arxiu_fasta" "$nom_carpeta"
cd "$nom_carpeta"

# blastn metaxa: metaxa2 -i input -o output
nom_output="output-$arxiu_sense_extensio"
echo $nom_output
echo "Metaxa blastn in progress..."
metaxa2 -i "$arxiu_fasta" -o "$nom_output"
cp "output-$arxiu_sense_extensio.summary.txt" ../
mv "$arxiu_fasta" ../
cd ..

# comprimir arxiu fasta
gzip "$arxiu_fasta"

bacteria=$(cat "$nom_output".summary.txt | grep Bacteria | sed -n 2p | awk '{print $2}')
echo "$bacteria sequences are related to bacterial 16S"

# Comptar nombre de resistències per cada família d'antibiòtics
res_count="res_count-$arxiu_sense_extensio.csv"
gawk '{print $2;}' "$arxiu_m8" | cut -d "|" -f1 | sort | uniq -c | column -t > "$res_count"

# Afegir capçalera
resistome="res-$arxiu_sense_extensio.csv"
echo -e 'count   family'|column -t > "$resistome" && cat "$res_count" >> "$resistome"




