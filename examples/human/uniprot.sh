#!/usr/bin/env bash

# Human
wget -q https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz
zless ./UP000005640_9606.fasta.gz > ./human.fa
rm ./UP000005640_9606.fasta.gz

# JUN
wget -q https://rest.uniprot.org/uniprotkb/P05412.fasta
mv ./P05412.fasta ./JUN_HUMAN.fa
