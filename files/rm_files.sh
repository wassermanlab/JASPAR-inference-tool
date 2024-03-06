#!/usr/bin/env bash

for TAXON in "fungi" "insects" "nematodes" "plants" "urochordates" "vertebrates"; do

    if [ -f .${TAXON}.uniaccs.pickle ]; then
        rm .${TAXON}.uniaccs.pickle
    fi

    if [ -f ${TAXON}.fa ]; then
        rm ${TAXON}.fa*
    fi

    if [ -f ${TAXON}.pfam.json ]; then
        rm ${TAXON}.pfam.json
    fi

    if [ -f ${TAXON}.profiles.json ]; then
        rm ${TAXON}.profiles.json
    fi

    if [ -f ${TAXON}.uniprot.json ]; then
        rm ${TAXON}.uniprot.json
    fi

    if [ -d ${TAXON} ]; then
        rm -r ${TAXON}
    fi

done
