#!/bin/bash

conda activate ldsc

#TODO: make sure these are arrays
FILES=$1
ANCESTRIES=$2
N=$3

LDSC_RESULTS=$RESULTS/ldsc
LDSC_DATA=$DATA/ldsc
SUMSTATS_FILES=""
mkdir -p "$LDSC_RESULTS"
mkdir -p "$LDSC_DATA"

for (( i = 0; i < ${#FILES[@]}; ++i )); do
  file=${FILES[i]}
  ancestry=${ANCESTRIES[i]}
  n=${N[i]}

  #python2.7 ~/ldsc/munge_sumstats.py \
  ./ldsc/munge_sumstats.py \
      --sumstats "${file}" \
      --N "$n" \
      --merge-alleles "$LDSC_DIR/w_hm3.snplist" \
      --snp RSID --a1 EA --a2 OA --p P --signed-sumstats BETA,0 \
      --chunksize 500000 \
      --out "$LDSC_DATA/${file}"

  SUMSTATS_FILES+=("$LDSC_DATA/${file}.sumstats.gz")
done

SUMSTATS_STRING=${SUMSTATS_FILES[@]}
SUMSTATS_STRING=${SUMSTATS_STRING// /,}

#TODO: ok, what are we actually trying to compare here?  Between GWASes of the same ancestry?
for ancestry in eur his afr
do
    #python2.7 ~/ldsc/ldsc.py \
    ./ldsc/ldsc.py \
        #--rg "$LDSC_RESULTS"/first_ais_${ancestry}.sumstats.gz,"$LDSC_RESULTS"/subsequent_ais_${ancestry}.sumstats.gz,$LDSC_RESULTS/subsequent_mace_${ancestry}.sumstats.gz \
        --rg "$SUMSTATS_STRING"
        --ref-ld-chr $LDSC_DIR/1000genomes/ldscores/$ancestry/ \
        --w-ld-chr $LDSC_DIR/1000genomes/ldscores/$ancestry/ \
        --out $LDSC_RESULTS/ais_${ancestry}_rg_results
done



#function parse_rg_file() {
    #h2_slope=$(grep "Total Observed scale h2" $LDSC_RESULTS/subsequent_ais_${ancestry}_ldsc.log | cut -d' ' -f5)
    #h2_se=$(grep "Total Observed scale h2" $LDSC_RESULTS/subsequent_ais_${ancestry}_ldsc.log | cut -d' ' -f6)
    #gc=$(grep "Lambda GC" $LDSC_RESULTS/subsequent_ais_${ancestry}_ldsc.log | cut -d' ' -f3)

    #TODO: do this later.  Think about how you might want to do this in the future

#}
