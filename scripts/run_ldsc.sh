#!/bin/bash
conda activate ldsc

#These neeed to be comma separated strings
GWASES=$1
NS=$2
ANCESTRY=$3
OUTPUT_PREFIX=$4

LDSC_DATA=$DATA/ldsc
SUMSTATS_FILES=""

mkdir -p "$LDSC_DATA"
mkdir -p $(dirname $OUTPUT_PREFIX)

N_LIST=(${NS//,/ })
ANCESTRIES_LIST=(${ANCESTRIES//,/ })
i=0
for gwas in ${GWASES//,/ } ; do
  n=${N_LIST[i]}

  #python2.7 ~/ldsc/munge_sumstats.py \
  ./ldsc/munge_sumstats.py \
      --sumstats "${gwas}" \
      --N "$n" \
      --merge-alleles "$LDSC_DIR/w_hm3.snplist" \
      --snp RSID --a1 EA --a2 OA --p P --signed-sumstats BETA,0 \
      --chunksize 500000 \
      --out "$LDSC_DATA/${gwas}"

  SUMSTATS="${LDSC_DATA}/${gwas}.sumstats.gz"
  SUMSTATS_FILES+=($SUMSTATS)
  i=$((i+1))
done

SUMSTATS_STRING=${SUMSTATS_FILES[@]}
SUMSTATS_STRING=${SUMSTATS_STRING// /,}

if [[ $GWASES =~ "," ]]; then
  ./ldsc/ldsc.py \
    --rg "$SUMSTATS_STRING" \
    --ref-ld-chr "$LDSC_DIR"/1000genomes/ldscores/"$ANCESTRY"/ \
    --w-ld-chr "$LDSC_DIR"/1000genomes/ldscores/"$ANCESTRY"/ \
    --out $OUTPUT_PREFIX

    ldsc_rg_log="$OUTPUT_PREFIX".log
    genetic_correlation=$(awk '/p1   /' RS= "$ldsc_rg_log" | tr -s ' ' '\t')
    echo "$genetic_correlation" > "${OUTPUT_PREFIX}".tsv
else
  ./ldsc/ldsc.py \
    --h2 "$SUMSTATS" \
    --ref-ld-chr "$LDSC_DIR"/1000genomes/ldscores/"$ANCESTRY"/ \
    --w-ld-chr "$LDSC_DIR"/1000genomes/ldscores/"$ANCESTRY"/ \
    --out $OUTPUT_PREFIX
fi

