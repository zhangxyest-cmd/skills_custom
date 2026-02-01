#!/usr/bin/env bash
set -euo pipefail

# Run DR_FU LOVO with full gene-based tests (no --joint)
# Outputs: results_gene/lovo/DR_FU_lovo_*_DR_FU.regenie.gz

jobs=${1:-/opt/notebooks/results_gene/lovo_jobs.tsv}
lovo_dir=${2:-/opt/notebooks/results_gene/lovo}
regenie=${3:-/opt/notebooks/regenie_v4.1.gz_x86_64_Linux_mkl}

mkdir -p "$lovo_dir"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

while IFS=$'\t' read -r trait setid mask aaf ID P chr outprefix completed; do
  if [ "$trait" = "trait" ]; then continue; fi
  if [ "$trait" != "DR_FU" ]; then continue; fi

  outbase="$lovo_dir/$outprefix"
  outcheck="${outbase}_DR_FU.regenie.gz"
  if [ -s "$outcheck" ]; then
    continue
  fi

  echo "[run] $trait $setid $mask $aaf chr${chr}"

  cond_file="/mnt/project/WGS_XY/conditional_list/DR_common_chr_pos_ref_alt.cleaned_chr${chr}.txt"
  cond_arg=()
  if [ -s "$cond_file" ]; then
    cond_arg=(--condition-list "$cond_file")
  else
    echo "[warn] empty or missing condition list for chr${chr}, skipping --condition-list"
  fi

  "$regenie" \
    --step 2 \
    --pred /opt/notebooks/dr_t2d_wgs_pc20_conditioned_pred.list \
    --bgen "/mnt/project/WGS_XY/bgen_t2d/ukb_wgs_chr${chr}_subset_corrected.bgen" \
    --sample "/mnt/project/WGS_XY/bgen_t2d/ukb_wgs_chr${chr}_subset_corrected.sample" \
    --ref-first \
    --keep "/mnt/project/PRS/Ancestry_Inference/cohort_EUR.txt" \
    --keep "/opt/notebooks/t2d_dr_phe_wgs_ver1126.txt" \
    --phenoFile "/opt/notebooks/t2d_dr_phe_wgs_ver1126.txt" \
    --phenoCol=DR_FU \
    --bt --firth --firth-se --approx \
    --vc-tests skat,skato,skato-acat,acatv,acato \
    --minMAC 20 \
    --covarFile "/opt/notebooks/t2d_dr_covar_wgs_ver1126.txt" \
    --covarCol=age,agesqr,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc13,pc14,pc15,pc16,pc17,pc18,pc19,pc20 \
    --catCovarList sex,wgs_provider,center \
    --maxCatLevels 30 \
    --anno-file "/mnt/project/WGS_XY/VEP_annotation/mask_ver3/final_map_transcript/chr${chr}_final_3col_map.v1120_VGT.tsv" \
    --mask-def /opt/notebooks/dmvas_6_full.masks \
    --set-list "/mnt/project/WGS_XY/VEP_annotation/mask_ver3/final_map_transcript_gene_set/gene_variant_set_list_chr${chr}.final.sorted.txt" \
    --aaf-bins 0.001 \
    --vc-maxAAF 0.001 \
    "${cond_arg[@]}" \
    --extract-setlist "$setid" \
    --mask-lovo "$setid,$mask,0.001" \
    --bsize 200 \
    --out "$outbase" --gz

done < "$jobs"

echo "DR_FU LOVO full-tests runs finished"
