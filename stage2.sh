STAGE2_INPUT_DIR=/GLOBALFS/sysu_hpcscc_2/zhj/projects/ASC25_simulation/rna_2/m5C-UBSseq/workspace/internal_files/count_sites
STAGE2_OUTPUT_DIR=./stage2_output
STAGE2_BIN_DIR=./bin
STAGE2_CASE1=SRR23538290
STAGE2_CASE2=SRR23538291
STAGE2_CASE3=SRR23538292

rm -rf ${STAGE2_OUTPUT_DIR}
mkdir -p ${STAGE2_OUTPUT_DIR}

# echo ${STAGE2_INPUT_DIR}/${STAGE2_CASE1}.genome.arrow
echo ${STAGE2_INPUT_DIR}/${STAGE2_CASE1}.genome.arrow

python ${STAGE2_BIN_DIR}/group_pileup.py \
    -i ${STAGE2_INPUT_DIR}/${STAGE2_CASE1}.genome.arrow \
       ${STAGE2_INPUT_DIR}/${STAGE2_CASE2}.genome.arrow \
       ${STAGE2_INPUT_DIR}/${STAGE2_CASE3}.genome.arrow \
    -o ${STAGE2_OUTPUT_DIR}/WT.arrow

python ${STAGE2_BIN_DIR}/select_sites.py \
    -i ${STAGE2_OUTPUT_DIR}/WT.arrow \
    -o ${STAGE2_OUTPUT_DIR}/WT.prefilter.tsv

python ${STAGE2_BIN_DIR}/filter_sites.py \
    -i  ${STAGE2_INPUT_DIR}/${STAGE2_CASE1}.genome.arrow \
    -m ${STAGE2_OUTPUT_DIR}/WT.prefilter.tsv \
    -b ${STAGE2_OUTPUT_DIR}/${STAGE2_CASE1}.bg.tsv \
    -o ${STAGE2_OUTPUT_DIR}/${STAGE2_CASE1}.filtered.tsv

python ${STAGE2_BIN_DIR}/filter_sites.py \
    -i  ${STAGE2_INPUT_DIR}/${STAGE2_CASE2}.genome.arrow \
    -m ${STAGE2_OUTPUT_DIR}/WT.prefilter.tsv \
    -b ${STAGE2_OUTPUT_DIR}/${STAGE2_CASE2}.bg.tsv \
    -o ${STAGE2_OUTPUT_DIR}/${STAGE2_CASE2}.filtered.tsv

python ${STAGE2_BIN_DIR}/filter_sites.py \
    -i  ${STAGE2_INPUT_DIR}/${STAGE2_CASE3}.genome.arrow \
    -m ${STAGE2_OUTPUT_DIR}/WT.prefilter.tsv \
    -b ${STAGE2_OUTPUT_DIR}/${STAGE2_CASE3}.bg.tsv \
    -o ${STAGE2_OUTPUT_DIR}/${STAGE2_CASE3}.filtered.tsv