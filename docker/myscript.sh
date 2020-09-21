read -r -d '' USAGE <<'EOF'
Perform differential expression analysis between cell type of interest and the rest of cells
Usage:
  myscript.sh [options] ARGS

Options:
-h: Print this help message
-O OUT_DIR: required, path to output folder
-R RDS_LIST: required, path of a file with the list of sample names and rds pathspath of a file with the list of sample names and rds paths
-C CELL_TYPE: required, cell type of interest
-G GTEX_TT: required, tissue type of interest on GTEX
-H HPA_RNA_TT: required, tissue type(s) of interest on Human Protein Atlas RNA data
-P HPA_PRO_TT: required, tissue type of interest on Human Protein Atlas Protein data

RDS_LIST is a text file with 2 columns. First column is the path to rds. Second column is the sample id.

EOF

function confirm {
    FN=$1
    WARN=$2
    if [ ! -s $FN ]; then
        if [ -z $WARN ]; then
            >&2 echo ERROR: $FN does not exist or is empty
            exit 1
        else
            >&2 echo WARNING: $FN does not exist or is empty.  Continuing
        fi
    fi
}

function test_exit_status {
    # Evaluate return value for chain of pipes; see https://stackoverflow.com/questions/90418/exit-shell-script-based-on-process-exit-code
    rcs=${PIPESTATUS[*]};
    for rc in ${rcs}; do
        if [[ $rc != 0 ]]; then
            >&2 echo Fatal error.  Exiting.
            exit $rc;
        fi;
    done
}

while getopts ":h:O:R:C:G:H:P:" opt; do
  case $opt in
    h)
      echo "$USAGE"
      exit 0
      ;;
    O) 
      OUT_DIR=$OPTARG
      ;;
    R) 
      RDS_LIST=$OPTARG
      confirm $RDS_LIST
      ;;
    C)
      CELL_TYPE=$OPTARG
      ;;
    G)
      GTEX_TT=$OPTARG
      ;;
    H)
      HPA_RNA_TT=$OPTARG
      ;;
    P)
      HPA_PRO_TT=$OPTARG
      ;;
    \?)
      >&2 echo "Invalid option: -$OPTARG" 
      >&2 echo "$USAGE"
      exit 1
      ;;
    :)
      >&2 echo "Option -$OPTARG requires an argument." 
      >&2 echo "$USAGE"
      exit 1
      ;;
  esac
done
shift $((OPTIND-1))

CMD="Rscript tumor_specific_markers_V1.0.R -r ${RDS_LIST} -o ${OUT_DIR} -c ${CELL_TYPE}  -s ./DOWNLOAD_DB/Cell_Surface_Protein_Atlas_S2_File.txt -p ./DOWNLOAD_DB/Human_Protein_Atlas_subcellular_location.txt"
>&2 echo Running $CMD
eval $CMD

CMD="perl step0.get_surface_genes_metadata.pl > ${OUT_DIR}/GTEX_metadata_subset_tumor_cells_DEG.txt"
>&2 echo Running $CMD
eval $CMD

CMD="perl step1.get_surface_genes_expression_data.pl ${OUT_DIR}/${CELL_TYPE}_specific_DEG_with_surface_annotations_from_3DB_gene_list.txt > ${OUT_DIR}/GTEX_expression_matrix_subset_tumor_cells_DEG.txt"
>&2 echo Running $CMD
eval $CMD

CMD="grep captured ${OUT_DIR}/GTEX_expression_matrix_subset_tumor_cells_DEG.txt > ${OUT_DIR}/not_captured_GTEX_DEG.txt"
>&2 echo Running $CMD
eval $CMD

CMD="sed -e s/#//g -i ${OUT_DIR}/not_captured_GTEX_DEG.txt"
>&2 echo Running $CMD
eval $CMD

CMD="perl step2.analyze_expression_change.pl ${OUT_DIR} > ${OUT_DIR}/GTEX_tissue_specific_DEG.txt"
>&2 echo Running $CMD
eval $CMD

CMD="perl HPA_RNA_find_outlier_expressions_t_test.pl ${OUT_DIR}/${CELL_TYPE}_specific_DEG_with_surface_annotations_from_3DB.txt DOWNLOAD_DB/rna_tissue_consensus.tsv > ${OUT_DIR}/HPA_RNA_tissue_specific_DEG.txt"
>&2 echo Running $CMD
eval $CMD

CMD="python3 HPA_Protein_filtering.py -p DOWNLOAD_DB/HPA_normal_tissue.tsv -o ${OUT_DIR} -t DOWNLOAD_DB/HPA_Tissue_type_matching.txt -g ${OUT_DIR}/${CELL_TYPE}_specific_DEG_with_surface_annotations_from_3DB.txt"
>&2 echo Running $CMD
eval $CMD

CMD="Rscript adding_tissue_specificity_annotations.R -o ${OUT_DIR} -c ${CELL_TYPE} -g ${GTEX_TT} -r ${HPA_RNA_TT} -p ${HPA_PRO_TT}"
>&2 echo Running $CMD
eval $CMD

>&2 echo COMPLETED SUCCESSFULLY
