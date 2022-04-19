# tumor specific marker discovery using sc/snRNA-seq data

## Requirements
- Seurat3 environment 
- list of seurat objects with cell type annotation 

## Running
```
usage: 
bash tumor_specific_marker.sh ${OUT_DIR} ${RDS_LIST} ${CELL_TYPE} {GTEX_TISSUE} {HPA_RNA_TISSUE} ${HPA_PROTEIN_TISSUE}
```

Example:
```
bash tumor_specific_marker.sh /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/test /diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/rds_list.txt Plasma lymphocytes bone_marrow,thymus,spleen,appendix,tonsil,lymph_node BM_and_lymphoid_tissues
```
