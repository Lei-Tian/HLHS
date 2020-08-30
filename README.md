# HLHS
Single-Cell RNA-Seq Reveals Endocardial Defects in Hypoplastic Left Heart Syndrome
https://singlecell.broadinstitute.org/single_cell/study/SCP1020/single-cell-rna-seq-of-ipsc-derived-endocardium-endothelium-from-hypoplastic-left-heart-syndrome-patient#study-visualize


1) HLHS.R: combine all (Left Ventricle CDH5+, Left Ventricle CDH5-, Right Ventricle CDH5+, Right Ventricle CDH5-) cells of undeveloped left heart (ULH) fetal heart tissue;
2) d83.R: combine all (Left Ventricle CDH5+, Left Ventricle CDH5-, Right Ventricle CDH5+, Right Ventricle CDH5-) cells of normal fetal heart tissue;
3) d83_HLHS_LVp.R: combine Left Ventricle CDH5+ cells from ULH and normal fetal heart tissue;
4) d83_HLHS_RVp.R: combine Right Ventricle CDH5+ cells from ULH and normal fetal heart tissue;
5) iPSCEC.R: combine iPSC-ECs from both Hypoplastic left heart syndrome (HLHS) and control;
6) FunctionalEnrichmentAnalysis.R: functional enrichment analysis, the input is a gene list;
7) detectDEGs_bulkRNAseq: codes detect Differentially Expressed Genes in bulk RNA-seq data;
8) heatmap.R: plot heatmap in bulk RNA-seq data;
9) crosstalk.R: detect crosstalks between CM and Endocarduim.
10) app.R: code for web app
