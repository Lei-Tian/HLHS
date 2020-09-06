• Project Title: Single-Cell RNA-Seq Reveals Endocardial Defects in Hypoplastic Left Heart Syndrome

• Project Overview: https://github.com/Lei-Tian/HLHS/wiki/Project-Overview

• Web apps for data visulization:

  1) https://singlecell.broadinstitute.org/single_cell/study/SCP1020/single-cell-rna-seq-of-ipsc-derived-endocardium-endothelium-from-hypoplastic-left-heart-syndrome-patient#study-visualize
  2) https://singlecell.broadinstitute.org/single_cell/study/SCP1021/single-cell-rna-seq-of-normal-human-fetal-heart#study-visualize

• Code:
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

Intrinsic Endocardial Defects Contribute to Hypoplastic Left Heart Syndrome (https://doi.org/10.1016/j.stem.2020.07.015)

Yifei Miao#, Lei Tian#, Marcy Martin, Sharon L. Paige, Francisco X. Galdos, Jibiao Li, Alyssa Klein, Hao Zhang, Ning Ma, Yuning Wei, Maria Stewart, Soah Lee, Jan-Renier Moonen, Bing Zhang, Paul Grossfeld, Seema Mital, David Chitayat, Joseph C. Wu, Marlene Rabinovitch, Timothy J. Nelson, Shuyi Nie, Sean M. Wu, Mingxia Gu

#These authors contributed equally
Contact person: Mingxia Gu, Mingxia.Gu@cchmc.org
