# BATF TNFR2 regulatory T-cell state in ulcerative colitis

This repository contains the analysis code supporting the study of a BATF-associated, TNFR2-high, IL10-low or undetected regulatory T-cell state in ulcerative colitis and its relationship to treatment responses in experimental colitis.

## Data availability

- Newly generated mouse colon bulk RNA-seq: GEO **GSE349603**.
- Public human single-cell and spatial datasets: GSE214695, GSE231993, GSE282122, GSE189184, GSE242087 and GSE273194.
- Public mouse BATF ChIP-seq: GSE267882.
- Bulk-expression validation datasets: GSE87473 and GSE75214.
- Processed scATAC-seq T-cell object: Zenodo DOI **10.5281/zenodo.14165540**.

The GEO record GSE349603 is currently held private by GEO and should be released when the accession is quoted in a publicly available article, preprint or thesis.

## Repository structure

```text
analysis/
  01_scRNA_seq/          scRNA-seq integration, composition, trajectory, hdWGCNA and machine learning
  02_scATAC_ChIP/        scATAC-seq, chromVAR, footprinting and BATF ChIP-seq processing
  03_bulk_RNA_seq/       Il10-deficient mouse colon bulk RNA-seq and Figure 6 analyses
  04_multiplex_imaging/  multiplex fluorescence line-profile analysis
```

## Reproducibility notes

The scripts are organized according to the figure-level analysis workflow. Large expression objects, FASTQ files, alignment files, model objects and figure outputs are intentionally excluded from GitHub. Download public datasets from their primary repositories and place input files under a local `data/` directory. Some scripts originated from exploratory analyses and may require adjustment of filenames and object names to match the downloaded data.

Patient- or sample-level proportions were used as the statistical unit for cell-composition and inflammation-correlation analyses. The `GSE242087` material represents pouch and pre-pouch samples and was used as supportive epigenomic evidence rather than as a native-colon UC discovery cohort. The BATF ChIP-seq dataset is mouse-derived.

## Primary software

The analyses use R and Python packages including Seurat, Signac, ArchR, CellRank, scFates, CytoTRACE2, hdWGCNA, SCENIC, CellOracle, chromVAR, DESeq2, clusterProfiler and standard machine-learning packages. Exact package versions should be recovered from the original analysis environments where available.

## Citation

Please cite the associated manuscript and GEO accession GSE349603 when using this code. Manuscript citation details will be added after publication.

