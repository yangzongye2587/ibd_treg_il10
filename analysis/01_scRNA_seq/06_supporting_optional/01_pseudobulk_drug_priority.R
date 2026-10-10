library(Seurat)
library(dplyr)
library(Matrix)
library(edgeR)
library(tibble)
library(tidyr)
setwd("data")  
obj <- TREG_SUBSET
Idents(obj) <- "treg_sub"

sample_col <- "orig.ident"
group_col  <- "treg_sub"

obj2 <- subset(obj, subset = treg_sub %in% c("TNF_TREG", "NAIVE_LIKE_TREG"))

meta <- obj2@meta.data %>% rownames_to_column("cell")
counts <- GetAssayData(obj2, assay = "RNA", slot = "counts")

meta2 <- meta %>%
  dplyr::select(cell, all_of(sample_col), all_of(group_col)) %>%
  dplyr::rename(sample = all_of(sample_col), group = all_of(group_col)) %>%
  dplyr::filter(group %in% c("TNF_TREG", "NAIVE_LIKE_TREG"))

cell_stat <- meta2 %>% dplyr::count(sample, group, name = "n_cells")
keep_sample_group <- cell_stat %>% dplyr::filter(n_cells >= 20)

meta2 <- meta2 %>%
  dplyr::inner_join(keep_sample_group, by = c("sample", "group"))

pb_list <- split(meta2$cell, paste(meta2$sample, meta2$group, sep = "__"))
pb_mat <- sapply(pb_list, function(cells) Matrix::rowSums(counts[, cells, drop = FALSE]))
pb_mat <- as.matrix(pb_mat)

pb_meta <- tibble(pb_id = colnames(pb_mat)) %>%
  tidyr::separate(pb_id, into = c("sample", "group"), sep = "__", remove = FALSE)

pb_meta$group <- factor(pb_meta$group, levels = c("NAIVE_LIKE_TREG", "TNF_TREG"))

dge <- DGEList(counts = pb_mat, samples = pb_meta)
keep <- filterByExpr(dge, group = pb_meta$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ group, data = pb_meta)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)

deg <- topTags(qlf, n = Inf)$table %>%
  rownames_to_column("gene") %>%
  as_tibble()

write.csv(deg, "TNF_TREG_vs_NAIVE_pseudobulk_DEG.csv", row.names = FALSE)

up_genes <- deg %>% dplyr::filter(logFC > 0.5, FDR < 0.05) %>% dplyr::pull(gene)
down_genes <- deg %>% dplyr::filter(logFC < -0.5, FDR < 0.05) %>% dplyr::pull(gene)

write.table(up_genes, "TNF_TREG_up_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(down_genes, "TNF_TREG_down_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
dir.create("scDrugPrio_TNF_TREG", showWarnings = FALSE)
dir.create("scDrugPrio_TNF_TREG/network_distances", showWarnings = FALSE)
dir.create("scDrugPrio_TNF_TREG/translated_DEGs", showWarnings = FALSE)
dir.create("scDrugPrio_TNF_TREG/FC_criteria_checking", showWarnings = FALSE)
dir.create("scDrugPrio_TNF_TREG/intracellular_centrality", showWarnings = FALSE)
dir.create("scDrugPrio_TNF_TREG/Final_ranking", showWarnings = FALSE)
# install.packages("BiocManager")
# BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"))
# install.packages("devtools")
# devtools::install_github("SDTC-CPMed/scDrugPrio")

library(scDrugPrio)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
deg_use <- deg %>%
  dplyr::filter(!is.na(gene), !is.na(logFC), !is.na(FDR)) %>%
  dplyr::mutate(
    ENTREZID = AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = gene,
      keytype = "SYMBOL",
      column = "ENTREZID",
      multiVals = "first"
    )
  ) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

head(deg_use)
disease_genes_entrez <- as.numeric(unique(deg_use$ENTREZID))
length(disease_genes_entrez)
deg_for_scdrugprio <- deg_use %>%
  dplyr::transmute(
    gene = as.character(ENTREZID),
    p_val_adj = FDR,
    avg_log2FC = logFC
  )

write.table(
  deg_for_scdrugprio,
  file = "scDrugPrio_TNF_TREG/translated_DEGs/Cluster_TNF_TREG.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
load_sample_data()
ppin <- ppin_formatting(lit_ppi)
dim(ppin)
head(ppin)
drug_bank_use <- drug_bank_example_data

drug_bank_use <- drug_bank_use[grepl("Humans", drug_bank_use[,8]), ]
drug_bank_use <- drug_bank_use[grepl("approved", drug_bank_use[,3]), ]
drug_bank_use <- drug_bank_use[!is.na(drug_bank_use[,6]), ]

# vignette 里这里把 drug target 转成 Entrez ID
# 如果你的对象列位置和官方一致，直接照抄
drug_bank_use[,6] <- translation_mouse_human[
  match(drug_bank_use[,6], translation_mouse_human[,1]), 2
]
drug_bank_use <- drug_bank_use[!is.na(drug_bank_use[,6]), ]

drug_target_matrix <- create_drug_target_matrix(
  drugID = drug_bank_use[,1],
  target = as.numeric(drug_bank_use[,6])
)

drug_target_matrix <- prepare_drug_target_matrix_for_network_distance_calculation(
  ppin,
  drug_target_matrix,
  file_name = "TNF_TREG_lit_ppin",
  out_dir = "scDrugPrio_TNF_TREG"
)

same_drugs <- read.table(
  "scDrugPrio_TNF_TREG/SAME_DRUG_TARGETS_TNF_TREG_lit_ppin.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)
average_closest_distance_network_drug_screening(
  ppin = ppin,
  drug_target_matrix = drug_target_matrix,
  disease_genes = disease_genes_entrez,
  file_name = "Cluster_TNF_TREG",
  cores = 4,   # 你自己改
  out_dir = "scDrugPrio_TNF_TREG/network_distances"
)
lf <- list.files(
  path = "scDrugPrio_TNF_TREG/network_distances",
  pattern = "drug-disease_closest_distances_vs_random_bin_adjusted__"
)

separate_unique_drug_target_combinations_into_individual_drugs(
  files = lf,
  remove_part_of_file_name_for_output = "drug-disease_closest_distances_vs_random_bin_adjusted__",
  same_drugs = same_drugs,
  disease_specific_drugs = NULL,
  output_file_name_add_on = "INDIVIDUAL_DRUGS_",
  in_dir = "scDrugPrio_TNF_TREG/network_distances",
  out_dir = "scDrugPrio_TNF_TREG/network_distances"
)
drug_dists <- list.files(
  "scDrugPrio_TNF_TREG/network_distances",
  pattern = "drug-disease_closest_distance",
  full.names = TRUE
)

fc_evaluation <- pharma_effect_on_fold_change(
  drug_dist_files = drug_dists,
  deg_files = "scDrugPrio_TNF_TREG/translated_DEGs/Cluster_TNF_TREG.txt",
  pharma_effect = drug_bank_use[, c(1,6:7)],
  save_names = "Cluster_TNF_TREG",
  out_dir = "scDrugPrio_TNF_TREG/FC_criteria_checking"
)
fc_evaluation_filt <- fc_evaluation
fc_evaluation_filt <- fc_evaluation_filt[as.numeric(fc_evaluation_filt[,3]) < 1, ]
fc_evaluation_filt <- fc_evaluation_filt[as.numeric(fc_evaluation_filt[,7]) < 0.05, ]
fc_evaluation_filt <- fc_evaluation_filt[order(fc_evaluation_filt[,1]), ]

fc_evaluation_filt <- cbind(
  Drug_name = drug_bank_use[match(fc_evaluation_filt[,1], drug_bank_use[,1]), 2],
  fc_evaluation_filt
)

write.table(
  fc_evaluation_filt,
  file = "scDrugPrio_TNF_TREG/FC_criteria_checking/SUMMARY_only_drugs_passing_network_criteria.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
degs_mat <- matrix(
  as.numeric(unique(deg_use$ENTREZID)),
  ncol = 1
)
colnames(degs_mat) <- "Cluster_TNF_TREG"

intra_cent <- intracellular_drug_centrality(
  ppin = ppin,
  drug_target_matrix = drug_target_matrix,
  degs = degs_mat,
  file_name = "intracellular_centrality_TNF_TREG",
  centrality_alg = "eigenvector centralities",
  out_dir = "scDrugPrio_TNF_TREG/intracellular_centrality"
)

head(intra_cent)
# 你手工补完后的表，读回来
fc_evaluation_done <- read.table(
  "scDrugPrio_TNF_TREG/FC_criteria_checking/SUMMARY_only_drugs_passing_network_criteria_MANUAL_DONE.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

inter_cent_one <- c(Cluster_TNF_TREG = 1)

drug_rank <- final_drug_prioritization(
  fc_evaluation = fc_evaluation_done,
  pos_DrugID = 2,
  pos_clusterID = 19,
  keep = c(1:4, 18),
  inter_cent = inter_cent_one,
  intra_cent = intra_cent[, ncol(intra_cent)],
  file_name = "FINAL_drug_ranking_TNF_TREG",
  out_dir = "scDrugPrio_TNF_TREG/Final_ranking"
)

head(drug_rank)