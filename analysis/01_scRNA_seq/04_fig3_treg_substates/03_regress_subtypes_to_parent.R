major_meta <- major_annotation_ibd_subtype@meta.data
major_meta$cell_name_gse <- rownames(major_meta)

t_meta <- IBD_T_TREG_SUB@meta.data[, c("final_subtype", "cell_name")]

merge_t <- left_join(major_meta, t_meta, by = "cell_name")

merge_t$final_subtype <- ifelse(
  is.na(merge_t$final_subtype),
  merge_t$cell_subtype,
  merge_t$final_subtype
)

rownames(merge_t) <- merge_t$cell_name_gse
merge_t <- merge_t[colnames(major_annotation_ibd_subtype), , drop = FALSE]
major_annotation_ibd_subtype2 <- major_annotation_ibd_subtype
major_annotation_ibd_subtype2@meta.data <- merge_t

DimPlot(major_annotation_ibd_subtype, group.by = "final_subtype")
major_annotation_ibd_subtype<- major_annotation_ibd_subtype2
major_annotation_ibd_subtype@meta.data$final_subtype[
  major_annotation_ibd_subtype@meta.data$final_subtype == "T/NK"
] <- "other"
saveRDS(major_annotation_ibd_subtype,file = 'major_annotation_ibd_subtype_tregback.rds')


output_h5ad_path <- "data\\IBD_all_cell_treg_sub.h5ad"

major_annotation_ibd_subtype$cell_subtype.x <- major_annotation_ibd_subtype$final_subtype
head(major_annotation_ibd_subtype@meta.data)
kidneyH5(major_annotation_ibd_subtype, output_h5ad_path)
output_h5ad_path <- "data\\treg_sub.h5ad"

major_annotation_ibd_subtype$cell_subtype.x <- major_annotation_ibd_subtype$final_subtype
head(major_annotation_ibd_subtype@meta.data)
kidneyH5(TREG_SUBSET, output_h5ad_path)
