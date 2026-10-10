library(dplyr)
library(stringr)

meta <- sce.all.filt@meta.data

meta <- meta %>%
  mutate(
    inflammation_score = as.numeric(
      str_extract(characteristics_ch1.14, "[0-9.]+")
    )
  )
trm_prop <- meta %>%
  filter(marjor_type == "Immune", celltype_immune == "T/NK") %>%
  group_by(orig.ident) %>%
  summarise(
    total_T = n(),
    TRM_CD8_n = sum(T_subtype == 'TH17'),
    TRM_CD8_prop = TRM_CD8_n / total_T,
    inflammation_score = unique(inflammation_score)
  )
stopifnot(length(unique(trm_prop$inflammation_score)) ==
            nrow(trm_prop))
cor.test(
  trm_prop$TRM_CD8_prop,
  trm_prop$inflammation_score,
  method = "spearman"
)
library(ggplot2)

ggplot(trm_prop, aes(TRM_CD8_prop, inflammation_score)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Proportion of TRM_CD8 cells",
    y = "Inflammation score",
    title = "Correlation between TRM_CD8 proportion and inflammation"
  )
