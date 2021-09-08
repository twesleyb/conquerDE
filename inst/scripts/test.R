#!/usr/bin/env Rscript

root <- "~/projects/conquerDE"
devtools::load_all(root)

library(dplyr)
library(data.table)

# load the data
data(meta)
data(gene_map)
data(filt_counts)

# mean-reference model (intercept term embodies first level of condition, CTRL)
model_matrix <- model.matrix(~condition,data=meta)

# edgeR analysis
edgeR_data <- edgeR::DGEList(filt_counts)
edgeR_data <- edgeR::calcNormFactors(edgeR_data)
edgeR_data <- edgeR::estimateDisp(edgeR_data, model_matrix)
edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
edgeR_res <- edgeR::glmQLFTest(edgeR_fit)
edgeR_res <- edgeR::topTags(edgeR_res, n=Inf)
edgeR_df <- edgeR_res$table %>% 
		as.data.table(keep.rownames="ensembl") %>% 
		left_join(gene_map, by="ensembl") %>% 
		setcolorder(c("ensembl", "symbol", "entrez")) %>% 
		mutate(candidate=FDR<0.05)

sum(edgeR_df$candidate)

glmGamPoi_fit <- glmGamPoi::glm_gp(filt_counts, design = model_matrix, on_disk = FALSE)
glmGamPoi_res <- glmGamPoi::test_de(glmGamPoi_fit, contrast='conditionIgAN')
glmGamPoi_df <- glmGamPoi_res %>% rename(ensembl=name) %>%
		left_join(gene_map, by="ensembl") %>% 
		setcolorder(c("ensembl", "symbol", "entrez")) %>% 
		mutate(candidate=adj_pval<0.05) %>%
		arrange(adj_pval)

sum(glmGamPoi_df$candidate)

# compare results
df_a <- glmGamPoi_df %>% mutate(x=lfc * log(pval)) %>% select(ensembl,x)
df_b <- edgeR_df %>% mutate(y=logFC * log(PValue)) %>% select(ensembl,y)
dm <- df_a %>% left_join(df_b, by="ensembl") %>% as.data.table() %>% as.matrix(rownames="ensembl")
cor(dm)

# benchmark to compare timing
bench::mark(
			edgeR={
					edgeR_data <- edgeR::DGEList(filt_counts)
					edgeR_data <- edgeR::calcNormFactors(edgeR_data)
					edgeR_data <- edgeR::estimateDisp(edgeR_data, model_matrix)
					edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
					edgeR_res <- edgeR::glmQLFTest(edgeR_fit)
					edgeR_res <- edgeR::topTags(edgeR_res, n=Inf)
			},
			glmGamPoi={
					glmGamPoi_fit <- glmGamPoi::glm_gp(filt_counts, design = model_matrix, on_disk = FALSE)
					glmGamPoi_res <- glmGamPoi::test_de(glmGamPoi_fit, contrast='conditionIgAN')
			},
			check=FALSE, min_iterations=3
			)
