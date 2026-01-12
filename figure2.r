#figure 2a
#构建microeco 对象
safe_microtable <- function(otu_table, tax_table, sample_table) {
  # 1. 确保样本 ID 匹配
  common_samples <- intersect(colnames(otu_table), rownames(sample_table))
  otu_table <- otu_table[, common_samples, drop = FALSE]
  sample_table <- sample_table[common_samples, , drop = FALSE]
  
  # 2. 确保 taxa ID 匹配
  common_taxa <- intersect(rownames(otu_table), rownames(tax_table))
  otu_table <- otu_table[common_taxa, , drop = FALSE]
  tax_table <- tax_table[common_taxa, , drop = FALSE]
  
  # 3. 移除所有 abundance 全为 0 的 OTU
  nonzero_taxa <- rowSums(otu_table) > 0
  otu_table <- otu_table[nonzero_taxa, , drop = FALSE]
  tax_table <- tax_table[rownames(otu_table), , drop = FALSE]
  
  # 4. 构建 microtable 对象
  mt <- microtable$new(
    sample_table = sample_table,
    otu_table = otu_table,
    tax_table = tax_table
  )
  
  return(mt)
}

# Create new microtable object with matched data
dataset_virus_aggregated <- safe_microtable(
  otu_table = species_RPKM_RA,
  tax_table = tax_table_virus

#using microeco to output heatmap
t5 <- trans_abund$new(dataset = group_CD, taxrank = "genus", input_taxaname =diff_names)#input_taxaname =lefnames
g1 <- t5$plot_heatmap(facet ="enterotype", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
g1
g1 + theme(axis.text.y = element_text(face = 'italic'))


#figure2b 

dataset_bacteria_otu <- microtable$new(
  otu_table = bacteria_otu,
  tax_table = tax_table_bacteria,
  sample_table = sample_data
)

# Plot results
t2$plot_diff_bar(use_number = 1:40, width = 0.8, group_order = c("1", "2"))
t2$plot_diff_abund(use_number = 1:20, group_order = c("1", "2"))
t2$plot_diff_cladogram(
  use_taxa_num = 250,
  use_feature_num = 10,
  clade_label_level = 3,
  group_order = c("1", "2")
)
#for different analysis
lefnames <- rownames(t2$res_diff)[1:40] 
lefnames <- sub(".*\\|s__", "", lefnames)

t1 <- trans_abund$new(dataset = group_YN, taxrank = "Species", input_taxaname =lefnames)#input_taxaname =lefnames
g1 <- t1$plot_heatmap(facet =c("enterotype", "CD_HC"), xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
g1
g1 + theme(axis.text.y = element_text(face = 'italic'))

t3 <- trans_abund$new(dataset = group_GZ, taxrank = "Species", ntaxa = 40)#input_taxaname =lefnames
g1 <- t1$plot_heatmap(facet =c("enterotype", "CD_HC"), xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
g1
g1 + theme(axis.text.y = element_text(face = 'italic'))