###############################################################################################
##
## SMT_analysis.r
## Date: 23.05.2021
## Author: Zohar Meir
##
##
#
#
# Overview:
# This script summarizes main steps in analysis of single-meristem-transcriptome profiles, as presented in Meir, Aviezer, et al. (Nat. Plants, 2021)
#
################################################################################################
######
######		module 0 - Download data, generate SMT objects and compare between technical replicates.
###### 		module 1 - KNN-normalization of batch trends.
###### 		module 2 - in-silico ordering WT meristems.
######		module 3 - ordering of mutants based on WT trajectory.
######		module 4 - Additional analyses.
######
################################################################################################


#load  project functions and common gene names
#source("SMT_nat_plants_2021_functions.r")
source("/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/gitHub_code/SMT_nat_plants_2021_functions.r")
smt_load_common_nms()

#load dependency packages
require("RColorBrewer")
require("pheatmap")
require("zoo")
require("stringr")
require("magicaxis")
require("tidyverse")
require("Matrix")
require("tgstat")
require("svglite")










# Module 0
################################################################################################
######
###### module 0: KNN-normalization of batch trends
######
######	step 0.1 - Download data.
######	step 0.2 - Generate SMT objects.
######  step 0.3 - Compare between technical replicates.
######
################################################################################################

#
######	step 0.1 - Download data.
#
#

# download umi_tables
if(!dir.exists("Data/")) { dir.create("Data/") }
smt_download_umiTables()
			  
# download supplementary tables
if(!dir.exists("Supplementary_Tables/")) { dir.create("Supplementary_Tables/") }
smt_download_supTables()


#
######	step 0.2 - Generate SMT objects.
#
#
# you can either a) use the basic functions below to generate the objects from scratch - where you can have data on individual wells (technical replicates), filter by minimal coverage, control gene names etc.
# Or Alternatively, you can b) skip these loading steps and directly use the processed tables analyzed in the manuscript.
#
#
#
# But first, we can generate two pools of technical replicates (typically 12 per meristem), and see for example how good we quantify expression levels in two, physically separated pools of transcripts :
#
# generate two pools from replicates of each meristem (WT only in this example)
smt_load_technical_replicates(minimal_umi_per_replicate = 1e5,
								dst=FALSE,
								sft=FALSE,
								dstsft=FALSE,
								uf=FALSE,
								wt=TRUE,
								keep_empty_wells = FALSE,
								include_spike_ERCCs = FALSE,
								blacklist_meristems = "sMer_mix_mistake_plate29_empty_2",
								blacklist_genes = substring(mod_apoptosis,1,14))
#								
# we can plot examples of such comparison for some genes (note that estimates become naturally more robust as genes are expressed at higher levels)
#
example_genes = c("Solyc01g067540_YE_ERF_12_PUCHI",
				  "Solyc01g091420_YE_LOB_domain_30_CR_LBD30",
				  "Solyc02g077390_YE_WOX_WOX9_S",
				  "Solyc04g078460_YE_Florigen_asparaginase",
				  "Solyc09g090180_YE_TMF_(TFAM_founder)",
				  "Solyc03g118160_YE_Falsiflora")
#				  
example_genes_nms = c("PUCHI",
					  "LOB30",
					  "S-WOX9",
					  "ASPGB1",
					  "TMF",
					  "FA")
#
for (i in 1:length(example_genes)) {
smt_plot_by_replicates (gene_to_plot = example_genes[i],
						gene_nm = example_genes_nms[i],
                        repA_normalilzed_mat = smer_repA_n,
                        repB_normalilzed_mat = smer_repB_n,
                        svg = FALSE,
                        pp = 72,
                        h = 250,
                        w = 250,
                        log2_exp=FALSE,
                        log2_reg = 1,
                        add_equator = FALSE,
                        meristems_metadata = NULL,
                        spearman_cor = TRUE,
                        dots_size = 0.8,
                        cex_txt = 1.5,
                        base_dir = "./",
                        plot_margins = c(5,5,1,1),
                        show_cor = TRUE) 
}

# After having some trust in our technical replicates, we can aggregate the counted UMIs from all wells derived from the same meristem/sample (controlled for amplification biases) - and generate UMI counts per meristem
smt_load_meristems_season2 (minimal_umi_per_meristem = 1e5,
							dst = TRUE,
							sft = TRUE,
							dstsft = TRUE,
							wt = TRUE,
							uf = TRUE,
							annotation_names = TRUE,
							keep_empty_wells = FALSE,
							blacklist_genes  = substring(mod_apoptosis,1,14),
							blacklist_meristems = c("sMer_253", "sMer_288", "sMer_480", "sMer_156",
													"sMer_333", "sMer_111", "sMer_78", "sMer_116"),
							plates_design_dir = "Data/",
							umi_tab_dir = "Data/",
							annotation_csv_path = "Supplementary_Tables/tomato_gene_annotations.csv")



## b) load processed tables directly
# wt_umi = read.csv("Data/SMT_wt_umiTable.csv",row.names=1)
# wt_n = apply(wt_umi,2,function(x) { 1e5 * (x/sum(x)) } )
# wt_md = read.csv("Data/SMT_wt_metadata.csv",row.names=1)
# 
# sft_umi = read.csv("Data/SMT_sft_umiTable.csv",row.names=1)
# sft_n = apply(sft_umi,2,function(x) { 1e5 * (x/sum(x)) } )
# sft_md = read.csv("Data/SMT_sft_metadata.csv",row.names=1)
# 
# dst_umi = read.csv("Data/SMT_dst_umiTable.csv",row.names=1)
# dst_n = apply(dst_umi,2,function(x) { 1e5 * (x/sum(x)) } )
# dst_md = read.csv("Data/SMT_dst_metadata.csv",row.names=1)
# 
# sftdst_umi = read.csv("Data/SMT_sft_dst_umiTable.csv",row.names=1)
# sftdst_n = apply(sftdst_umi,2,function(x) { 1e5 * (x/sum(x)) } )
# sftdst_md = read.csv("Data/SMT_sft_dst_metadata.csv",row.names=1)
# 
# uf_umi = read.csv("Data/SMT_uf_umiTable.csv",row.names=1)
# uf_n = apply(uf_umi,2,function(x) { 1e5 * (x/sum(x)) } )
# uf_md = read.csv("Data/SMT_uf_metadata.csv",row.names=1)
# 
# smer_umi = cbind(wt_umi, sft_umi, dst_umi, sftdst_umi, uf_umi)
# smer_n = cbind(wt_n, sft_n, dst_n, sftdst_n, uf_n)
# smer_md = rbind(wt_md, sft_md, dst_md, sftdst_md, uf_md)



#Regardless whether you chose a) or b), we should now have 3 SMT objects -  that we will serve as basis for downstream analyses:
## smer_umi: Raw UMI counts per meristem  (ROWS: genes, COLUMNS: meristems).
## smer_n: Normalized UMI counts by sampling depth - to units of UMI-per-100k-UMIs (ROWS: genes, COLUMNS: meristems).
## smer_md: Corresponding metadata collected for each meristem - leaf counts, morphology score, etc (ROWS: meristems, COLUMNS: metadata fields).




################################################################################################
######
###### module 1: KNN-normalization of batch trends
######
######	step 1.1 - Find K nearest neighbours for each meristem based on expression of batchy gene signatures.
######  step 1.2 - normalization of gene expression by subtraction based on KNN model.
######  step 1.3 - plot examples (or, sanity-check) for the effect of such normalization on gene variance.
######
################################################################################################


######	step 1.1 - Find K nearest neighbours for each meristem based on expression of batchy gene signatures.
#
# We will first assign a score to each meristem (stratified by genotype), based on expression of the two main (anti-correlated) batch-dependent gene modules
wt_meristems = colnames(smer_n)[smer_md[colnames(smer_n),"genotype"]=="WT"]
batch_genes = read.csv("Supplementary_Tables/sup_table_lateral_gene_modules.csv")
mod_lateral_dissection_1 = batch_genes[batch_genes$lateral_module=="Non-proliferation","gene_annotation"]
mod_lateral_dissection_2 = batch_genes[batch_genes$lateral_module=="Proliferation","gene_annotation"]
batch_features = data.frame(non_proliferation = colSums(smer_n[mod_lateral_dissection_1, wt_meristems]),
						proliferation = colSums(smer_n[mod_lateral_dissection_2, wt_meristems])) %>%						
						t() %>% t() %>%
						scale(center=TRUE, scale=TRUE)
#						
# Then, for each meristem we can find K meristems that are most similar to it, and that will later serve to normalize its expression.			
# NOTE: you must have "tgstat" R package installed in order to perform this step. see: https://github.com/tanaylab/tgstat
knn_matrix_wt = smt_get_norm_knn_matrix(feats = batch_features,
									  k = 15)
#									  
# We can plot examples of such knn selection:
smt_plot_knn_selection (exp_n = smer_n[,wt_meristems],
						knn_matrix = knn_matrix_wt,
						svg = FALSE,
						pp = 72,
						cex_selected_meristem = 1.4,
						cex_knn = 0.7,
						col_selected_meristem = "blue",
						col_knn = "#ff000080",
						plot_margins = c(5,5,4,4),
						h = 280,
						w = 280)






######  step 1.2 - normalization of gene expression by subtraction based on KNN model.
# We do this by simply subtract mean expression of K neareast to each meristem (according to batch trends) from its own expression
wt_norm_smer_n = smt_batch_normalization(mat_n = smer_n, meristems = wt_meristems, k_meristems=15)
# In theory -  Genes that their expression is co-varying independently of the batch trends, we expect expression variance to be maintained
# (as the mean expression of their K neighbours should be similar to the overall mean)
# However, for genes that are associated with the batch trend - we expect to see a strong reduction in their variance following normalization.


			 
######  step 1.3 - plot examples (or, sanity-check) for the effect of such normalization on gene variance.
example_batch_genes = c("Solyc05g055440_SG_Histone_H2B", "Solyc01g099830_SG_60S_ribosomal")
example_non_batch_genes = c("Solyc01g067540_YE_ERF_12_PUCHI", "Solyc07g040890_SG_Gibberellin_receptor")
example_non_batch_genes_nms = c("PUCHI","CXE20")

# before KNN normalization							
smt_plot_xy (genes_a = example_batch_genes[1],
			 genes_b = example_batch_genes[2],
			 genes_a_label = paste0(substring(example_batch_genes[1],1,14)," (before KNN-norm.)"),
			 genes_b_label = paste0(substring(example_batch_genes[2],1,14)," (before KNN-norm.)"),						 
			 normalized_expression_mat = smer_n,
			 ignore_meristems = setdiff(colnames(smer_n),wt_meristems),
			 scatter_margins = c(5,5,1,1),
			 base_dir = "./",
			 show_legend=FALSE,
			 scatter_height = 420,
			 scatter_width = 420,
			 legend_pos="topright",
			 save_to_file=TRUE,
			 test_meristems=FALSE)
			 
smt_plot_xy (genes_a = example_non_batch_genes[1],
			 genes_b = example_non_batch_genes[2],
			 genes_a_label = paste0(example_non_batch_genes_nms[1]," (before KNN-norm.)"),
			 genes_b_label = paste0(example_non_batch_genes_nms[2]," (before KNN-norm.)"),						 
			 normalized_expression_mat = smer_n,
			 ignore_meristems = setdiff(colnames(smer_n),wt_meristems),
			 scatter_margins = c(5,5,1,1),
			 base_dir = "./",
			 scatter_height = 420,
			 scatter_width = 420,
			 show_legend=FALSE,
			 legend_pos="topright",
			 save_to_file=TRUE,
			 test_meristems=FALSE)			 

# and after KNN normalization							
smt_plot_xy (genes_a = example_batch_genes[1],
			 genes_b = example_batch_genes[2],
			 genes_a_label = paste0(substring(example_batch_genes[1],1,14)," (after KNN-norm.)"),
			 genes_b_label = paste0(substring(example_batch_genes[2],1,14)," (after KNN-norm.)"),						 
			 normalized_expression_mat = wt_norm_smer_n,
			 ignore_meristems = setdiff(colnames(wt_norm_smer_n),wt_meristems),
			 scatter_margins = c(5,5,1,1),
			 base_dir = "./",
			 scatter_height = 420,
			 scatter_width = 420,			 
			 show_legend=FALSE,
			 legend_pos="topright",
			 save_to_file=TRUE,
			 test_meristems=FALSE)
			 
smt_plot_xy (genes_a = example_non_batch_genes[1],
			 genes_b = example_non_batch_genes[2],
			 genes_a_label = paste0(example_non_batch_genes_nms[1]," (after KNN-norm.)"),
			 genes_b_label = paste0(example_non_batch_genes_nms[2]," (after KNN-norm.)"),						 
			 normalized_expression_mat = wt_norm_smer_n,
			 ignore_meristems = setdiff(colnames(wt_norm_smer_n),wt_meristems),
			 scatter_margins = c(5,5,1,1),
			 scatter_height = 420,
			 scatter_width = 420,
			 base_dir = "./",
			 show_legend=FALSE,
			 legend_pos="topright",
			 save_to_file=TRUE,
			 test_meristems=FALSE)	


		
	




	
################################################################################################
######
###### module 2: in-silico ordering WT meristems
######
######	step 2.1 - Generation of similarity matrix between meristems
######			
######  step 2.2 - Ordering of wild-type by slanting
######
################################################################################################


######	step 2.1 - Generation of similarity matrix between meristems
# After selection of feature genes (as described in Meir, Aviezer et al. 2021), we simply compute Spearman's correlation between meristems based only on the KNN-normalized expression of these stage marker genes.
mrkr_genes = read.csv("Supplementary_Tables/sup_table_stage_gene_modules.csv")

# For At fans in the crowd, you can translate some of this to Arabidopsis language using Ensembl One-to-One comparison
mrkr_genes_arabidopsis = smt_findAnnotations_from_arabidopsis(solycIDs = substring(as.character(mrkr_genes$gene_annotation),1,14),
															  source_annotation_csvFile = "Supplementary_Tables/Ensembl_6652_tomato_arabidopsis_oneToOne_orthologs.csv")

# computing meristem-meristem correlations
wt_wt_corr = cor(wt_norm_smer_n[as.character(mrkr_genes$gene_annotation), wt_meristems], method="spearman")
							   

######  step 2.2 - Ordering of wild-type by slanting
# see vignette explaining principles of ordering-samples-by-slanting on this exact example (207 WT meristems collected during season #2) here:
# https://tanaylab.github.io/slanter/articles/meristems.html

# first we can keep only meaningful positive correlations between meristems
wt_wt_corr_pos = wt_wt_corr; 
wt_wt_corr_pos[wt_wt_corr_pos<0.25] = 0; 

# optionally: initialize order by morphology, just because ordering wouldn't be able to define the direction (two options: early-->late or late-->early)
# (Ordering isn't relying on that and can perform well without it in this case at least)
initial_order = case_when(smer_md[wt_meristems,"dev_stage"]=="VM" ~ 1,
						  smer_md[wt_meristems,"dev_stage"]=="TM0" ~ 2,
						  smer_md[wt_meristems,"dev_stage"]=="TM1" ~ 3,
						  smer_md[wt_meristems,"dev_stage"]=="TM2" ~ 4,
						  smer_md[wt_meristems,"dev_stage"]=="LTM" ~ 5,
						  smer_md[wt_meristems,"dev_stage"]=="EFM" ~ 6)

# then we find solution that optimize proximity of this values to the diagonal
wt_wt_reordered = slanter::slanted_reorder(data = wt_wt_corr[order(initial_order),order(initial_order)],
   										   order_data = wt_wt_corr_pos [order(initial_order),order(initial_order)],
   										   order_rows = TRUE,
   										   order_cols = TRUE,
   										   squared_order = TRUE)
wt_oclust_reorder = slanter::oclust(dist(wt_wt_corr[rownames(wt_wt_reordered),colnames(wt_wt_reordered)]))

# plot order
mm = quantile(abs(wt_wt_corr), 0.99, na.rm=TRUE)
wtwtcor_breaks = sort(unique (c( min(wt_wt_corr), seq(-mm, mm, l=99),max(pmin(pmax(wt_wt_corr,(-1)),1)))))
wt_wt_corr_color =   rev(colorRampPalette( (brewer.pal(n = 9, name ="RdGy")))(100))
mers_df = data.frame(stage = smer_md$dev_stage)
					 #batch = smer_md$date_diesct)
rownames(mers_df) = rownames(smer_md)
Var1        <- c("green", "green3","yellow","orange","grey","red","purple","blue")
names(Var1) <- c("VM","TM0","TM1","TM2","TM2E","LTM","EFM","FM")
anno_colors <- list(stage = Var1)
									
smt_mypheatmap(wt_wt_reordered,
				h=1600,w=1600,
				center0=TRUE,outdir="./",
				fname = paste0(myDate,"WT_WT_correlations.png"),
				svg = FALSE,
				cluster_cols=wt_oclust_reorder,
				cluster_rows=wt_oclust_reorder,			
				treeheight_row = 0,
				treeheight_col=0,
				cutree_rows = 5,
				cutree_cols = 5,
				show_rownames=FALSE,
				show_colnames=FALSE,
				cellheight=6/2,
				cellwidth=6/2,
				fontsize=12,
				legend=FALSE,
				border_color=NA,
				my_color=wt_wt_corr_color,
				annotation_legend=FALSE,
				annotation_names_col=FALSE,annotation_names_row=FALSE,
				breaks = wtwtcor_breaks,
				annotation_col=mers_df,annotation_row=mers_df,annotation_colors = anno_colors)		
#
# in the plot above we also partitioned the WT meristems into five groups/phases (veg., tr-I, tr-II, tr-III, FI)
# These are not "real" states but just a simplified view on the CONTINOUS process, that can support for example plotting of specific genes:
phase_colors_wt = c("#00ff0035","#ffff0035","#ffA50035","#ff000035","#8A2BE235")
wt_5phases = cutree(wt_oclust_reorder, k=5)
phase_fractions_wt = table(wt_5phases)/length(wt_5phases)
#
#
example_genes = c("Solyc01g067540_YE_ERF_12_PUCHI",
				  "Solyc01g091420_YE_LOB_domain_30_CR_LBD30",
				  "Solyc02g077390_YE_WOX_WOX9_S",
				  "Solyc04g078460_YE_Florigen_asparaginase",
				  "Solyc09g090180_YE_TMF_(TFAM_founder)",
				  "Solyc03g118160_YE_Falsiflora")
#				  
example_genes_nms = c("PUCHI",
					  "LOB30",
					  "S-WOX9",
					  "ASPGB1",
					  "TMF",
					  "FA")
#
example_ylab_color = c(rep("orange",3),"red","green3","yellow3")
#
for (i in 1:length(example_genes)){
	smt_plot_by_ord_season2(
		mat_n = smer_n[, wt_meristems],
		gene_to_plot = example_genes[i],
		log2_exp = FALSE,
		remove_frame = TRUE,
		mat = smer_umi [,wt_meristems],
		gene_nm = example_genes_nms[i],
		log_reg = 1,
		values_las = 2,
		bg_vec = NULL,
		bg_group_colors = phase_colors_wt,
		bg_group_fracs = phase_fractions_wt,
		base_dir = paste0(base_dir,substring(Sys.time(),1,10),"_WT_expressionKinetics_"),
		k_roll = 15,
		ylab_dist = 3.7,
		xlab_dist = 1.2,
		plot_margins = c(6,6,1,1),
		my_xlab = "",
		ylab_col = example_ylab_color[i],
		trend_col = "black",
		k_cex = 1.2,
		my_height=280,
		my_width=560,
		ordered_meristems = rownames(wt_wt_reordered) )
}







################################################################################################
######
###### module 3 - ordering of mutants based on WT trajectory
######
######	step 3.1 - Generation of WT-mutant chimeric similarity matrix between meristems
######			
######  step 3.2 - Ordering of mutant by the wild-type trajectory
######
######  step 3.3 - Temporal-Alignment: pairing of WT and mutants pseudo-times
######			
######  step 3.4 - Examination of deviation of gene expression in mutants over aligned pseudo-time.
######
################################################################################################

######	step 3.1 - Generation of WT-mutant chimeric similarity matrix between meristems
#
# This example will follow comparison of the WT process to two mutants in main flowering pathwats - sft (florigen) and dst (florigen-independent)
#
# a) We will first normalize the two mutant dataset to control for largest batch effect, just as we did for WT before:
sft_meristems = colnames(smer_n) [smer_md [colnames(smer_n),"genotype"] =="sft"]
dst_meristems = colnames(smer_n) [smer_md [colnames(smer_n),"genotype"] =="dst"]
#
sft_norm_smer_n = smt_batch_normalization(mat_n = smer_n, meristems = sft_meristems, k_meristems=15)
dst_norm_smer_n = smt_batch_normalization(mat_n = smer_n, meristems = dst_meristems, k_meristems=15)
#
# we will then compute meristem-meristem correlations based on expression of feature genes, but this time between WT and mutant meristems (rather then WT to WT)
# After generating this chimeric (WT-mutant) correlation matrix, we will use "slanter" to find the order that maximize proximity of positive values to the diagonal.
#
######  step 3.2 - Ordering of mutant by the wild-type trajectory
# (Note: there is no actual execution of any slanting algoeithm in this case, as the WT order is fixed, and there is simply one mutant order with minimal score to our goal function.)
order_sft = smt_order_mutant_by_wt (wt_norm_ordered_mat = wt_norm_smer_n[,rownames(wt_wt_reordered)], 
									mut_norm_mat = sft_norm_smer_n, 
									corr_thresh = 0.25, 
									mrkr_genes = as.character(mrkr_genes$gene_annotation))
									
order_dst = smt_order_mutant_by_wt (wt_norm_ordered_mat = wt_norm_smer_n[,rownames(wt_wt_reordered)], 
									mut_norm_mat = dst_norm_smer_n, 
									corr_thresh = 0.25, 
									mrkr_genes = as.character(mrkr_genes$gene_annotation))		
								
#
# we can now plot the order, and see how morphology-based scores given to mutant distribute along it (we expect it to be generally in line, but not perfect)
# when ploting the chimeric meristem-meristem correlation matrix, we can also partition mutant meristems as we did for WT (again - not for modeling this continuous process, but for visualization purposes only)
smt_mypheatmap(order_sft$corr_mat,
				h=1600,w=1600,
				center0=TRUE,outdir="./",
				fname = paste0(substring(Sys.time(),1,10),"_WT_sft_correlations.png"),
				svg = FALSE,
				cluster_cols = wt_oclust_reorder,
				cluster_rows = order_sft$mutTree,			
				treeheight_row = 0,
				treeheight_col=0,
				cutree_rows = 5,
				cutree_cols = 5,
				show_rownames=FALSE,
				show_colnames=FALSE,
				cellheight=6/2,
				cellwidth=6/2,
				fontsize=12,
				legend=FALSE,
				border_color=NA,
				my_color=wt_wt_corr_color,
				annotation_legend=FALSE,
				annotation_names_col=FALSE,annotation_names_row=FALSE,
				breaks = wtwtcor_breaks,
				annotation_col=mers_df,annotation_row=mers_df,annotation_colors = anno_colors)

smt_mypheatmap(order_dst$corr_mat,
				h=1600,w=1600,
				center0=TRUE,outdir="./",
				fname = paste0(substring(Sys.time(),1,10),"_WT_dst_correlations.png"),
				svg = FALSE,
				cluster_cols = wt_oclust_reorder,
				cluster_rows = order_dst$mutTree,			
				treeheight_row = 0,
				treeheight_col=0,
				cutree_rows = 5,
				cutree_cols = 5,
				show_rownames=FALSE,
				show_colnames=FALSE,
				cellheight=6/2,
				cellwidth=6/2,
				fontsize=12,
				legend=FALSE,
				border_color=NA,
				my_color=wt_wt_corr_color,
				annotation_legend=FALSE,
				annotation_names_col=FALSE,annotation_names_row=FALSE,
				breaks = wtwtcor_breaks,
				annotation_col=mers_df,annotation_row=mers_df,annotation_colors = anno_colors)
				
######  step 3.3 - Temporal-Alignment: pairing of WT and mutants pseudo-times.
# To reduce gene expression deviations that are stemming from mis-alignment of the processes (e.g, a later stage that wasn't collected in one of the datasets),
# and to distinguish between genes that are globally eliminated or induced in mutant from those that are dis-synchronized in time, or expressed at the "right" timing at lower amplitude - 
# we implelemented a methdology we denote as "temporal alignment", where each mutant meristem is compared to WT meristem that show highest similarity to it ("WT sister") - 
# by using this approach we can compare expression of mutant pseudo-time to an aligned WT pseudo-time, comprised exclusively of matching "WT sister" meristems.

# first, we will assign WT-sister to each mutant meristem
wt_sisters_sft = smt_find_wt_sisters (chimeric_ordered_correlation_matrix = order_sft$corr_mat,
									  n_top_meristems = 3)
									  
wt_sisters_dst = smt_find_wt_sisters (chimeric_ordered_correlation_matrix = order_dst$corr_mat,
									  n_top_meristems = 3)
									  
# This is how this selection looks like (plotting 5 examples)
# Note that we pick the median of top meristems (by default: 3), so we impose some smoothing on the selection process, ultimately this parameter should be 1 (single, top hit) - and in any case not very high.
smt_plot_wt_sister_assigment (chimeric_ordered_correlation_matrix = order_sft$corr_mat,
							  n_top_meristems = 3,
							  assigned_wt_sisters = wt_sisters_sft,
							  meristems_metadata = smer_md,
							  wt_phase_colors = phase_colors_wt,
							  wt_phase_assigment = wt_5phases,
							  mutant_genotype_name = "sft",
							  base_dir = "./",
							  svg=FALSE,
							  example_meristems_to_plot = 5,
							  myDate = substring(Sys.time(),1,10))
							  
smt_plot_wt_sister_assigment (chimeric_ordered_correlation_matrix = order_dst$corr_mat,
							  n_top_meristems = 3,
							  assigned_wt_sisters = wt_sisters_dst,
							  meristems_metadata = smer_md,
							  wt_phase_colors = phase_colors_wt,
							  wt_phase_assigment = wt_5phases,
							  mutant_genotype_name = "dst",
							  base_dir = "./",
							  svg=FALSE,
							  example_meristems_to_plot = 5,
							  myDate = substring(Sys.time(),1,10))
							  							  

# Lets look also if the suggested aligned pseudo-time makes sense (e.g, WT-sister progress and are -/+ covering most of the process evenly)
smt_compare_pseudotime_ranks (chimeric_ordered_correlation_matrix = order_sft$corr_mat,
										  assigned_wt_sisters = wt_sisters_sft,
  										  wt_phase_assigment = wt_5phases,
										  wt_phase_colors = phase_colors_wt,
										  meristems_metadata = smer_md,
										  svg = FALSE,
										  base_dir = "./",
										  mutant_genotype_name = "sft",
										  w = 600,
										  h = 600) 		

smt_compare_pseudotime_ranks (chimeric_ordered_correlation_matrix = order_dst$corr_mat,
										  assigned_wt_sisters = wt_sisters_dst,
  										  wt_phase_assigment = wt_5phases,
										  wt_phase_colors = phase_colors_wt,
										  meristems_metadata = smer_md,
										  svg = FALSE,
										  base_dir = "./",
										  mutant_genotype_name = "dst",
										  w = 600,
										  h = 600) 											  

										  
######  step 3.4 - Examination of deviation of gene expression in mutants over aligned pseudo-time.
#
# To get a global view - we can do differential expression analysis of matching phases 
# (Note this will be more sensitive than crude differential expression analysis between all mutant-WT meristems, and yet less sensitive than one-one temporal matching of single meristems
phases_nms = c("veg","tr-I","tr-II","tr-III","FI")
phases_colors = c("green3","yellow3","orange","red","purple")

for (i in 1:length(phases_nms)) {
wt_meristems_to_compare = rownames(wt_wt_reordered)[cutree(wt_oclust_reorder,5)==i]
sft_meristems_to_compare = rownames(order_sft$corr_mat)[cutree(order_sft$mutTree,5)==i]

smt_compare_meristems (group1 = wt_meristems_to_compare,
					   group2 = sft_meristems_to_compare,
					   mat_umi = smer_umi)
					   
smt_plot_meristems_scatter(tot1 = gene_tot1,
							tot2 = gene_tot2,
							base_dir = "./",
							n_reg = 20, 
							fig_nm = paste0("WT_sft_",phases_nms[i],"_comparison.png"),
							top_genes = 10,
							bad_genes = NULL,
							highlight_genes = NULL,
							highlight_genes_nms = NULL,
							highlight_genes_x = NULL,
							highlight_genes_y = NULL,
							highlight_genes_colors = "cyan",
							dots_color = phases_colors[i],
							highlight_cex_bonus = 0.5,
							fig_h=1600, fig_w=1600,
							segments_calibration_right = 1.05,
							segments_calibration_left = 0.75,							
							lab1 = paste0("WT (",phases_nms[i], ")"),
							lab2 = paste0("sft (",phases_nms[i], ")"),
							main=paste0("WT-sft ",phases_nms[i], " comparison"),
							pt_col1 = "darkred", pt_col2 = "darkblue",
							show_gene_ticks = T)
}


										  
# But ultimately, we can just plot expression kinetics of individual genes in mutants VS. their wt-sister counterparts
# (This time, it makes more sense that background demarcating phases will be defined by the partitioning of mutants and not of WT)
sft_5phases = cutree(order_sft$mutTree, k=5)
phase_fractions_sft = table(sft_5phases)/length(sft_5phases)

dst_5phases = cutree(order_dst$mutTree, k=5)
phase_fractions_dst = table(dst_5phases)/length(dst_5phases)

example_genes = c("Solyc01g067540_YE_ERF_12_PUCHI",
				  "Solyc01g091420_YE_LOB_domain_30_CR_LBD30",
				  "Solyc02g077390_YE_WOX_WOX9_S",
				  "Solyc04g078460_YE_Florigen_asparaginase",
				  "Solyc07g040890_SG_Gibberellin_receptor",
				  "Solyc01g095100_SG_WRKY_transcription",
				  "Solyc02g092240_SG_DUF4228_domain",
				  "Solyc02g062750_SG_sugar_facilitator__Solyc02g062760",
				  "Solyc03g118160_YE_Falsiflora",
				  "Solyc06g069430_YE_FUL/AP1_FUL1",
				  "Solyc11g072600",
				  "Solyc04g008330_SG_Glycosyltransferase",
				  "Solyc10g078700_SG_Squamosa_promoter")
#				  
example_genes_nms = c("PUCHI",
					  "LOB30",
					  "S-WOX9",
					  "ASPGB1",
					  "CXE20",
					  "WRKY22",
					  "DUF4228",
					  "HAN",
					  "FA",
					  "FUL1",
					  "TOE1",
					  "UGT74D1",
					  "SPL15")
					  
for (i in 1:length(example_genes))	{
smt_compare_trends_by_ord_season2 (
							mat1 = smer_n[,rownames(order_sft$corr_mat)],
							mat2 = smer_n[,colnames(order_sft$corr_mat)[wt_sisters_sft]],
							base_dir = paste0("./",substring(Sys.time(),1,10),"_sft_aligned_expression_"),
							gene_to_plot = example_genes[i],
							gene_nm = example_genes_nms[i],
							k_roll = 15,
							log2_exp = FALSE,
							log_reg = 1,
							trend_col1 = "black",
							trend_col2 = "blue",
							my_xlab = "",
							ylab_col = "black",
							ylab_dist = 2.5,
							txt_size = 1.5,
							plot_margins = c(1,5,1,1),
							bg_group_colors = phase_colors_wt,
							bg_group_fracs = phase_fractions_sft,
							remove_frame = TRUE,
							trend_lwd = 10,
							my_ylim = NULL,
							values_las = 2,
							xlab_dist = 3.5,
							my_height = 200,
							my_width = 500,
							svg = FALSE)	
}

for (i in 1:length(example_genes))	{
smt_compare_trends_by_ord_season2 (
							mat1 = smer_n[,rownames(order_dst$corr_mat)],
							mat2 = smer_n[,colnames(order_dst$corr_mat)[wt_sisters_dst]],
							base_dir = paste0("./",substring(Sys.time(),1,10),"_dst_aligned_expression_"),
							gene_to_plot = example_genes[i],
							gene_nm = example_genes_nms[i],
							k_roll = 15,
							log2_exp = FALSE,
							log_reg = 1,
							trend_col1 = "black",
							trend_col2 = "red",
							my_xlab = "",
							ylab_col = "black",
							ylab_dist = 2.5,
							txt_size = 1.5,
							plot_margins = c(1,5,1,1),
							bg_group_colors = phase_colors_wt,
							bg_group_fracs = phase_fractions_dst,
							remove_frame = TRUE,
							trend_lwd = 10,
							my_ylim = NULL,
							values_las = 2,
							xlab_dist = 3.5,
							my_height = 200,
							my_width = 500,
							svg = FALSE)	
}


# There are many more data, analyses and fun with meristems to explore!
# Check out this always-updating webpage, where you can currently navigate within transcriptomes of nearly 3,000 single meristems:
# https://tanaylab.weizmann.ac.il/SMT/

















################################################################################################
######
###### module 4 - Additional analyses
######
######	4.1 - Identification of early-acting transition genes.
######			
######  4.2 - Isolation of genomic modifiers to DST by SNP-chip analysis.
######
######	4.3 - Analysis of age-dependent signatures in shoot apical meristems.
######
################################################################################################




######	4.1 - Identification of early-acting transition genes.
#
# This example will run on wild-type meristems, but could be implemented in principle on any genotype.
#
# to isolate candidate TFs that regulate the transient programs we observed (see fig. 1, Meir, Aviezer et al) - 
# We will first divide the early-most vegetative meristems into two groups, based on their correlations with the ever-vegetative genotype (dst sft - double mutant)
wt_early = rownames(wt_wt_reordered)[ wt_5phases < 3 ]
wt_early1 = wt_early[1:50]
wt_early2 = setdiff(wt_early,wt_early1)
#
#we can limit UMIs obtained from a single meristem - to reduce chances of get a signal from spurious expression in only few meristems
wt_early_squeezed = smt_squeeze_meristems(smer_umi[,wt_early], min_umi=1e5, max_umi=1e6)
#
#
# Let's take only genes that at least 90% of their UMIs on the later group
early_acting_genes_wt = smt_filter_genes_by_two_groups (mat = wt_early_squeezed,
														group1 = wt_early1,
														group2 = wt_early2,
														min_gene_umi = 25,
														min_umi_fraction_group2 = 0.9)
#														
# and order their temporal induction based on UMI cumulative distributions
smt_plot_genes_cummulative  (mat = wt_early_squeezed[early_acting_genes_wt,],
 							 group1 = wt_early1,
 							 group2 = wt_early2,
							 add_common_names = TRUE,
							 base_dir = "./",
							 plot_height = 350,
							 plot_width = 500,
							 file_name = paste0(substring(Sys.time(),1,10),"_WT_early_acting_transition_genes.png"),
							 umi_fraction_threshold = 0.3,
							 font_size = 14)

							 
							 
							 
							 
							 
							 
							 
							 
######  4.2 - Isolation of genomic modifiers to DST by SNP-chip analysis.
#
# This analysis is based on data derived from SNP-chip assay detailed in:
# Sim S, Durstewitz G, Plieske J, Wieseke R, Ganal M, et al. (2012) Development of a large SNP genotyping array and generation of high-density genetic maps in tomato. PLoS ONE 7: e40563.
#
# Counting SNPs per defined genomic bins (compared to M82 sample that served as a reference)
dst_snp_analysis = smt_load_snpChip_data (bin_size_SNP_summation = 5e4,
										  snp_csv_path = "Supplementary_Tables/sup_table_dst_SNPchip.csv",
										  chromosome_length_table_path = "Supplementary_Tables/SL2_4_chromosomeLength.csv",
										  myDate = substring(Sys.time(),1,10))
#
#
# Plotting number of SNPs detected in dst (e137, dst-1 allele)
# Optionally - we can point out at genes that are globally induced when DST gene is impaired, and that reside within a genomic island on chromosome 4 which is associated with dst fertility
#
smt_plot_snpChip_data (snp_data = dst_snp_analysis$snp_counts[,1:3],
					   background_data = dst_snp_analysis$background_counts,
					   markers_data = dst_snp_analysis$markers,
   					   highlight_dst_modifiers_chrom4 = FALSE,
					   dst_modifiers_chrom4_path = "Supplementary_Tables/chrom4_dstGenesAnnotation_SL2_4_assembly.csv",
					   text_cnv=2,
					   sample_name = "e137",
					   factor_img=0.225,
					   fc_factor=0.025,
   					   snp_color = "lightsalmon3",
					   width_factor=2,
					   ppi = 72,
					   svg = FALSE)
#
# Plotting number of SNPs detected in dst "modifiers" only (same genetic background of e137, dst-1 allele - but after crossing out the dst mutation itself)
smt_plot_snpChip_data (snp_data = dst_snp_analysis$snp_counts[,c(1,2,4)],
					   background_data = dst_snp_analysis$background_counts,
					   markers_data = dst_snp_analysis$markers,
   					   highlight_dst_modifiers_chrom4 = FALSE,
					   dst_modifiers_chrom4_path = "Supplementary_Tables/chrom4_dstGenesAnnotation_SL2_4_assembly.csv",
					   text_cnv=2,
					   sample_name = "modifiers",
					   snp_color = "deepskyblue2",
					   factor_img=0.225,
					   fc_factor=0.025,
					   width_factor=2,
					   ppi = 72,
					   svg = FALSE)									





					   
######	4.3 - Analysis of age-dependent signatures in shoot apical meristems.
# general note: our data is strongly biased to span the process of floral transition - therefore, we believe that a more neutral analysis should be carried out
# in order to focus properly on transcriptional signatures which are associated with SAM maturation by age (time / leaves)

# And yet - to fish out genes that are changing as SAM matures we used the small cohort of ever-vegetative meristems we collected (48 overall, more to be analyzed soon)
dstsft_n = smer_n[,smer_md$genotype=="dst sft"]
# Also, we performed same screen on vegetative meristems of WT, dst or sft meristems
wt_veg =  rownames(wt_wt_reordered)[cutree(wt_oclust_reorder,5)==1]
sft_veg = rownames(order_sft$corr_mat)[cutree(order_sft$mutTree,5)==1]
dst_veg = rownames(order_dst$corr_mat)[cutree(order_dst$mutTree,5)==1]
veg_n = smer_n[,c(wt_veg,sft_veg,dst_veg)]

# We looked for genes that show strongest correlation with number of leaves produced by those SAMs:						
age_genes = smt_find_age_genes (veg_mat = veg_n,
								dstsft_mat = smer_n[,smer_md$genotype=="dst sft"],
								corr_method = "pearson",
								corr_thresh = 0.25,
								min_gene_umi = 200)
								
# and plot their expression along all meristems (ordered by transition dynamics, or for ever-vegetative genotype - by # leaves produced)
smt_plot_age_genes (mat = smer_n[,setdiff(colnames(smer_n),rownames(smer_md)[smer_md$genotype=="uf"])],
					age_genes_to_plot = rownames(age_genes),
					meristems_metadata = smer_md,
					scale_expression = TRUE,
					add_common_names = TRUE,
					base_dir = paste0("./",substring(Sys.time(),1,10),"_"),
					file_name = "age_expression_heatmap.png",
					show_legend = FALSE,
					k_clusters_genes = 6,
					font_size=10,
					gaps_width = 2,
					wt_meristems_order = rownames(wt_wt_reordered),
					sft_meristems_order = rownames(order_sft$corr_mat),
					dst_meristems_order = rownames(order_dst$corr_mat),
					heatmap_colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
					stage_genes_csv_path = "Supplementary_Tables/sup_table_stage_gene_modules.csv")
					
					
# We can get the strongest markers and plot them individually
example_age_genes = smt_find_age_genes (veg_mat = veg_n,
								dstsft_mat = smer_n[,smer_md$genotype=="dst sft"],
								corr_method = "pearson",
								corr_thresh = 0.45,
								min_gene_umi = 200)

example_gene_nms = substring(rownames(example_age_genes),1,14)
example_gene_nms = gsub ("Solyc10g078700","SPL15", example_gene_nms)
example_gene_nms = gsub ("Solyc02g077920","CNR-SPL3b", example_gene_nms)

for (i in 1:length(example_gene_nms)) {
print (paste0("Plotting gene ", example_gene_nms [i]));
smt_plot_panel_by_order (genes = rownames(example_age_genes)[i], 
						  label = example_gene_nms [i], 
						  trendline_color = "black",
						  trendline_k = 11,
						  dot_size = 3.5,
						  unified_ylim = TRUE,
						  w = 2000,
						  h = 240,
						  meristems_metadata = smer_md,
						  wt_tree = wt_oclust_reorder,
						  sft_tree = order_sft$mutTree,
						  dst_tree = order_dst$mutTree,
						  wt_ordered_meristems = rownames(wt_wt_reordered),
						  sft_ordered_meristems = rownames(order_sft$corr_mat),
						  dst_ordered_meristems = rownames(order_dst$corr_mat),
						  ylab_distance = 4.4,
						  show_genotype_title = TRUE,
						  panel_margins = c(5,5,4,1),
						  plots_margins = c(6,6,4,1),
						  text_size = 2.2,
						  log2_expression = FALSE,
						  svg=FALSE,pp=72)				
}


#
## END SMT_analysis.r
#########################################################################################
