#
# SMT - list of functions:
#####################
#
# 
# 1. smt_download_umiTables <-
#	 function (wis_url = "http://www.wisdom.weizmann.ac.il/~zoharme/SMT_nat_plants_2021_github/")
#
#
# 2. smt_download_supTables <- 
#	 function (wis_url = "http://www.wisdom.weizmann.ac.il/~zoharme/SMT_nat_plants_2021_github/")
#
#
# 3. smt_load_meristems_season1 <- 
#  	function (minimal_umi_per_meristem = 1e4,
#	 		  nightDay = TRUE,
#	 		  annotation_names = TRUE,
#	 		  keep_empty_wells = FALSE,
#	 		  blacklist_genes = NULL,
#	 		  batch2 = TRUE,
#	 		  plates_design_dir = "Data/",
#	 		  umi_tab_dir = "Data/",
#	 		  annotation_csv_path = "Supplementary_Tables/tomato_gene_annotations.csv")
# 		 
# 
# 4. smt_load_meristems_season2 <-
#	 function (minimal_umi_per_meristem = 1e5,
#			  dst = TRUE,
#			  sft = TRUE,
#			  dstsft = TRUE,
#			  wt = TRUE,
#			  uf = TRUE,
#			  annotation_names = TRUE,
#			  keep_empty_wells = FALSE,
#			  blacklist_genes  = NULL,
#			  blacklist_meristems = NULL,
#			  plates_design_dir = "Data/",
#			  umi_tab_dir = "Data/",
#			  annotation_csv_path = "Supplementary_Tables/tomato_gene_annotations.csv")
# 
#
# 5. smt_findAnnotations_from_solycID <- 
#	 function (solycIDs,
#			   source_annotation_csvFile = "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/desc_tomato_short.csv")	
#
#
# 6. smt_findAnnotations_from_arabidopsis <- 
#	 function (solycIDs,
#			  source_annotation_csvFile = "Supplementary_Tables/Ensembl_6652_tomato_arabidopsis_oneToOne_orthologs.csv")
#
#
# 7. smt_agg_UMIs <- 
#	 function (object_name, 
#			   umiTab_dir = "Data/",
#			   pd_dir = "Data/",
#			   b1_umiTab, b2_umiTab)
#
#
# 8. smt_agg_UMIs_single_batch <- 
#	 function (object_name, 
# 			   umiTab_dir = "Data/",
# 			   pd_dir="Data/",
#			   batch_umiTab)
#
#
# 9. smt_technical_replicates_generator <- 
#	 function (object_name, 
#			   umiTab_dir = "Data/",
#			   pd_dir="Data/",
#			   batch_umiTab,
#			   remove_ERCCs = FALSE,
#			   min_replicate = 1000)
#
# 10. smt_load_technical_replicates <- 
#	function (minimal_umi_per_replicate = 1e3,
#			  dst = TRUE,
#			  sft = TRUE,
#			  dstsft = TRUE,
#			  wt = TRUE,
#			  uf = TRUE,
#			  annotation_names = TRUE,
#			  keep_empty_wells = FALSE,
#			  blacklist_genes  = NULL,
#			  umi_tables_dir = "Data/",
#			  plates_design_dir = "Data/",
#			  blacklist_meristems = NULL,
#			  include_spike_ERCCs = TRUE,								
#			  annotation_csv_path = "Supplementary_Tables/tomato_gene_annotations.csv")
#
# 11. smt_plot_by_replicates <- 
#		function(gene_to_plot = NULL,
#				 gene_nm = NULL,
#				 svg=FALSE,
#				 repA_normalilzed_mat = smer_repA_n,
#				 repB_normalilzed_mat = smer_repB_n,
#				 pp=72,
#				 h = 250,
#				 w = 250,
#				 log2_exp=FALSE,
#				 log2_reg = 1,
#				 add_equator = FALSE,
#				 meristems_metadata = NULL,
#				 spearman_cor = TRUE,
#				 dots_size = 0.8,
#				 cex_txt = 1.5,
#				 myDate = substring(Sys.time(),1,10),
#				 base_dir = "./",
#				 plot_margins = c(5,5,1,1),
#				 show_cor = TRUE)
#
# 12. smt_ds <- 
#		function (umi_tab, 
#				  ds_val)
#
#
# 13. smt_plot_knn_selection <- 
#		function(exp_n = smer_n,
# 	    knn_matrix = NULL,
# 	    svg = FALSE,
# 	    pp = 72,
# 	    cex_selected_meristem = 1.4,
# 	    cex_knn = 0.7,
# 	    col_selected_meristem = "blue",
# 	    col_knn = "#ff000080",
# 	    plot_margins = c(5,5,4,4),
# 	    h = 280,
# 	    w = 280)
# 
# 14. smt_batch_normalization <- 
#		function (mat_n = NULL, 
# 				  meristems = NULL, 
# 				  k_meristems = NULL) 
#
# 15. smt_get_norm_knn_matrix <- 
#		function(feats, 
#				 k)
#
# 16. smt_normalize_knn <- 
#		function(raw, 
#				 knn_mat)
#
# 17. smt_order_mutant_by_wt <- 
#		function(wt_norm_ordered_mat, 
#				 mut_norm_mat, 
#				 corr_thresh, 
#				 mrkr_genes)
#
# 18. smt_plot_xy <- 
#		function(genes_a = NULL,
#				 genes_b = NULL,
#				 genes_a_label = "label A",
#				 genes_b_label = "label B",						 
#				 log2_exp = FALSE,
#				 color_by_dev=TRUE,
#				 svg=FALSE,
#				 pp=72,
#				 reg_log=1,
#				 normalized_expression_mat = smer_n,
#				 ignore_meristems=NULL,
#				 scatter_margins = c(5,5,1,1),
#				 scatter_width=280,
#				 scatter_height=280,
#				 base_dir = "./",
#				 myDate = substr(Sys.time(),1,10),
#				 show_legend=FALSE,
#				 legend_pos="topright",
#				 legend_cex=1.5,
#				 my_xlim=NULL,
#				 my_ylim=NULL,
#				 my_cex=0.8,
#		 		 save_to_file=FALSE,
#				 test_meristems=FALSE)
#
#
# 19. smt_mypheatmap <- 
#		function (mat, 
#				  w=1000, 
#				  h=1000, 
#				  fname, outdir, center0 = F, 
#				  svg = FALSE, svg_ppi=72,
# 				  my_color = colorRampPalette( rev(brewer.pal(n = 7, name ="RdBu")))(100),
# 				  external_breaks_for_center0 = NULL,
# 				  color_quantile = 0.99, ...)
# 
# 20. smt_plot_by_ord_season2 <- 
#		function(mat_n = smer_n,
# 				 mat = smer_umi,
# 				 base_dir,gene_to_plot,gene_nm,
# 				 k_roll=11,
# 				 log2_exp=FALSE,log_reg=1,
# 				 trend_col="black",
# 				 manual_colors = NULL,
# 				 ignore_meristems = NULL,
# 				 my_xlab="ordered meristems",
# 				 ylab_col = "black",
# 				 ylab_dist = 1.8,
# 				 txt_size = 1.8,
# 				 plot_axs = "r",
# 				 plot_margins = c(8,5,1,1),
# 				 max_umi_meristem=NULL,
# 				 bg_vec=NULL,
# 				 bg_group_colors=NULL,
# 				 bg_group_fracs=NULL,
# 				 remove_frame = FALSE,
# 				 k_cex=0.7,
# 				 my_cex=1,
# 				 save_to_file=TRUE,
# 				 break_by_groups = NULL,
# 				 my_ylim=NULL,
# 				 values_las = 2,
# 				 svg=FALSE,
# 				 pp=72,
# 				 plot_grid = FALSE,
# 				 figure_main="",
# 				 ordered_meristems = NULL,
# 				 dev_colorbar=FALSE,xlab_dist=3.5,my_height=300,my_width=500,
# 				 show_colorbar=FALSE,image_mult_cols=0.1)
# 
#
# 21. smt_find_wt_sisters <- 
#		function (chimeric_ordered_correlation_matrix = NULL,
#				  n_top_meristems = 3)
#
# 22. smt_plot_wt_sister_assigment <- 
#		function (chimeric_ordered_correlation_matrix = NULL,
#				  n_top_meristems = 3,
#				  assigned_wt_sisters = NULL,
#				  meristems_metadata = smer_md,
#				  wt_phase_colors = NULL,
#				  w = 600,
#				  h = 200,
#				  wt_phase_assigment = NULL,
#				  mutant_genotype_name = "",
#				  base_dir = "./",
#				  svg=FALSE,
#				  ppi=72,
#				  example_meristems_to_plot = 5,
#				  myDate = substring(Sys.time(),1,10))
#
#
# 23. smt_compare_pseudotime_ranks <- 
#		function (chimeric_ordered_correlation_matrix = NULL,
#				  assigned_wt_sisters = NULL,
#  				  wt_phase_assigment = NULL,
#				  wt_phase_colors = NULL,
#				  meristems_metadata = smer_md,
#				  svg = FALSE,
#				  ppi = 72,
#				  base_dir = "./",
#				  myDate = substring(Sys.time(),1,10),
#				  mutant_genotype_name = "",
#				  w = 600,
#				  h = 600)
#
#
#
# 24. smt_compare_trends_by_ord_season2 <- 
#		function(mat1 = NULL,
#				 mat2 = NULL,
#				 base_dir = "./",
#				 gene_to_plot,
#				 gene_nm,
#				 k_roll = 11,
#				 log2_exp = FALSE,
#				 log_reg = 0.1,
#				 trend_col1 = "black",
#				 trend_col2 = "red",
#				 ignore_meristems = NULL,
#				 my_xlab = "ordered meristems",
#				 ylab_col = "black",
#				 ylab_dist = 1.8,
#				 txt_size = 1.8,
#				 plot_margins = c(8,5,1,1),
#				 bg_group_colors = NULL,
#				 bg_group_fracs = NULL,
#				 remove_frame = FALSE,
#				 trend_lwd = 0.7,
#				 my_ylim = NULL,
#				 values_las = 2,
#				 xlab_dist = 3.5,
#				 my_height = 300,my_width=500,
#				 ppi = 72,
#				 svg = FALSE)
#
#
# 25. smt_squeeze_meristems <- 
#		function (mat = smer_umi,
# 				  min_umi = 1e4,
# 				  max_umi = 2e5)								 
# 
# 26. smt_filter_genes_by_two_groups <- 
#		function (mat = NULL,
# 				  group1 = NULL,
# 				  group2 = NULL,
# 				  min_gene_umi = 25,
# 				  min_umi_fraction_group2 = 0.9) 
#
# 27. smt_plot_genes_cummulative <- 
#		function (mat = NULL,
#  				  group1 = NULL,
#  				  group2 = NULL,
# 				  add_common_names = TRUE,
# 				  base_dir = "./",
# 				  plot_height = 250,
# 				  plot_width = 500,
# 				  show_legend = TRUE,
# 				  file_name = "early_acting_transition_genes.png",
# 				  heatmap_colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)[15:100],
# 				  stage_genes_csv_path = "Supplementary_Tables/sup_table_stage_gene_modules.csv",
#  				  umi_fraction_threshold = 0.3,
# 				  plot_svg = FALSE,
# 				  font_size = 14)
# 
#
# 28. smt_load_snpChip_data <- 
#		function (snp_csv_path = "Supplementary_Tables/sup_table_dst_SNPchip.csv",
# 				  chromosome_length_table_path = "Supplementary_Tables/SL2_4_chromosomeLength.csv",
# 				  bin_size_SNP_summation = 5e4,
# 				  myDate = substring(Sys.time(),1,10))
# 
# 
# 29. smt_plot_snpChip_data <- 
#		function (snp_data = NULL,
# 				  background_data = NULL,
# 				  markers_data = NULL,
# 				  sample_name="",
#    			  highlight_dst_modifiers_chrom4 = FALSE,
# 				  dst_modifiers_chrom4_path = "Supplementary_Tables/chrom4_dstGenesAnnotation_SL2_4_assembly.csv",
# 				  text_cnv=2,
# 				  base_dir = "./",
# 				  snp_color = "blue",
# 				  myDate = substring(Sys.time(),1,10),
# 				  factor_img=0.225,
# 				  fc_factor=0.025,
# 				  width_factor=2,
# 				  ppi = 72,
# 				  svg = FALSE)
#
#
# 30. smt_find_age_genes <- 
#		function (veg_mat = smer_n,
# 				  dstsft_mat = smer_n[,smer_md$genotype=="dst sft"],
# 				  corr_method = "pearson",
# 				  corr_thresh = 0.25,
# 				  meristems_metadata = smer_md,
# 				  min_gene_umi = 200) 
# 
# 
# 
# 
# 31. smt_plot_age_genes <- 
#		function(mat = smer_n,
# 				 age_genes_to_plot = NULL,
# 				 meristems_metadata = smer_md,
# 				 scale_expression = TRUE,
# 				 add_common_names = TRUE,
# 				 base_dir = paste0("./",substring(Sys.time(),1,10),"_"),
# 				 file_name = "age_expression_heatmap.png",
# 				 show_legend = FALSE,
# 				 font_size = 14,
# 				 k_clusters_genes = 2,
# 				 gaps_width = 2,
# 				 wt_meristems_order = NULL,
# 				 sft_meristems_order = NULL,
# 				 dst_meristems_order = NULL,
# 				 dst_sft_order = colnames(dstsft_n)[order(smer_md[colnames(dstsft_n),"Leaf_Num"])],
# 				 plot_height=1200,
# 				 plot_width=2000,
# 				 plot_cellheight = 3.5,
# 				 plot_cellwidth = 2,
# 				 heatmap_colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
# 				 stage_genes_csv_path = "Supplementary_Tables/sup_table_stage_gene_modules.csv")
# 							   
# 							   
# 							   
# 32. smt_findAnnotations_from_arabidopsis <- 
#		function (solycIDs,
# 					  source_annotation_csvFile = "Supplementary_Tables/Ensembl_6652_tomato_arabidopsis_oneToOne_orthologs.csv")								   
# 												   
# 												   
# 												   
# 33. smt_plot_panel_by_order <- 
#		function(genes, 
#				 label, 
# 				 trendline_color=NA,
# 				 trendline_k = 15,
# 				 trendline_cex = 2,
# 				 dot_size=3.5,unified_ylim=FALSE,
# 				 w = 2000,
# 				 h = 240,
# 				 base_dir = "./",
# 				 myDate = substring(Sys.time(),1,10),
# 				 meristems_metadata = smer_md,
# 				 wt_tree = NULL,
# 				 sft_tree = NULL,
# 				 dst_tree = NULL,
# 				 wt_ordered_meristems =  NULL,
# 				 sft_ordered_meristems = NULL,
# 				 dst_ordered_meristems = NULL,
# 				 ylab_distance=4.4,
# 				 show_genotype_title=TRUE,
# 				 panel_margins = c(5,5,4,1),
# 				 plots_margins = c(6,6,4,1),
# 				 text_size = 2.2,
# 				 log2_expression = FALSE,
# 				 svg=FALSE,ppi=72)
# 
# 
# 34. smt_compare_meristems <- 
#		function (group1 = NULL,
#				  group2 = NULL,
#				  mat_umi=smer_umi,
#				  max_umi_per_meristem = 2e6,
#				  min_umi_per_meristem = 1e5,
#				  ds=FALSE)
#
# 
# 35. smt_plot_meristems_scatter <- 
#		function(tot1 = NULL, 
# 				 tot2 = NULL,
# 				 dots_cex = 0.5, 
# 				 text_cex=2,
# 				 base_dir = "./",
# 				 myDate = substring(Sys.time(),1,10),
# 				 n_reg = 5, 
# 				 fig_nm = NULL,
# 				 top_genes = 10,
# 				 bad_genes = NULL,
# 				 highlight_genes = NULL,
# 				 highlight_genes_nms = NULL,
# 				 highlight_genes_x = NULL,
# 				 highlight_genes_y = NULL,
# 				 highlight_genes_colors = "cyan",
# 				 dots_color = "black",
# 				 gene_characters_to_show = 28,
# 				 highlight_cex_bonus = 0.5,
# 				 fig_h=900, fig_w=1200,
# 				 segments_calibration_right = 1.1,
# 				 segments_calibration_left = 1.1,										
# 				 figure_margins = c(25,25,25,25),
# 				 lab1="grp1", lab2="grp2",
# 				 main="compare bulk",
# 				 plot_grid = FALSE,
# 				 pt_col1 = "darkred", pt_col2 = "darkblue",
# 				 show_gene_ticks = T)
# 
# 
# 
# 36. smt_load_common_nms <- 
#		function() 
#
#
#
#
###############################################################################





smt_download_umiTables <- function (wis_url = "http://www.wisdom.weizmann.ac.il/~zoharme/SMT_nat_plants_2021_github/")
{
umiTable_files = c(
"ABZM0269.txt",
"ABZM0270.txt",
"ABZM0271.txt",
"ABZM0272.txt",
"ABZM0273.txt",
"ABZM0274.txt",
"ABZM0275.txt",
"ABZM0276.txt",
"ABZM0277.txt",
"ABZM0278.txt",
"ABZM0279.txt",
"ABZM0280.txt",
"ABZM0303.txt",
"ABZM0304.txt",
"ABZM0305.txt",
"ABZM0306.txt",
"ABZM0307.txt",
"ABZM0308.txt",
"ABZM0309.txt",
"ABZM0310.txt",
"ABZM0311.txt",
"ABZM0312.txt",
"ABZM0313.txt",
"ABZM0314.txt",
"ABZM0315.txt",
"ABZM0316.txt",
"ABZM0370.txt",
"ABZM0371.txt",
"ABZM0372.txt",
"ABZM0373.txt",
"ABZM0374.txt",
"ABZM0375.txt",
"ABZM0376.txt",
"ABZM0377.txt",
"ABZM0378.txt",
"ABZM0403.txt",
"ABZM0404.txt",
"ABZM0405.txt",
"ABZM0406.txt",
"ABZM0407.txt",
"ABZM0408.txt",
"ABZM0409.txt",
"ABZM0410.txt",
"ABZM0411.txt",
"ABZM0413.txt",
"ABZM0414.txt",
"ABZM0415.txt",
"ABZM0416.txt",
"ABZM0417.txt",
"ABZM0418.txt",
"ABZM0419.txt",
"ABZM0420.txt",
"ABZM0421.txt",
"ABZM0422.txt")

metadata_files  = c(
"SMT_dst_metadata.csv",
"SMT_dst_umiTable.csv",
"SMT_sft_dst_metadata.csv",
"SMT_sft_dst_umiTable.csv",
"SMT_sft_metadata.csv",
"SMT_sft_umiTable.csv",
"SMT_uf_metadata.csv",
"SMT_uf_umiTable.csv",
"SMT_wt_lowDepth_metadata.csv",
"SMT_wt_lowDepth_umiTable.csv",
"SMT_wt_metadata.csv",
"SMT_wt_umiTable.csv",
"SMT_metadata_season1.txt",
"SMT_metadata_season2.csv")

platesDesign_files = c(
"plateDesign_singleMers_WT_45to74.txt",
"plateDesign_singleMers_WT_75to106.txt",
"plateDesign_singleMers_WT_107to137.txt",
"plateDesign_singleMers_WT_150to181.txt",
"plateDesign_singleMers_WT_182to209.txt",
"plateDesign_singleMers_WT_210to221_plus11_dst.txt",
"plateDesign_singleMers_plate11_wt_night_1to28.txt",
"plateDesign_singleMers_plate12_wt_day_29to60.txt",
"plateDesign_singleMers_plate13_wt_day_61to92.txt",
"plateDesign_singleMers_plate14_wt_day_night_93to117.txt",
"plateDesign_singleMers_plate15_wt_night_5mads_118to141.txt",
"plateDesign_singleMers_plate16_wt_night_13mads_142to157.txt",
"plateDesign_singleMers_plate17_mads_ful1ful2cr17_9to34.txt",
"plateDesign_singleMers_plate22_sft_dst_dynabeads.txt",
"plateDesign_singleMers_plate23_sft_dst_dynabeads.txt",
"plateDesign_singleMers_plate24_sft_dst_dynabeads.txt",
"plateDesign_singleMers_plate25_sft_dst_dynabeads.txt",
"plateDesign_singleMers_plate26_sft_dst_dynabeads.txt",
"plateDesign_singleMers_plate27_wt_dynabeads.txt",
"plateDesign_singleMers_plate28_wt_dynabeads.txt",
"plateDesign_singleMers_plate29_wt_dynabeads.txt",
"plateDesign_singleMers_plate30_wt_dynabeads.txt",
"plateDesign_singleMers_plate31_dst_dynabeads.txt",
"plateDesign_singleMers_plate32_dst_dynabeads.txt",
"plateDesign_singleMers_plate33_dst_dynabeads.txt",
"plateDesign_singleMers_plate34_sft_dynabeads.txt",
"plateDesign_singleMers_plate35_sft_dynabeads.txt",
"plateDesign_singleMers_plate36_wt_dynabeads.txt",
"plateDesign_singleMers_plate37_wt_dynabeads.txt",
"plateDesign_singleMers_plate38_sft_dynabeads.txt",
"plateDesign_singleMers_plate39_sft_dynabeads.txt",
"plateDesign_singleMers_plate40_dst_dynabeads.txt",
"plateDesign_singleMers_plate41_dst_dynabeads.txt",
"plateDesign_singleMers_plate42_uf_dynabeads.txt",
"plateDesign_singleMers_plate43_uf_dynabeads.txt",
"plateDesign_singleMers_plate44_uf_dynabeads.txt",
"plateDesign_singleMers_plate45_wt_dynabeads.txt",
"plateDesign_singleMers_plate46_wt_dynabeads.txt",
"plateDesign_singleMers_plate47_wt_dynabeads.txt",
"plateDesign_singleMers_plate48_dstsft_dynabeads.txt",
"plateDesign_singleMers_plate49_dstsft_dynabeads.txt",
"wellsCells_singleMers_WT_45to74.txt",
"wellsCells_singleMers_WT_75to106.txt",
"wellsCells_singleMers_WT_107to137.txt",
"wellsCells_singleMers_WT_150to181.txt",
"wellsCells_singleMers_WT_182to209.txt",
"wellsCells_singleMers_WT_210to221_plus11_dst.txt",
"wellsCells_singleMers_plate11_wt_night_1to28.txt",
"wellsCells_singleMers_plate12_wt_day_29to60.txt",
"wellsCells_singleMers_plate13_wt_day_61to92.txt",
"wellsCells_singleMers_plate14_wt_day_night_93to117.txt",
"wellsCells_singleMers_plate15_wt_night_5mads_118to141.txt",
"wellsCells_singleMers_plate16_wt_night_13mads_142to157.txt",
"wellsCells_singleMers_plate17_mads_ful1ful2cr17_9to34.txt",
"wellsCells_singleMers_plate22_sft_dst_dynabeads.txt",
"wellsCells_singleMers_plate23_sft_dst_dynabeads.txt",
"wellsCells_singleMers_plate24_sft_dst_dynabeads.txt",
"wellsCells_singleMers_plate25_sft_dst_dynabeads.txt",
"wellsCells_singleMers_plate26_sft_dst_dynabeads.txt",
"wellsCells_singleMers_plate27_wt_dynabeads.txt",
"wellsCells_singleMers_plate28_wt_dynabeads.txt",
"wellsCells_singleMers_plate29_wt_dynabeads.txt",
"wellsCells_singleMers_plate30_wt_dynabeads.txt",
"wellsCells_singleMers_plate31_dst_dynabeads.txt",
"wellsCells_singleMers_plate32_dst_dynabeads.txt",
"wellsCells_singleMers_plate33_dst_dynabeads.txt",
"wellsCells_singleMers_plate34_sft_dynabeads.txt",
"wellsCells_singleMers_plate35_sft_dynabeads.txt",
"wellsCells_singleMers_plate36_wt_dynabeads.txt",
"wellsCells_singleMers_plate37_wt_dynabeads.txt",
"wellsCells_singleMers_plate38_sft_dynabeads.txt",
"wellsCells_singleMers_plate39_sft_dynabeads.txt",
"wellsCells_singleMers_plate40_dst_dynabeads.txt",
"wellsCells_singleMers_plate41_dst_dynabeads.txt",
"wellsCells_singleMers_plate42_uf_dynabeads.txt",
"wellsCells_singleMers_plate43_uf_dynabeads.txt",
"wellsCells_singleMers_plate44_uf_dynabeads.txt",
"wellsCells_singleMers_plate45_wt_dynabeads.txt",
"wellsCells_singleMers_plate46_wt_dynabeads.txt",
"wellsCells_singleMers_plate47_wt_dynabeads.txt",
"wellsCells_singleMers_plate48_dstsft_dynabeads.txt",
"wellsCells_singleMers_plate49_dstsft_dynabeads.txt")

cat("\nDownloading raw UMI tables..\n\n")
for (i in umiTable_files){
if (file.exists(paste0("Data/",i))) {
cat(paste0("Data/",i," already exists :)\n") ) } else {
download.file(url = paste0(wis_url,"Data/",i), 
			  destfile = paste0("Data/",i))}}
			  
cat("\nDownloading SMT plates design..\n\n")			  	  
for (i in platesDesign_files){
if (file.exists(paste0("Data/",i))) {
cat(paste0("Data/",i," already exists :)\n")) } else {
download.file(url = paste0(wis_url,"Data/",i), 
			  destfile = paste0("Data/",i))}}
			  
				  
cat("\nDownloading processed UMI tables and metadata..\n\n")			  
for (i in metadata_files){
if (file.exists(paste0("Data/",i))) {
cat(paste0("Data/",i," already exists :)\n") ) } else {
download.file(url = paste0(wis_url,"Data/",i), 
			  destfile = paste0("Data/",i))}	}	  
}			  
			  
			  
smt_download_supTables <- function (wis_url = "http://www.wisdom.weizmann.ac.il/~zoharme/SMT_nat_plants_2021_github/")
{	
cat("\nDownloading supplementary tables..\n\n")

sup_tables_files = c(
"sup_table_dst_SNPchip.csv",
"sup_table_lateral_gene_modules.csv",
"sup_table_stage_gene_modules.csv",
"tomato_gene_annotations.csv",
"Ensembl_6652_tomato_arabidopsis_oneToOne_orthologs.csv",
"SL2_4_chromosomeLength.csv",
"chrom4_dstGenesAnnotation_SL2_4_assembly.csv")

for (i in sup_tables_files){
if (file.exists(paste0("Supplementary_Tables/",i))) {
cat(paste0("Supplementary_Tables/",i," already exists :)\n") ) } else {
download.file(url = paste0(wis_url,"Supplementary_Tables/",i), 
			  destfile = paste0("Supplementary_Tables/",i))}}		  		  
}



smt_load_meristems_season1 <- 
function(minimal_umi_per_meristem = 1e4,
		 nightDay = TRUE,
		 annotation_names = TRUE,
		 keep_empty_wells = FALSE,
		 blacklist_genes = NULL,
		 batch2 = TRUE,
		 plates_design_dir = "Data/",
		 umi_tab_dir = "Data/",
		 annotation_csv_path = "Supplementary_Tables/tomato_gene_annotations.csv")
{

print("loading plate3 (WT)..")
mer45_74 = smt_agg_UMIs (object_name="singleMers_WT_45to74",
umiTab_dir = umi_tab_dir,
pd_dir = plates_design_dir,
b1_umiTab ="ABZM0269", b2_umiTab="ABZM0270")

print("loading plate4 (WT)..")
mer75_106 = smt_agg_UMIs (object_name="singleMers_WT_75to106",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0271", b2_umiTab="ABZM0272")

print("loading plate5 (WT)..")
mer107_137 = smt_agg_UMIs (object_name="singleMers_WT_107to137",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0273", b2_umiTab="ABZM0274")

print("loading plate6 (WT)..")
mer150_181 = smt_agg_UMIs (object_name="singleMers_WT_150to181",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0275", b2_umiTab="ABZM0276")

print("loading plate7 (WT)..")
mer182_209 = smt_agg_UMIs (object_name="singleMers_WT_182to209",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0277", b2_umiTab="ABZM0278")

print("loading plate8 (WT)..")
mer210_221 = smt_agg_UMIs (object_name="singleMers_WT_210to221_plus11_dst",
umiTab_dir= umi_tab_dir,
pd_dir=plates_design_dir,
b1_umiTab="ABZM0279", b2_umiTab="ABZM0280")

if (nightDay) {
print("loading plate11 (evening harvest)..")
wt_night_mer1_28 = smt_agg_UMIs (object_name="singleMers_plate11_wt_night_1to28",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0303", b2_umiTab="ABZM0304")

print("loading plate12 (morning harvest)..")
wt_day_mer29_60 = smt_agg_UMIs (object_name="singleMers_plate12_wt_day_29to60",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0305", b2_umiTab="ABZM0306")

print("loading plate13 (morning harvest)..")
wt_day_mer61_92 = smt_agg_UMIs (object_name="singleMers_plate13_wt_day_61to92",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0307", b2_umiTab="ABZM0308")

print("loading plate14 (morning + evening harvest)..")
wt_day_mer93_117 = smt_agg_UMIs (object_name="singleMers_plate14_wt_day_night_93to117",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0309", b2_umiTab="ABZM0310")

print("loading plate15 (evening harvest + ful1/ful2/ap1_65 mutants)..")
mads_night_mer118_141 = smt_agg_UMIs (object_name="singleMers_plate15_wt_night_5mads_118to141",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0311", b2_umiTab="ABZM0312")

print("loading plate16 (evening harvest + ful1/ful2/ap1_65 mutants)..")
mads_night_mer142_157 = smt_agg_UMIs (object_name="singleMers_plate16_wt_night_13mads_142to157",
umiTab_dir= umi_tab_dir,
pd_dir=plates_design_dir,
b1_umiTab="ABZM0313", b2_umiTab="ABZM0314")

print("loading plate17 (morning harvest + ful1/ful2/ap1_65 mutants)..")
mads_mer9_34 = smt_agg_UMIs (object_name="singleMers_plate17_mads_ful1ful2cr17_9to34",
umiTab_dir= umi_tab_dir,
pd_dir= plates_design_dir,
b1_umiTab="ABZM0315", b2_umiTab="ABZM0316")
}



print("merging.. loading metadata per meristem")

if (nightDay){
wt_night_mer1_28_spikes=wt_night_mer1_28[grepl("ERCC-00",rownames(wt_night_mer1_28)),]
wt_night_mer1_28 = wt_night_mer1_28[!grepl("ERCC-00",rownames(wt_night_mer1_28)),]
wt_day_mer29_60_spikes=wt_day_mer29_60[grepl("ERCC-00",rownames(wt_day_mer29_60)),]
wt_day_mer29_60 = wt_day_mer29_60[!grepl("ERCC-00",rownames(wt_day_mer29_60)),]
wt_day_mer61_92_spikes=wt_day_mer61_92[grepl("ERCC-00",rownames(wt_day_mer61_92)),]
wt_day_mer61_92 = wt_day_mer61_92[!grepl("ERCC-00",rownames(wt_day_mer61_92)),]

wt_day_mer93_117_spikes=wt_day_mer93_117[grepl("ERCC-00",rownames(wt_day_mer93_117)),]
wt_day_mer93_117 = wt_day_mer93_117[!grepl("ERCC-00",rownames(wt_day_mer93_117)),]

mads_night_mer118_141_spikes=mads_night_mer118_141[grepl("ERCC-00",rownames(mads_night_mer118_141)),]
mads_night_mer118_141 = mads_night_mer118_141[!grepl("ERCC-00",rownames(mads_night_mer118_141)),]

mads_night_mer142_157_spikes=mads_night_mer142_157[grepl("ERCC-00",rownames(mads_night_mer142_157)),]
mads_night_mer142_157 = mads_night_mer142_157[!grepl("ERCC-00",rownames(mads_night_mer142_157)),]

mads_mer9_34_spikes=mads_mer9_34[grepl("ERCC-00",rownames(mads_mer9_34)),]
mads_mer9_34 = mads_mer9_34[!grepl("ERCC-00",rownames(mads_mer9_34)),]
}



mer45_74_spikes=mer45_74[grepl("ERCC-00",rownames(mer45_74)),]
mer45_74 = mer45_74[!grepl("ERCC-00",rownames(mer45_74)),]
mer75_106_spikes=mer75_106[grepl("ERCC-00",rownames(mer75_106)),]
mer75_106 = mer75_106[!grepl("ERCC-00",rownames(mer75_106)),]
mer107_137_spikes=mer107_137[grepl("ERCC-00",rownames(mer107_137)),]
mer107_137 = mer107_137[!grepl("ERCC-00",rownames(mer107_137)),]
mer150_181_spikes=mer150_181[grepl("ERCC-00",rownames(mer150_181)),]
mer150_181 = mer150_181[!grepl("ERCC-00",rownames(mer150_181)),]
mer182_209_spikes=mer182_209[grepl("ERCC-00",rownames(mer182_209)),]
mer182_209 = mer182_209[!grepl("ERCC-00",rownames(mer182_209)),]
mer210_221_spikes=mer210_221[grepl("ERCC-00",rownames(mer210_221)),]
mer210_221 = mer210_221[!grepl("ERCC-00",rownames(mer210_221)),]

### combine tables
if (nightDay)
{
smer_umi = cbind(mer45_74,mer75_106,mer107_137,mer150_181,mer182_209,mer210_221,wt_night_mer1_28,wt_day_mer29_60,wt_day_mer61_92,wt_day_mer93_117,mads_night_mer118_141,mads_night_mer142_157,mads_mer9_34) 
} else { 
smer_umi = cbind(mer45_74,mer75_106,mer107_137,mer150_181,mer182_209,mer210_221) }



### filter empty wells and bad genes
if (!keep_empty_wells){print("filtering empty wells.."); smer_umi = smer_umi[,!grepl("mpty|MPTY", colnames(smer_umi))]}
if (!is.null(blacklist_genes)){print("filtering blacklisted genes..");smer_umi = smer_umi[rownames(smer_umi)[!rownames(smer_umi) %in% blacklist_genes],]}

print(paste0("removing ",sum(colSums(smer_umi)<minimal_umi_per_meristem), " meristems.."));
smer_umi = smer_umi[,colSums(smer_umi)>=minimal_umi_per_meristem]





if (annotation_names) {
print("mishmashing gene names from different sources...")
rownames(smer_umi) = smt_findAnnotations_from_solycID(solycIDs = rownames(smer_umi),
											          source_annotation_csvFile = annotation_csv_path)
}

print("loading metadata..")
md_table=read.table("Data/SMT_metadata_season1.txt",sep="\t",header=TRUE)
md_table_filt  = md_table[match(colnames(smer_umi),md_table[,1]),]


###############################################
## add metadata per meristem
print ("labeling according to metadata..")
#
# developmental stage colors
md_table_filt$color=NA
md_table_filt$color[md_table_filt$dev_stage=="TM"]="grey"
md_table_filt$color[md_table_filt$dev_stage=="VM"]="green"
md_table_filt$color[md_table_filt$dev_stage=="TM0"]="darkgreen"
md_table_filt$color[md_table_filt$dev_stage=="TM1"]="yellow"
md_table_filt$color[md_table_filt$dev_stage=="TM2"]="darkgoldenrod1"
md_table_filt$color[md_table_filt$dev_stage=="FM"]="darkorange2"
md_table_filt$color[md_table_filt$dev_stage=="LTM"]="darkorange2"
md_table_filt$color[md_table_filt$dev_stage=="EFM"]="red"


md_table_filt$total_umi = colSums(smer_umi)
#
# dissection batch colors
md_table_filt$color_batch=NA
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 0"]="yellow"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 1"]="gold"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 2"]="orange"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 3"]="green"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 4"]="red"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 5"]="purple"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 6"]="grey"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 7"]="black"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 8"]="lightblue"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 9"]="blue"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 10"]="pink"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 11"]="navajowhite3"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 12"]="cyan"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 13"]="darkgrey"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 14"]="brown"
md_table_filt$color_batch[md_table_filt$Plant_Bulk=="bulk 15"]="lightsalmon1"


# time of dissection colors
md_table_filt$time_of_dissection = "7AM"
md_table_filt$time_of_dissection[grepl("_pm_",md_table_filt$eppendorf_Num)] = "7PM"
md_table_filt$color_time = "lightsalmon2"
md_table_filt$color_time[md_table_filt$time_of_dissection == "7PM"] = "deepskyblue3"
#

# total UMI per meristem
md_table_filt$total_umi = colSums(smer_umi)


# genotype
md_table_filt$genotype = "WT"
md_table_filt$genotype[grepl("dst",md_table_filt$eppendorf_Num)] = "dst"
md_table_filt$genotype[grepl("tcs",md_table_filt$eppendorf_Num)] = "tcsVENUS"
md_table_filt$genotype[grepl("mads",md_table_filt$eppendorf_Num)] = "mads"
md_table_filt$genotype[grepl("sft",md_table_filt$eppendorf_Num)] = "sft"

## add general batch
md_table_filt$batch_general = 2
md_table_filt$batch_general[md_table_filt$MARS_plate >=3 & md_table_filt$MARS_plate<=8] = 1

#save metadata into object
smer_md = md_table_filt;

print ("filtering genotypes/batches..")
if (!batch2){
 smer_umi = smer_umi[,!smer_md$batch_general==2]
 smer_md = smer_md[!smer_md$batch_general==2,]
}



print ("normalizing..")
smer_n = apply(smer_umi,2,function(x) {(x/sum(x))*1e5 })

print ("Downsampling..")
smer_downsampled = smt_ds(smer_umi,minimal_umi_per_meristem)

if (!keep_empty_wells) {rownames(smer_md) = smer_md$eppendorf_Num}
	
smer_umi <<- smer_umi;
smer_n <<- smer_n;
smer_downsampled <<- smer_downsampled;
smer_md <<- smer_md;
base_dir <<- "./"
myDate <<- substr(Sys.time(),1,10); 


cat(paste0("\nloaded objects:","\n",
			  "  smer_umi: raw expression table","\n",
			  "  smer_n: expression table, normalized to 100K UMIs","\n",
			  "  smer_downsampled: expression table, down-sampled to minimal threshold","\n",
			  "  smer_md: metadata per meristem","\n\n",
			  ncol(smer_umi)," total meristems above ",
			  minimal_umi_per_meristem," UMIs\n\n"))
}

		 

smt_load_meristems_season2 <-
function (minimal_umi_per_meristem = 1e5,
		  dst = TRUE,
		  sft = TRUE,
		  dstsft = TRUE,
		  wt = TRUE,
		  uf = TRUE,
		  annotation_names = TRUE,
		  keep_empty_wells = FALSE,
		  blacklist_genes  = NULL,
		  blacklist_meristems = NULL,
		  plates_design_dir = "Data/",
		  umi_tab_dir = "Data/",
		  annotation_csv_path = "Supplementary_Tables/tomato_gene_annotations.csv")
{


if (wt) {
print("loading plate27 (wild-type)..")
p27 = smt_agg_UMIs_single_batch(
object_name ="singleMers_plate27_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir = "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab = "ABZM0375")

print("loading plate28 (wild-type)..")
p28 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate28_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0376")

print("loading plate29 (wild-type)..")
p29 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate29_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0377")

print("loading plate30 (wild-type)..")
p30 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate30_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0378")

print("loading plate36 (wild-type)..")
p36 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate36_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0408")

print("loading plate37 (wild-type)..")
p37 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate37_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0409")

print("loading plate45 (wild-type)..")
p45 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate45_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0418")

print("loading plate46 (wild-type)..")
p46 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate46_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0419")

print("loading plate47 (wild-type)..")
p47 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate47_wt_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0420")
}


if (sft | dst){
print("loading plate22 (sft n7187 mutants + dst mutants)..")
p22 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate22_sft_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0370")

print("loading plate23 (sft n7187 mutants + dst mutants)..")
p23 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate23_sft_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0371")

print("loading plate24 (sft n7187 mutants + dst mutants)..")
p24 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate24_sft_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0372")

print("loading plate25 (sft n7187 mutants + dst mutants)..")
p25 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate25_sft_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0373")

print("loading plate26 (sft n7187 mutants + dst mutants)..")
p26 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate26_sft_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0374")
}

if (dst) {
print("loading plate31 (dst mutants)..")
p31 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate31_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0403")

print("loading plate32 (dst mutants)..")
p32 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate32_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0404")

print("loading plate33 (dst mutants)..")
p33 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate33_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0405")

print("loading plate40 (dst mutants)..")
p40 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate40_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0413")

print("loading plate41 (dst mutants)..")
p41 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate41_dst_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0414")

}


if (sft) {
print("loading plate34 (sft n7187 mutants)..")
p34 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate34_sft_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0406")

print("loading plate35 (sft n7187 mutants)..")
p35 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate35_sft_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0407")

print("loading plate38 (sft n7187 mutants)..")
p38 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate38_sft_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0410")

print("loading plate39 (sft n7187 mutants)..")
p39 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate39_sft_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0411")

}

if (uf) {
print("loading plate42 (uf mutants)..")
p42 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate42_uf_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0415")

print("loading plate43 (uf mutants)..")
p43 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate43_uf_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0416")

print("loading plate44 (uf mutants)..")
p44 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate44_uf_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0417")
}

if (dstsft) {
print("loading plate48 (dst;sft double mutants)..")
p48 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate48_dstsft_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0421")

print("loading plate49 (dst;sft double mutants)..")
p49 = smt_agg_UMIs_single_batch(
object_name="singleMers_plate49_dstsft_dynabeads",
umiTab_dir = umi_tab_dir,
pd_dir= "/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/plates_design/",
batch_umiTab="ABZM0422")
}


print("merging.. ")
print("loading metadata per meristem..")

if (wt){
plate27_spikes=p27[grepl("ERCC-00",rownames(p27)),]
plate27 = p27[!grepl("ERCC-00",rownames(p27)),]

plate28_spikes=p28[grepl("ERCC-00",rownames(p28)),]
plate28 = p28[!grepl("ERCC-00",rownames(p28)),]


plate29_spikes=p29[grepl("ERCC-00",rownames(p29)),]
plate29 = p29[!grepl("ERCC-00",rownames(p29)),]

plate30_spikes=p30[grepl("ERCC-00",rownames(p30)),]
plate30 = p30[!grepl("ERCC-00",rownames(p30)),]

plate36_spikes=p36[grepl("ERCC-00",rownames(p36)),]
plate36 = p36[!grepl("ERCC-00",rownames(p36)),]

plate37_spikes=p37[grepl("ERCC-00",rownames(p37)),]
plate37 = p37[!grepl("ERCC-00",rownames(p37)),]

plate45_spikes=p45[grepl("ERCC-00",rownames(p45)),]
plate45 = p45[!grepl("ERCC-00",rownames(p45)),]

plate46_spikes=p46[grepl("ERCC-00",rownames(p46)),]
plate46 = p46[!grepl("ERCC-00",rownames(p46)),]

plate47_spikes=p47[grepl("ERCC-00",rownames(p47)),]
plate47 = p47[!grepl("ERCC-00",rownames(p47)),]
}


if (sft | dst){
plate22_spikes=p22[grepl("ERCC-00",rownames(p22)),]
plate22 = p22[!grepl("ERCC-00",rownames(p22)),]

plate23_spikes=p23[grepl("ERCC-00",rownames(p23)),]
plate23 = p23[!grepl("ERCC-00",rownames(p23)),]

plate24_spikes=p24[grepl("ERCC-00",rownames(p24)),]
plate24 = p24[!grepl("ERCC-00",rownames(p24)),]

plate25_spikes=p25[grepl("ERCC-00",rownames(p25)),]
plate25 = p25[!grepl("ERCC-00",rownames(p25)),]

plate26_spikes=p26[grepl("ERCC-00",rownames(p26)),]
plate26 = p26[!grepl("ERCC-00",rownames(p26)),]
}

if (dst){
plate31_spikes=p31[grepl("ERCC-00",rownames(p31)),]
plate31 = p31[!grepl("ERCC-00",rownames(p31)),]

plate32_spikes=p32[grepl("ERCC-00",rownames(p32)),]
plate32 = p32[!grepl("ERCC-00",rownames(p32)),]

plate33_spikes=p33[grepl("ERCC-00",rownames(p33)),]
plate33 = p33[!grepl("ERCC-00",rownames(p33)),]

plate40_spikes=p40[grepl("ERCC-00",rownames(p40)),]
plate40 = p40[!grepl("ERCC-00",rownames(p40)),]

plate41_spikes=p41[grepl("ERCC-00",rownames(p41)),]
plate41 = p41[!grepl("ERCC-00",rownames(p41)),]
}

if (sft){
plate34_spikes=p34[grepl("ERCC-00",rownames(p34)),]
plate34 = p34[!grepl("ERCC-00",rownames(p34)),]

plate35_spikes=p35[grepl("ERCC-00",rownames(p35)),]
plate35 = p35[!grepl("ERCC-00",rownames(p35)),]

plate38_spikes=p38[grepl("ERCC-00",rownames(p38)),]
plate38 = p38[!grepl("ERCC-00",rownames(p38)),]

plate39_spikes=p39[grepl("ERCC-00",rownames(p39)),]
plate39 = p39[!grepl("ERCC-00",rownames(p39)),]
}

if (uf){
plate42_spikes=p42[grepl("ERCC-00",rownames(p42)),]
plate42 = p42[!grepl("ERCC-00",rownames(p42)),]

plate43_spikes=p43[grepl("ERCC-00",rownames(p43)),]
plate43 = p43[!grepl("ERCC-00",rownames(p43)),]

plate44_spikes=p44[grepl("ERCC-00",rownames(p44)),]
plate44 = p44[!grepl("ERCC-00",rownames(p44)),]
}

if (dstsft){
plate48_spikes=p48[grepl("ERCC-00",rownames(p48)),]
plate48 = p48[!grepl("ERCC-00",rownames(p48)),]

plate49_spikes=p49[grepl("ERCC-00",rownames(p49)),]
plate49 = p49[!grepl("ERCC-00",rownames(p49)),]
}




smer_umi=c();

if (wt) {
smer_umi=cbind(smer_umi,plate27,plate28,plate29,plate30,plate36,plate37,plate45,plate46,plate47)
}

if (sft | dst) {
smer_umi=cbind(smer_umi,plate22,plate23,plate24,plate25,plate26)
}

if (dst) {
smer_umi=cbind(smer_umi,plate31,plate32,plate33,plate40,plate41)
}

if (sft) {
smer_umi=cbind(smer_umi,plate34,plate35,plate38,plate39)
}

if (uf) {
smer_umi=cbind(smer_umi,plate42,plate43,plate44)
}

if (dstsft) {
smer_umi=cbind(smer_umi,plate48,plate49)
}

###################


### filter empty wells and bad genes
if (!keep_empty_wells){print("filtering empty wells.."); smer_umi = smer_umi[,!grepl("mpty|MPTY", colnames(smer_umi))]}
if (!is.null(blacklist_genes)){print("filtering blacklisted genes..");smer_umi = smer_umi[rownames(smer_umi)[!rownames(smer_umi) %in% blacklist_genes],]}

print(paste0("removing ",sum(colSums(smer_umi)<minimal_umi_per_meristem), " small meristems.."));
smer_umi = smer_umi[,colSums(smer_umi)>=minimal_umi_per_meristem]

if (annotation_names) {
print("mishmashing gene names from different sources...")
gene_names_extended = smt_findAnnotations_from_solycID (solycIDs = rownames(smer_umi),
											   source_annotation_csvFile = annotation_csv_path)

rownames(smer_umi)=gene_names_extended;
}

print("loading .csv combined meristems metadata table..")
md_table = read.csv("Data/SMT_metadata_season2.csv",header=TRUE)
md_table_filt = md_table[match(colnames(smer_umi),md_table[,1]),]

###############################################
## add metadata per meristem
print ("labeling according to metadata..")
#
# developmental stage colors
md_table_filt$color=NA
md_table_filt$color[md_table_filt$dev_stage=="TM"]="grey"
md_table_filt$color[md_table_filt$dev_stage=="VM"]="green"
md_table_filt$color[md_table_filt$dev_stage=="TM0"]="green3"
md_table_filt$color[md_table_filt$dev_stage=="TM1"]="yellow"
md_table_filt$color[md_table_filt$dev_stage=="TM2"]="orange"
md_table_filt$color[md_table_filt$dev_stage=="TM2E"]="orange"
md_table_filt$color[md_table_filt$dev_stage=="FM"]="red"
md_table_filt$color[md_table_filt$dev_stage=="LTM"]="red"
md_table_filt$color[md_table_filt$dev_stage=="EFM"]="purple"


#new stages for season3
md_table_filt$color[md_table_filt$dev_stage=="FM1" | md_table_filt$dev_stage=="FM2" | md_table_filt$dev_stage=="FM3" | md_table_filt$dev_stage=="FM4"] = "blue"


# time of dissection colors
md_table_filt$time_of_dissection = "7AM"
md_table_filt$time_of_dissection[grepl("_pm_",md_table_filt$eppendorf_Num)] = "7PM" 	
md_table_filt$color_time = "lightsalmon2"
md_table_filt$color_time[md_table_filt$time_of_dissection == "7PM"] = "deepskyblue3"
#

# total UMI per meristem
md_table_filt$total_umi = colSums(smer_umi)
#

#
#save metada into object
smer_md = md_table_filt;

print ("filtering genotypes..")

if (!dst)
 {
 smer_umi = smer_umi[,!smer_md$genotype=="dst"]
 smer_md = smer_md[!smer_md$genotype=="dst",]
 }
 
if (!sft)
 {
 smer_umi = smer_umi[,!smer_md$genotype=="sft"]
 smer_md = smer_md[!smer_md$genotype=="sft",]
 }

print ("normalizing..")
smer_n = apply(smer_umi,2,function(x) {(x/sum(x))*1e5 })

print ("Downsampling..")
smer_downsampled = smt_ds(smer_umi,minimal_umi_per_meristem)

if (!keep_empty_wells) {rownames(smer_md) = smer_md$eppendorf_Num}

## storing data in permanent objects	
if (!is.null(blacklist_meristems)) {
	print(paste0("removing ",sum(rownames(smer_md) %in% blacklist_meristems), " blacklisted meristems.."));
	nonbl_meristems = setdiff(rownames(smer_md),blacklist_meristems);
	smer_umi <<- smer_umi[,nonbl_meristems];
	smer_n <<- smer_n[,nonbl_meristems];
	smer_downsampled <<- smer_downsampled[,nonbl_meristems];
	smer_md <<- smer_md[nonbl_meristems,];
		} else {
	smer_umi <<- smer_umi;
	smer_n <<- smer_n;
	smer_downsampled <<- smer_downsampled;
	smer_md <<- smer_md; }
	
base_dir <<- "./"
myDate <<- substr(Sys.time(),1,10); 

cat(paste0("\nloaded objects:","\n",
			  "  smer_umi: raw expression table","\n",
			  "  smer_n: expression table, normalized to 100K UMIs","\n",
			  "  smer_downsampled: expression table, down-sampled to minimal threshold","\n",
			  "  smer_md: metadata per meristem","\n\n",
			  ncol(smer_umi)," total meristems above ",
			  minimal_umi_per_meristem," UMIs\n\n"))
}
##############################################################################################
##############################################################################################		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
### this function ->
### gets: a vector of solycs (possibly, fusion of more than one Solyc, separated by ";"
### returns: vector of same length, with additional annotations from YE lab and SOL databases
smt_findAnnotations_from_solycID <- function (solycIDs,
											   source_annotation_csvFile = "Supplementary_Tables/tomato_gene_annotations.csv")								   
											   {
	anots = read.csv(source_annotation_csvFile)
    names_per_gene = sapply(solycIDs,function(x){str_count(x,";")+1})
	new_names=rep("aaa",length(names_per_gene))
	for (i in 1:length(names_per_gene))
	{
		if (names_per_gene[i]==1)
		{
			if (!solycIDs[i] %in% as.character(anots[,2])) 
			{new_names[i]=solycIDs[i]
			} else{
			new_names[i] = as.character(anots[ as.character(anots[,2]) %in% solycIDs[i] ,5])
			}
		} else if (names_per_gene[i]>1) {
			temp_name="";
			for (j in 1:names_per_gene[i]){
			temp_gene = unlist(strsplit(solycIDs[i],";"))[j]
				if (!temp_gene %in% as.character(anots[,2])) {temp_gene_name=temp_gene
				} else {
				temp_gene_name = as.character(anots[ as.character(anots[,2]) %in% temp_gene ,5])
				}
			if (j==1){
				temp_name=temp_gene_name
				}else{
				temp_name=paste0(temp_name,"__",temp_gene_name)
				}
			new_names[i] = temp_name
			}
		}
	}
return(new_names)
}
##############################################################################################
##############################################################################################



### this function ->
### gets: a vector of solycs (possibly, fusion of more than one Solyc, separated by ";"
### returns: vector of same length, with additional annotations from arabidopsis 6,652 one-to-one orthologues (derived from the ENSEMBL database)
smt_findAnnotations_from_arabidopsis <- function (solycIDs,
												   source_annotation_csvFile = "Supplementary_Tables/Ensembl_6652_tomato_arabidopsis_oneToOne_orthologs.csv")								   
												   {
	anots = read.csv(source_annotation_csvFile)
	rownames(anots) = anots[,"Sl_ID"]
    names_per_gene = sapply(solycIDs,function(x){str_count(x,";")+1})
	new_names=rep("aaa",length(names_per_gene))
	for (i in 1:length(names_per_gene))
	{
		if (names_per_gene[i]==1)
		{
			if (!solycIDs[i] %in% as.character(anots[,"Sl_ID"])) 
			{new_names[i]=solycIDs[i]
			} else{
			new_names[i] = paste0(as.character(anots[ as.character(anots[,"Sl_ID"]) %in% solycIDs[i] ,"At_ID"]),"_", as.character(anots[ as.character(anots[,"Sl_ID"]) %in% solycIDs[i] ,"At_name"]))
			}
		} else if (names_per_gene[i]>1) {
			temp_name="";
			for (j in 1:names_per_gene[i]){
			temp_gene = unlist(strsplit(solycIDs[i],";"))[j]
				if (!temp_gene %in% as.character(anots[,"Sl_ID"])) {temp_gene_name=temp_gene
				} else {
				temp_gene_name = paste0(as.character(anots[ as.character(anots[,"Sl_ID"]) %in% temp_gene ,"At_ID"]),"_", as.character(anots[ as.character(anots[,"Sl_ID"]) %in% temp_gene ,"At_name"]))
				}
			if (j==1){
				temp_name=temp_gene_name
				}else{
				temp_name=paste0(temp_name,"__",temp_gene_name)
				}
			new_names[i] = temp_name
			}
		}
	}
return(new_names)
}
##############################################################################################
##############################################################################################


################################################################################################
# This function sums up UMI counts for custom 384 MARS-Seq plate.
# In addition to umi.tab, the function needs 3 tables
#	• relevant wells_cells.txt file (contains relevant plate)
#	• conversion of well_coordinate ("A1,O22..") to object name("cloneDKO246","clone115","megaclone17"...)
#	• conversion of well_names in umi.tab ("PBZM0010382","PBZM0021007"..) to well_coordinate ("A1,O22..")
smt_agg_UMIs <- function (object_name, 
						   umiTab_dir = "Data/",
						   pd_dir = "Data/",
						   b1_umiTab, b2_umiTab)
{
wc_obj = read.table(paste0(pd_dir,"wellsCells_",object_name,".txt"),
				  sep="\t",header=T,colClasses=c(rep("character",2),rep("NULL",11)))		  
pd_obj = read.table(paste0(pd_dir,"plateDesign_",object_name,".txt"),sep="\t",header=TRUE)
b1_obj = read.table(paste0(umiTab_dir,b1_umiTab,".txt"),sep="\t")
b2_obj = read.table(paste0(umiTab_dir,b2_umiTab,".txt"),sep="\t")
umis_obj = cbind(b1_obj,b2_obj)
wc_pd = data.frame(well=wc_obj$Well_ID,
				 clone=pd_obj$ID[match(wc_obj$well_coordinates,pd_obj$well)])
order_clones = match(colnames(umis_obj),wc_pd$well)
return(t(apply(umis_obj,1,function(x){tapply(x,wc_pd$clone[order_clones],sum)})))
}
################################################################################################
################################################################################################



################################################################################################
# This function sums up UMI counts for custom 384 MARS-Seq plate.
# In addition to umi.tab, the function needs 3 tables
#	• relevant wells_cells.txt file (contains relevant plate)
#	• conversion of well_coordinate ("A1,O22..") to object name("cloneDKO246","clone115","megaclone17"...)
#	• conversion of well_names in umi.tab ("PBZM0010382","PBZM0021007"..) to well_coordinate ("A1,O22..")
smt_agg_UMIs_single_batch <- function (object_name, 
									   umiTab_dir = "Data/",
									   pd_dir="Data/",
									   batch_umiTab)
{
wc_obj = read.table(paste0(pd_dir,"wellsCells_",object_name,".txt"),
				  sep="\t",header=T,colClasses=c(rep("character",2),rep("NULL",11)))		  
pd_obj = read.table(paste0(pd_dir,"plateDesign_",object_name,".txt"),sep="\t",header=TRUE)
batch_obj = read.table(paste0(umiTab_dir,batch_umiTab,".txt"),sep="\t")
umis_obj = batch_obj
wc_pd = data.frame(well=wc_obj$Well_ID,
				 clone=pd_obj$ID[match(wc_obj$well_coordinates,pd_obj$well)])
order_clones = match(colnames(umis_obj),wc_pd$well)
return(t(apply(umis_obj,1,function(x){tapply(x,wc_pd$clone[order_clones],sum)})))
}
################################################################################################
################################################################################################
		  

		  
		  
		  
		  
		  
################################################################################################
# This function split all UMIs into two tehcnical replicates (summing up half of the wells of each meristem)
# In addition to umi.tab, the function needs 3 tables
#	• relevant wells_cells.txt file (contains relevant plate)
#	• conversion of well_coordinate ("A1,O22..") to object name("cloneDKO246","clone115","megaclone17"...)
#	• conversion of well_names in umi.tab ("PBZM0010382","PBZM0021007"..) to well_coordinate ("A1,O22..")
smt_technical_replicates_generator <- function (object_name, 
												umiTab_dir = "Data/",
												pd_dir="Data/",
												batch_umiTab,
												remove_ERCCs = FALSE,
												min_replicate = 1000)
{
wc_obj = read.table(paste0(pd_dir,"wellsCells_",object_name,".txt"),
				  sep="\t",header=T,colClasses=c(rep("character",2),rep("NULL",11)))		  
pd_obj = read.table(paste0(pd_dir,"plateDesign_",object_name,".txt"),sep="\t",header=TRUE)
umis_obj = read.table(paste0(umiTab_dir,batch_umiTab,".txt"),sep="\t")
wc_pd = data.frame(well = wc_obj$Well_ID,
				   meristem = pd_obj$ID[match(wc_obj$well_coordinates,pd_obj$well)])

grp_generator =rep(-1,nrow(wc_pd))
for (i in 1:length(grp_generator)){
grp_generator[i] = sum(wc_pd$meristem[1:i]==wc_pd$meristem[i])
}
wc_pd$grp_replicate = case_when(grp_generator%%2 != 0 ~ "replicate_group_A", grp_generator%%2 == 0 ~ "replicate_group_B")

order_meristems = match(colnames(umis_obj),wc_pd$well)
umis_obj_ordered = umis_obj[,order_meristems]

if (remove_ERCCs) {
umis_obj_ordered = umis_obj_ordered[!grepl("ERCC-00",rownames(umis_obj_ordered)),]
}
umisA = t(apply(umis_obj_ordered[,wc_pd$grp_replicate=="replicate_group_A"],1,function(x){
						tapply(x,wc_pd[wc_pd$grp_replicate=="replicate_group_A","meristem"],sum)}))
umisB = t(apply(umis_obj_ordered[,wc_pd$grp_replicate=="replicate_group_B"],1,function(x){
						tapply(x,wc_pd[wc_pd$grp_replicate=="replicate_group_B","meristem"],sum)}))
f_reps = colSums(umisA)>=min_replicate & colSums(umisB)>=min_replicate
return(list(repA = umisA[,f_reps], repB = umisB[,f_reps]))
}
################################################################################################
################################################################################################











############################################################################################################################
####### This function uses smt_technical_replicates_generator to load 2 groups of technical replicates for each meristem
smt_load_technical_replicates <- function (
								minimal_umi_per_replicate = 1e3,
								dst = TRUE,
								sft = TRUE,
								dstsft = TRUE,
								wt = TRUE,
								uf = TRUE,
								annotation_names = TRUE,
								keep_empty_wells = FALSE,
								blacklist_genes  = NULL,
								umi_tables_dir = "Data/",
								plates_design_dir = "Data/",
								blacklist_meristems = NULL,
								include_spike_ERCCs = TRUE,								
								annotation_csv_path = "Supplementary_Tables/tomato_gene_annotations.csv")


{

if (wt) {

print("loading plate27 replicates (wild-type)..")
p27_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate27_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0375",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)
	
print("loading plate28 replicates (wild-type)..")
p28_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate28_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0376",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate29 replicates (wild-type)..")
p29_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate29_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0377",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate30 replicates (wild-type)..")
p30_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate30_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0378",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)


print("loading plate36 replicates (wild-type)..")
p36_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate36_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0408",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate37 replicates (wild-type)..")
p37_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate37_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0409",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)


print("loading plate45 replicates (wild-type)..")
p45_rep = smt_technical_replicates_generator(
object_name="singleMers_plate45_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0418",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate46 replicates (wild-type)..")
p46_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate46_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0419",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate47 replicates (wild-type)..")
p47_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate47_wt_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0420",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)
} #wt



if (sft | dst){
print("loading plate22 replicates (sft n7187 mutants + dst mutants)..")
p22_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate22_sft_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0370",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate23 replicates (sft n7187 mutants + dst mutants)..")
p23_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate23_sft_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0371",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate24 replicates (sft n7187 mutants + dst mutants)..")
p24_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate24_sft_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0372",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate25 replicates (sft n7187 mutants + dst mutants)..")
p25_rep = smt_technical_replicates_generator(
object_name="singleMers_plate25_sft_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0373",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate26 replicates (sft n7187 mutants + dst mutants)..")
p26_rep = smt_technical_replicates_generator(
object_name="singleMers_plate26_sft_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0374",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)
} # dst|sft

if (dst) {
print("loading plate31 replicates (dst mutants)..")
p31_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate31_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0403",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate32 replicates (dst mutants)..")
p32_rep = smt_technical_replicates_generator(
object_name="singleMers_plate32_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0404",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate33 replicates (dst mutants)..")
p33_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate33_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0405",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate40 replicates (dst mutants)..")
p40_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate40_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0413",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate41 replicates (dst mutants)..")
p41_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate41_dst_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0414",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)} #dst


if (sft) {
print("loading plate34 replicates (sft n7187 mutants)..")
p34_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate34_sft_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0406",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate35 replicates (sft n7187 mutants)..")
p35_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate35_sft_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0407",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate38 replicates (sft n7187 mutants)..")
p38_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate38_sft_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0410",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate39 replicates (sft n7187 mutants)..")
p39_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate39_sft_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0411",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)} #sft



if (uf) {
print("loading plate42 replicates (uf mutants)..")
p42_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate42_uf_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0415",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate43 replicates (uf mutants)..")
p43_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate43_uf_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0416",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate44 replicates (uf mutants)..")
p44_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate44_uf_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0417",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)} #uf

if (dstsft) {
print("loading plate48 replicates (dst;sft double mutants)..")
p48_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate48_dstsft_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0421",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate)

print("loading plate49 replicates (dst;sft double mutants)..")
p49_rep = smt_technical_replicates_generator(
	object_name="singleMers_plate49_dstsft_dynabeads",
	umiTab_dir = umi_tables_dir,
	pd_dir = plates_design_dir,
	batch_umiTab="ABZM0422",
	remove_ERCCs=!include_spike_ERCCs,
	min_replicate = minimal_umi_per_replicate) #dst;sft
}
###############################################


print("merging.. ")
smer_repA=c();
smer_repB=c();

if (wt) {
smer_repA=cbind(smer_repA,p27_rep$repA,p28_rep$repA,p29_rep$repA,p30_rep$repA,p36_rep$repA,p37_rep$repA,p45_rep$repA,p46_rep$repA,p47_rep$repA)
smer_repB=cbind(smer_repB,p27_rep$repB,p28_rep$repB,p29_rep$repB,p30_rep$repB,p36_rep$repB,p37_rep$repB,p45_rep$repB,p46_rep$repB,p47_rep$repB)
}

if (sft | dst) {
smer_repA=cbind(smer_repA,p22_rep$repA,p23_rep$repA,p24_rep$repA,p25_rep$repA,p26_rep$repA)
smer_repB=cbind(smer_repB,p22_rep$repB,p23_rep$repB,p24_rep$repB,p25_rep$repB,p26_rep$repB)
} 

if (dst) {
smer_repA=cbind(smer_repA,p31_rep$repA,p32_rep$repA,p33_rep$repA,p40_rep$repA,p41_rep$repA)
smer_repB=cbind(smer_repB,p31_rep$repB,p32_rep$repB,p33_rep$repB,p40_rep$repB,p41_rep$repB)
}

if (sft) {
smer_repA=cbind(smer_repA,p34_rep$repA,p35_rep$repA,p38_rep$repA,p39_rep$repA)
smer_repB=cbind(smer_repB,p34_rep$repB,p35_rep$repB,p38_rep$repB,p39_rep$repB)
}

if (uf) {
smer_repA=cbind(smer_repA,p42_rep$repA,p43_rep$repA,p44_rep$repA)
smer_repB=cbind(smer_repB,p42_rep$repB,p43_rep$repB,p44_rep$repB)
}

if (dstsft) {
smer_repA=cbind(smer_repA,p48_rep$repA,p49_rep$repA)
smer_repB=cbind(smer_repB,p48_rep$repB,p49_rep$repB)
}

###################



### filter empty wells and bad genes
if (!keep_empty_wells){
	print("filtering empty wells.."); 
	f_empty =grepl("mpty|MPTY", colnames(smer_repA));
	smer_repA = smer_repA[,!f_empty];
	smer_repB = smer_repB[,!f_empty];
	}
if (!is.null(blacklist_genes)){print("filtering blacklisted genes..");
	smer_repA = smer_repA[rownames(smer_repA)[!rownames(smer_repA) %in% blacklist_genes],];
	smer_repB = smer_repB[rownames(smer_repB)[!rownames(smer_repB) %in% blacklist_genes],];
	}

#print(paste0("removing ",sum(colSums(smer_umi)<minimal_umi_per_meristem), " small meristems.."));
#smer_umi = smer_umi[,colSums(smer_umi)>=minimal_umi_per_meristem]

if (annotation_names) {
print("mishmashing gene names from different sources...")
gene_names_extended = smt_findAnnotations_from_solycID (solycIDs = rownames(smer_repA),
														 source_annotation_csvFile = annotation_csv_path)

rownames(smer_repA)=gene_names_extended;
rownames(smer_repB)=gene_names_extended;
}


########## removing blacklisted meristems (such as a pre-defined pipetting mixture)
if (!is.null(blacklist_meristems)) {
smer_repA = smer_repA[,!colnames(smer_repA) %in% blacklist_meristems]
smer_repB = smer_repB[,!colnames(smer_repB) %in% blacklist_meristems]}


########## normalizing and store objects in memory
print ("normalizing..")
smer_repA_n = apply(smer_repA,2,function(x) {(x/sum(x))*1e5 })
smer_repB_n = apply(smer_repB,2,function(x) {(x/sum(x))*1e5 })


print ("Downsampling..")
smer_repA_ds = smt_ds(smer_repA,minimal_umi_per_replicate)
smer_repB_ds = smt_ds(smer_repB,minimal_umi_per_replicate)
#browser()

#if (!keep_empty_wells) {rownames(smer_md) = smer_md$eppendorf_Num}

## storing data in permanent objects	
smer_repA <<- smer_repA;
smer_repB <<- smer_repB;
smer_repA_n <<- smer_repA_n;
smer_repB_n <<- smer_repB_n;
smer_repA_ds <<- smer_repA_ds;
smer_repB_ds <<- smer_repB_ds;

base_dir <<- "./"
myDate <<- substr(Sys.time(),1,10); 

cat(paste0("\nloaded objects:","\n",
			  "  smer_repA: raw expression table of odd technical replicates","\n",
			  "  smer_repA: raw expression table of even technical replicates","\n",
			  "  smer_repA/B_n: expression tables, normalized to 100K UMIs","\n",
			  "  smer_repA/B_ds: expression table, down-sampled to minimal threshold per replicate","\n",
			  ncol(smer_repA)," total meristems with both replicates above ",
			  minimal_umi_per_replicate," UMIs\n\n"))


} 
############################################################################################################################
############################################################################################################################




############################################################################################################################
# this function plots specific gene in technical replicates and show corr value
smt_plot_by_replicates <- function(gene_to_plot = NULL,
									gene_nm = NULL,
									svg=FALSE,
									repA_normalilzed_mat = smer_repA_n,
									repB_normalilzed_mat = smer_repB_n,
									pp=72,
									h = 250,
									w = 250,
									log2_exp=FALSE,
									log2_reg = 1,
									add_equator = FALSE,
									meristems_metadata = NULL,
									spearman_cor = TRUE,
									dots_size = 0.8,
									cex_txt = 1.5,
									myDate = substring(Sys.time(),1,10),
									base_dir = "./",
									plot_margins = c(5,5,1,1),
									show_cor = TRUE) {

if (is.null(gene_to_plot) | length(gene_to_plot)!=1) { stop("One gene must be specified for this replicates plot") }
if (is.null(gene_nm) | length(gene_nm)!=1) { stop("One gene label/name be specified for this replicates plot") }
									
expA = repA_normalilzed_mat[gene_to_plot,]
expB = repB_normalilzed_mat[gene_to_plot,]

if (log2_exp) {
expA = log2(expA+log2_reg)
expB = log2(expB+log2_reg)
}

#load metadata for coloring of meristems:
if (is.null(meristems_metadata)) {
md_table = read.csv("Data/SMT_metadata_season2.csv",header=TRUE,row.names=1)
md_table_filt = md_table[match(colnames(repA_normalilzed_mat),rownames(md_table)),]
md_table_filt$color=NA
md_table_filt$color[md_table_filt$dev_stage=="TM"]="grey"
md_table_filt$color[md_table_filt$dev_stage=="VM"]="green"
md_table_filt$color[md_table_filt$dev_stage=="TM0"]="green3"
md_table_filt$color[md_table_filt$dev_stage=="TM1"]="yellow"
md_table_filt$color[md_table_filt$dev_stage=="TM2"]="orange"
md_table_filt$color[md_table_filt$dev_stage=="TM2E"]="orange"
md_table_filt$color[md_table_filt$dev_stage=="FM"]="red"
md_table_filt$color[md_table_filt$dev_stage=="LTM"]="red"
md_table_filt$color[md_table_filt$dev_stage=="EFM"]="purple"
meristems_metadata = md_table_filt}

plot_lims = range(c(expA,expB))
replicate_cors = round(cor(expB,expA,method=ifelse(spearman_cor,"spearman","pearson")),digits=2)
relicate_sign = ifelse(spearman_cor,intToUtf8(0x03C1),"r")
cat(paste0("Plotting replicates of ",gene_nm,"\n"))
if (!svg){
png(paste0(base_dir,myDate,"_",gene_nm,"_tech_replicates.png"),width=w,height=h)} else {
svglite(paste0(base_dir,myDate,"_",gene_nm,"_tech_replicates.svg"),width=w/pp,height=h/pp)} 

par(mar=plot_margins,bty="L")
plot(expA,expB,
	 col=meristems_metadata[names(expA),"color"],
	 pch=19,
	 cex=dots_size,
	 cex.axis = cex_txt,
	 cex.lab = cex_txt,
	 xlim = plot_lims,
	 ylim = plot_lims,
	 xlab = paste0(gene_nm," (replicate #1)"),
	 ylab = paste0(gene_nm," (replicate #2)"))
if (show_cor) {
legend(x="topleft",
	  legend=paste0(relicate_sign,"=",replicate_cors),
	   bty="n",cex=cex_txt,x.intersp=0.2) 
	   }
if (add_equator){abline(a=0,b=1,lty=2,col="#202020",lwd=1)}
dev.off()
}	 
############################################################################################################################
############################################################################################################################



#################################################################################################################
## This function downsample a umi_tab to a specified value
###########################################################
smt_ds <- function (umi_tab, ds_val)
{
umi_tab_filt = umi_tab[,colSums(umi_tab)>=ds_val]
umi_tab_ds = apply(umi_tab_filt, 2, .downsamp_one, ds_val)
rownames(umi_tab_ds)=rownames(umi_tab);
return(umi_tab_ds)
}
#################################################################################################################
.downsamp_one=function(v,n)
{
  hist(sample(rep(1:length(v),times=v),replace=F,size=n),0.5+0:length(v),plot=F)$counts
}
##############################################################################################
##############################################################################################



##############################################################################################
smt_plot_knn_selection <- function(exp_n = smer_n,
								   knn_matrix = NULL,
								   svg = FALSE,
								   pp = 72,
								   cex_selected_meristem = 1.4,
								   cex_knn = 0.7,
								   col_selected_meristem = "blue",
								   col_knn = "#ff000080",
								   plot_margins = c(5,5,4,4),
								   h = 280,
								   w = 280)
								   {
								   
								   
if (is.null(knn_matrix)) {stop ("Please provide a KNN matrix in order to plot example")}
if (nrow(knn_matrix)!=ncol(exp_n)) {stop ("Provide a valid KNN matrix in order to plot example")}
 

if (svg) {
svglite(paste0(base_dir,myDate,"_KNN_selection_example.svg"),width = w/pp, height = h/pp) 
	} else {
png(paste0(base_dir,myDate,"_KNN_selection_example.png"),width = w, height = h) }

par(mar=plot_margins,bty="L")
example_mer = rownames(knn_matrix)[sample(1:nrow(knn_matrix),size=1)]
plot(colSums(exp_n[mod_lateral_dissection_1,]),
	 colSums(exp_n[mod_lateral_dissection_2,]),
	 col="grey63",
	 xlab="non-proliferation",
	 ylab="proliferation",
	 cex=cex_knn,
	 cex.lab=1.5,cex.axis=1.5,
	 pch=19)
points(sum(exp_n[mod_lateral_dissection_1,example_mer]),
	   sum(exp_n[mod_lateral_dissection_2,example_mer]),
	   col="blue",
	   cex=cex_selected_meristem,
	   pch=19)
points(colSums(smer_n[mod_lateral_dissection_1,colnames(knn_matrix)[which(knn_matrix[example_mer,]==1)]]),
	   colSums(smer_n[mod_lateral_dissection_2,colnames(knn_matrix)[which(knn_matrix[example_mer,]==1)]]),
	   col="#ff000080",
	   cex=cex_knn,
	   pch=19)
legend(x="topright",
		cex=1.5,
		legend=paste0("K=",as.character(sum(knn_matrix[example_mer,]))),
		bty="n",
		text.col="red",
		x.intersp=0.2,
		y.intersp=0.8)
dev.off()	   
}
##############################################################################################
##############################################################################################





##############################################################################################
smt_batch_normalization <<- function (mat_n = NULL, 
									  meristems = NULL, 
									  k_meristems = NULL) {
smt_load_common_nms()
if (is.null(mat_n)) {stop ("Please provide expression matrix that will serve to KNN-based batch normalization")}
if (is.null(meristems)) {stop ("Please provide a subset of meristems you wish to normalize by one another")}
if (is.null(k_meristems)) {stop ("Please provide K - how many neighboring meristems to normalize by")}
prol_feats = data.frame(non_proliferation = colSums(mat_n[mod_lateral_dissection_1,meristems]),
						proliferation = colSums(mat_n[mod_lateral_dissection_2,meristems])) %>%						
						t() %>% t() %>%
						scale(center=TRUE, scale=TRUE)
knn_matrix = smt_get_norm_knn_matrix(feats = prol_feats,
						k = k_meristems)	
norm_smer_n = smt_normalize_knn (raw = mat_n[,meristems],
								  knn_mat = knn_matrix)

return(norm_smer_n)
}
##############################################################################################
##############################################################################################

###### this are the two base functions that substract expression by knn
# select k meristems most similar to each one along specific featsn
##############################################################################################
smt_get_norm_knn_matrix <- function(feats, k){
    dist_mat <- dist(feats)
    knn_df <- tgs_knn(100-as.matrix(dist_mat), k)
    knn_df <- knn_df %>% mutate(col1 = factor(col1), col2 = factor(col2, levels=levels(col1)))
    
    knn_mat <- sparseMatrix(as.numeric(knn_df$col1), as.numeric(knn_df$col2), x=1)
    rownames(knn_mat) <- levels(knn_df$col1)
    colnames(knn_mat) <- levels(knn_df$col2)
    return(knn_mat)
}
##############################################################################################
##############################################################################################

# normalize expression of each gene in each meristems by substracting the average expression in its k nearest meristems
##############################################################################################
smt_normalize_knn <- function(raw, knn_mat){
    raw <- raw[, colnames(knn_mat)]
    raw_filt_na <- raw
    raw_filt_na[is.na(raw)] <- 0
    met_exp <- as.matrix(raw_filt_na) %*% as.matrix(t(knn_mat))
    not_na_mat <- !is.na(raw)
    met_exp_n <- as.matrix(not_na_mat) %*% as.matrix(t(knn_mat))
    met_exp_norm <- met_exp / met_exp_n
    met_oe <- raw - met_exp_norm
    return(met_oe)
}
##############################################################################################
##############################################################################################





##############################################################################################
smt_order_mutant_by_wt <<- function(wt_norm_ordered_mat, 
									mut_norm_mat, 
									corr_thresh, 
									mrkr_genes) {
chimeric_mat = t(cor(wt_norm_ordered_mat[mrkr_genes,],
				   mut_norm_mat[mrkr_genes,],
				  method = "spearman"))
chimeric_mat_positive = chimeric_mat
chimeric_mat_positive[chimeric_mat_positive<(corr_thresh)] = min(c(0,corr_thresh))
mut_orders <- slanted_orders(chimeric_mat_positive, order_cols=FALSE)
chimeric_mat_reord = chimeric_mat[mut_orders$rows, mut_orders$cols]  				 

chimeric_mat_mutTree <- slanter::oclust(dist(chimeric_mat_reord))

return(list (corr_mat = chimeric_mat_reord, mutTree = chimeric_mat_mutTree))
}
##############################################################################################
##############################################################################################











### this function gets two vectors of at least one gene in  each,
### and plot their normalized expression per 100K UMIs in each meristem, colored by its stage
##############################################################################################
smt_plot_xy <- function(genes_a = NULL,
						 genes_b = NULL,
						 genes_a_label = "label A",
						 genes_b_label = "label B",						 
						 log2_exp = FALSE,
						 color_by_dev=TRUE,
						 svg=FALSE,
						 pp=72,
						 reg_log=1,
						 normalized_expression_mat = smer_n,
						 ignore_meristems=NULL,
						 scatter_margins = c(5,5,1,1),
						 scatter_width=280,
						 scatter_height=280,
						 base_dir = "./",
						 myDate = substr(Sys.time(),1,10),
						 show_legend=FALSE,
						 legend_pos="topright",
						 legend_cex=1.5,
						 my_xlim=NULL,
						 my_ylim=NULL,
						 my_cex=0.8,
		 				 save_to_file=FALSE,
						 test_meristems=FALSE){
						 
if (is.null (genes_a) | is.null (genes_b)) {print ("no genes inserted!!\n"); break;}




### keep this option to replace normalized mat
smer_n = normalized_expression_mat;
###

	if (!is.null(ignore_meristems)) {
	smer_n=smer_n[,!colnames(smer_n) %in% ignore_meristems]
	smer_umi=smer_umi[,!colnames(smer_umi) %in% ignore_meristems]
	smer_md=smer_md[!smer_md$eppendorf_Num %in% ignore_meristems,]
	}

if (!log2_exp){
	#if (!is.null(my_xlim)){xlim_cur=my_xlim} else{
	#					   xlim_cur = max(log2()
	if (save_to_file){
	if (svg) { svg(paste0(base_dir,myDate,paste0(genes_a_label,"_",genes_b_label,".svg")),width=scatter_width/pp, height=scatter_height/pp)
		} else { png(paste0(base_dir,myDate,paste0(genes_a_label,"_",genes_b_label,".png")),width=scatter_width, height=scatter_height) }
	par(mar=scatter_margins)}
	
	if(length(genes_a)==1 & length(genes_b)==1) {
	plot(
		smer_n[genes_a,],
		smer_n[genes_b,],
		col=list(smer_md$color,smer_md$color_batch)[[ifelse(color_by_dev,1,2)]],
		pch=19,
		cex.axis=1.5,cex.lab=1.5,
		xlab=genes_a_label,
		ylab=genes_b_label,
		xlim=list( c(min(c(smer_n[genes_a,])), max(c(smer_n[genes_a,]))),
					my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
		ylim=list( c(min(c(smer_n[genes_b,])), max(c(smer_n[genes_b,]))),
					my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
					
		cex=my_cex)
	} else if (length(genes_a) > 1 & length(genes_b) == 1){
	plot(
		colSums(smer_n[genes_a,]),
		smer_n[genes_b,],
		col=list(smer_md$color,smer_md$color_batch)[[ifelse(color_by_dev,1,2)]],
		pch=19,
		cex.axis=1.5,cex.lab=1.5,
		xlab=genes_a_label,
		ylab=genes_b_label,
		xlim=list( c(min(c(colSums(smer_n[genes_a,]))), max(c(colSums(smer_n[genes_a,])))),
			my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
		ylim=list( c(min(c(smer_n[genes_b,])), max(c(smer_n[genes_b,]))),
			my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
		cex=my_cex)
	} else if (length(genes_a) == 1 & length(genes_b) > 1){
		plot(
		smer_n[genes_a,],
		colSums(smer_n[genes_b,]),
		col=list(smer_md$color,smer_md$color_batch)[[ifelse(color_by_dev,1,2)]],
		pch=19,
		cex.axis=1.5,cex.lab=1.5,
		xlab=genes_a_label,
		ylab=genes_b_label,
		xlim=list( c(min(c(smer_n[genes_a,])), max(c(smer_n[genes_a,]))),
			my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
		ylim=list( c(min(c(colSums(smer_n[genes_b,]))), max(c(colSums(smer_n[genes_b,])))),
			my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
		cex=my_cex)
	} else {
		plot(
		colSums(smer_n[genes_a,]),
		colSums(smer_n[genes_b,]),
		col=list(smer_md$color,smer_md$color_batch)[[ifelse(color_by_dev,1,2)]],
		pch=19,
		cex.axis=1.5,cex.lab=1.5,
		xlab=genes_a_label,
		ylab=genes_b_label,
		xlim=list( c(min(c(colSums(smer_n[genes_a,]))), max(c(colSums(smer_n[genes_a,])))),
					my_xlim)[[ifelse(is.null(my_xlim),1,2)]], 
		ylim = list( c(min(c(colSums(smer_n[genes_b,]))), max(c(colSums(smer_n[genes_b,])))),
					my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
		cex=my_cex)	
	}
	grid()
	if (show_legend & color_by_dev) {
		legend(x=legend_pos,
				x.intersp=0.2,
				y.intersp=0.8,
				bty="n",
				fill=list(c("green","green3","yellow","orange","red","purple"),c("green","green3","yellow","orange","grey","red"))[[ifelse(test_meristems,2,1)]],
				legend=list(c("VM","TM0","TM1","TM2","LTM","EFM"),c("VM","TM0","TM1","TM2","TM","FM"))[[ifelse(test_meristems,2,1)]],
				cex=legend_cex)
					
				
		}else if (show_legend & !color_by_dev) {
		legend(x=legend_pos,
				x.intersp=0.2,
				y.intersp=0.8,
				bty="n",
				fill=c("yellow","gold","orange","green","red","purple","grey","black","lightblue","blue","pink","navajowhite3","cyan","darkgrey","brown"),
				legend=c(paste0("batch",c(0:14))),
				cex=legend_cex)		
		}
	if (save_to_file){
	dev.off()}

} else {					 

	if (save_to_file){
	if (svg) { svg(paste0(base_dir,myDate,paste0("_log2_",genes_a_label,"_",genes_b_label,".svg")), width=scatter_width/pp, height=scatter_height/pp)
		} else { png(paste0(base_dir,myDate,paste0("_log2_",genes_a_label,"_",genes_b_label,".png")),width=scatter_width,height=scatter_height) }
	par(mar=scatter_margins)}
	if(length(genes_a)==1 & length(genes_b)==1) {
	plot(
		log2(smer_n[genes_a,]+reg_log),
		log2(smer_n[genes_b,]+reg_log),
		col=list(smer_md$color,smer_md$color_batch)[[ifelse(color_by_dev,1,2)]],
		pch=19,
		cex.axis=1.5,cex.lab=1.5,
		xlab=genes_a_label,
		ylab=genes_b_label,
		cex=my_cex)
	} else if (length(genes_a) > 1 & length(genes_b) == 1){
	plot(
		log2(colSums(smer_n[genes_a,])+reg_log),
		log2(smer_n[genes_b,]+reg_log),
		col=list(smer_md$color,smer_md$color_batch)[[ifelse(color_by_dev,1,2)]],
		pch=19,
		cex.axis=1.5,cex.lab=1.5,
		xlab=genes_a_label,
		ylab=genes_b_label,
		cex=my_cex)	
	} else if (length(genes_a) == 1 & length(genes_b) > 1){
		plot(
		log2(smer_n[genes_a,]+reg_log),
		log2(colSums(smer_n[genes_b,])+reg_log),
		col=list(smer_md$color,smer_md$color_batch)[[ifelse(color_by_dev,1,2)]],
		pch=19,
		cex.axis=1.5,cex.lab=1.5,
		xlab=genes_a_label,
		ylab=genes_b_label,
		cex=my_cex)
	} else {
		plot(
		log2(colSums(smer_n[genes_a,])+reg_log),
		log2(colSums(smer_n[genes_b,])+reg_log),
		col=list(smer_md$color,smer_md$color_batch)[[ifelse(color_by_dev,1,2)]],
		pch=19,
		cex.axis=1.5,cex.lab=1.5,
		xlab=genes_a_label,
		ylab=genes_b_label,
		cex=my_cex)
	}
	grid()

	if (show_legend & color_by_dev) {
		legend(x=legend_pos,
				x.intersp=0.2,
				y.intersp=0.8,
				bty="n",
				fill=list(c("green","darkgreen","yellow","orange","red","purple"),c("green","darkgreen","yellow","orange","grey","red"))[[ifelse(test_meristems,2,1)]],
				legend=list(c("VM","TM0","TM1","TM2","LTM","EFM"),c("VM","TM0","TM1","TM2","TM","FM"))[[ifelse(test_meristems,2,1)]],
				cex=legend_cex)
				
				
		} else if (show_legend & !color_by_dev) {
		legend(x=legend_pos,
				x.intersp=0.2,
				y.intersp=0.8,
				bty="n",
				fill=c("yellow","gold","orange","green","red","purple","grey","black","lightblue","blue","pink","navajowhite3","cyan","darkgrey","brown"),
				legend=c(paste0("batch",c(0:14))),
				cex=legend_cex)		
		}
	if(save_to_file){	
	dev.off()}

}

}
##############################################################################################
##############################################################################################





#### this function is a modification for pheatmap, where center0 makes 0 values as white
##############################################################################################
smt_mypheatmap = function (mat, w=1000, h=1000, fname,outdir, center0 = F, svg = FALSE, svg_ppi=72,
							my_color = colorRampPalette( rev(brewer.pal(n = 7, name ="RdBu")))(100),
							external_breaks_for_center0 = NULL,
							color_quantile = 0.99, ...) {
  require(svglite)
  parms = list(... )
  if (svg) { svglite(sprintf("%s%s", outdir, fname), width=w/svg_ppi, height=h/svg_ppi) } else {
  png(sprintf("%s%s", outdir, fname), width=w, height=h)}
  if (center0) {
    if (is.null(parms$color)) {
      parms$color =  my_color
    }
    mm = quantile(abs(mat), color_quantile,na.rm=T)
    parms$breaks = list(sort(unique (c( min(mat), seq(-mm, mm, l=99),max(mat)))),external_breaks_for_center0)[[ifelse(is.null(external_breaks_for_center0),1,2)] ]
  }
   parms$color =  my_color
  parms$mat = mat
  do.call (pheatmap, parms) 
 dev.off()
}
##############################################################################################
##############################################################################################



#### This function plot expression of one gene (or sum of several genes) along a suggested pseudo-time provided to it (ordered_meristems)
##############################################################################################
smt_plot_by_ord_season2 <- function(
							mat_n = smer_n,
							mat = smer_umi,
							base_dir,gene_to_plot,gene_nm,
							k_roll=11,
							log2_exp=FALSE,log_reg=1,
							trend_col="black",
							manual_colors = NULL,
							ignore_meristems = NULL,
							my_xlab="ordered meristems",
							ylab_col = "black",
							ylab_dist = 1.8,
							txt_size = 1.8,
							plot_axs = "r",
							plot_margins = c(8,5,1,1),
							max_umi_meristem=NULL,
							bg_vec=NULL,
							bg_group_colors=NULL,
							bg_group_fracs=NULL,
							remove_frame = FALSE,
							k_cex=0.7,
							my_cex=1,
							save_to_file=TRUE,
							break_by_groups = NULL,
							my_ylim=NULL,
							values_las = 2,
							svg=FALSE,
							pp=72,
							plot_grid = FALSE,
							figure_main="",
							ordered_meristems = NULL,
							dev_colorbar=FALSE,xlab_dist=3.5,my_height=300,my_width=500,
							show_colorbar=FALSE,image_mult_cols=0.1) {

if (is.null(ordered_meristems)){ stop ("insert ordered vector of meristem names (ordered_meristems parameter")}
							
if(!is.null(max_umi_meristem)) { 
mat_n = apply(smer_squeeze_meristems(mat=mat,max_umi=max_umi_meristem),2,function(x){(x/sum(x))*1e5})
}

if(!is.null(ignore_meristems)) { 
mat_n = mat_n[,colnames(mat_n)[!colnames(mat_n) %in% ignore_meristems]]
smer_md = smer_md[!smer_md$eppendorf_Num %in% ignore_meristems,]
ordered_meristems = ordered_meristems[!ordered_meristems %in% ignore_meristems]
}


if (!is.null(manual_colors)) {
	comb_colors = manual_colors } else {
	comb_colors = smer_md$color[match(ordered_meristems,smer_md$eppendorf_Num)] 
}

if (!is.null(bg_vec)){
bg_vec = bg_vec[ordered_meristems]
bg_color = unique(bg_vec)
bg_frac = rep(-1,length(unique(bg_vec)))
	for (i in 1:length(bg_frac)){
	bg_frac[i] = sum( (bg_vec %in% unique(bg_vec)[i]) /length(bg_vec)) 
	}
}

if (!is.null(bg_group_colors) & !is.null(bg_group_fracs)){
bg_color = bg_group_colors
bg_frac = bg_group_fracs
}

if (is.null(break_by_groups)) {
xs = 1:length(ordered_meristems)
	} else {
	xgroups = length(break_by_groups)+1
	space_group = length(ordered_meristems) / xgroups
	xpos = c();
		for (i in 1:xgroups){
	
			if (i==1) {group_n = break_by_groups[1]; current_range= c(0,space_group)}
			else if (i==xgroups) {group_n = length(ordered_meristems) - break_by_groups[xgroups-1];
								current_range = c((i*space_group)-space_group,i*space_group); 
			} else {
			current_range = c((i*space_group)-space_group,i*space_group); 
			group_n = break_by_groups[i] - break_by_groups[i-1]
			}
		xpos = c(xpos,seq(current_range[1],current_range[2],length.out = group_n))
		}
	xs = xpos;
	}


if (length(gene_to_plot)==1){
if (!dev_colorbar){

if (save_to_file) {
if (!svg) {
png(paste0(base_dir,gsub("/","_",gene_nm),"_sorted_meristems.png"),width=my_width,height=my_height) 
	} else {
svglite(paste0(base_dir,gsub("/","_",gene_nm),"_sorted_meristems.svg"),width=my_width/pp,height=my_height/pp) 
	}
}
	
par(mar = plot_margins)
ordered_exp = list(c(mat_n[gene_to_plot,ordered_meristems]),c(log2(mat_n[gene_to_plot,ordered_meristems]+log_reg)) )[[ifelse(log2_exp,2,1)]]
plot(xs,ordered_exp,pch=19,cex=my_cex,
	 col=comb_colors,
	  main=figure_main,
	  cex.main = txt_size,
	  xlab="",
	 xaxs = plot_axs,
	 las = values_las,
	 frame.plot = !remove_frame,
	 ylab="",
	 			ylim=list( c(min(ordered_exp), max(ordered_exp)),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
	 cex.axis=txt_size,
	 cex.lab=txt_size,
	 xaxt="n",
	 yaxt="n")	 
title(xlab=my_xlab,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=xlab_dist)
title(ylab=gene_nm,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=ylab_dist, col.lab = ylab_col)

if (!remove_frame){	 
axis(1,labels=FALSE,cex.axis=txt_size)}
yticks =seq(min(ordered_exp),max(ordered_exp),length.out=5)
yticks_round = round(yticks,digits=(-3))
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=(-2))}
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=(-1))}
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=0)}
magaxis(side=2,majorn=4,minorn="auto",
		frame.plot=!remove_frame,cex.axis=txt_size,cex=txt_size,las=values_las)

points(xs,rollmean(ordered_exp,k=k_roll,na.pad=TRUE),pch=19,cex=k_cex,col=trend_col) 
if (plot_grid) {grid()}

	if(!is.null(bg_vec) | (!is.null(bg_group_colors) & !is.null(bg_group_fracs))){
	bg_lims= cumsum(c(par("usr")[1],cumsum(par("usr")[2]-par("usr")[1])*bg_frac))
			for (i in 1:length(bg_color)){
			rect(bg_lims[i], par("usr")[3], bg_lims[i+1], par("usr")[4], col = bg_color[i], border =NA)
			}
	}

if (save_to_file) {	dev.off()     }

} else {

if (save_to_file) {
if (!svg) {
	png(paste0(base_dir,gsub("/","_",gene_nm),"_sorted_meristems.png"),width=my_width,height=my_height) 
	} else {
	svglite(paste0(base_dir,gsub("/","_",gene_nm),"_sorted_meristems.svg"),width=my_width/pp,height=my_height/pp) 
	}
}

par(mar=plot_margins)
plot(xs,
	rollmean(mat_n[gene_to_plot,ordered_meristems],k=k_roll,na.pad=TRUE),
	type="l",lwd=4,col=trend_col, 
	xaxs = plot_axs,
	 xlab="",
	 frame.plot = !remove_frame,
	  main=figure_main,
	  cex.main = txt_size,
	  las = values_las,
	 ylab="",
	 cex.axis=txt_size,
	 cex.lab=txt_size,
	 xaxt="n",
	 	 			ylim=list( c(min(mat_n[gene_to_plot,ordered_meristems]), max(mat_n[gene_to_plot,ordered_meristems])),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
	 yaxt="n")
	if(!is.null(bg_vec) | (!is.null(bg_group_colors) & !is.null(bg_group_fracs))){
bg_lims= cumsum(c(par("usr")[1],cumsum(par("usr")[2]-par("usr")[1])*bg_frac))
		for (i in 1:length(bg_color)){
		rect(bg_lims[i], par("usr")[3], bg_lims[i+1], par("usr")[4], col = bg_color[i], border = bg_color[i])
		}
}

if (!remove_frame){
axis(side=1,at=1:length(comb_colors),labels=FALSE,cex.axis=txt_size,tck=-0.02)}
yticks =seq(min(ordered_exp),max(ordered_exp),length.out=5)
yticks_round = round(yticks,digits=(-3))
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=(-2))}
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=(-1))}
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=0)}
magaxis(side=2,majorn=4,minorn="auto",
		frame.plot=!remove_frame,cex.axis=txt_size,cex=txt_size,las=values_las)


title(xlab=my_xlab,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=xlab_dist)
title(ylab=gene_nm,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=ylab_dist, col.lab = ylab_col)

if (plot_grid) {grid()}
	if (show_colorbar) {
	counter_pos = 1;
		for (i in 1:length(comb_colors)){
		segments(xpd=TRUE,
				x0=counter_pos-1,
				 y0=ifelse(comb_colors[i]=="green",-5*image_mult_cols,
					ifelse(comb_colors[i]=="darkgreen",-8*image_mult_cols,
					ifelse(comb_colors[i]=="yellow",-11*image_mult_cols,
					ifelse(comb_colors[i]=="orange",-14*image_mult_cols,-17*image_mult_cols)))),
				 x1=counter_pos,
				 y1=ifelse(comb_colors[i]=="green",-5*image_mult_cols,
					ifelse(comb_colors[i]=="darkgreen",-8*image_mult_cols,
					ifelse(comb_colors[i]=="yellow",-11*image_mult_cols,
					ifelse(comb_colors[i]=="orange",-14*image_mult_cols,-17*image_mult_cols)))),
				 col=comb_colors[i],
					lwd=8)
		counter_pos=counter_pos+1;
		}
	}

	
if (save_to_file) {	dev.off()     }
}



## if there's more than one gene to plot along the ordered meristems (so we sum their output)
}else{

if (!dev_colorbar){

if (save_to_file) {
if (!svg) {
png(paste0(base_dir,gsub("/","_",gene_nm[1],"and_others"),"_sorted_meristems.png"),width=my_width,height=my_height) 
	} else {
svglite(paste0(base_dir,gsub("/","_",gene_nm[1],"and_others"),"_sorted_meristems.svg"),width=my_width/pp,height=my_height/pp) 	
	}
}	

par(mar=plot_margins)
ordered_exp = list(c(colSums(mat_n[gene_to_plot,ordered_meristems])),c(log2(colSums(mat_n[gene_to_plot,ordered_meristems])+log_reg)) )[[ifelse(log2_exp,2,1)]]
plot(xs,
	ordered_exp,pch=19,cex=my_cex,col=comb_colors,
	 xlab="",
	 ylab="",
	  main=figure_main,
	  cex.main = txt_size,
	 xaxs = plot_axs,
	 frame.plot = !remove_frame,
 	 las = values_las,
	 cex.axis=txt_size,
	 cex.lab=txt_size,
	 ylim=list( c(min(ordered_exp), max(ordered_exp)),
	 my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
	 xaxt="n",
	 yaxt="n")
title(xlab=my_xlab,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=xlab_dist)
title(ylab=gene_nm,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=ylab_dist, col.lab = ylab_col)	 

	if(!is.null(bg_vec) | (!is.null(bg_group_colors) & !is.null(bg_group_fracs))){
bg_lims= cumsum(c(par("usr")[1],cumsum(par("usr")[2]-par("usr")[1])*bg_frac))
		for (i in 1:length(bg_color)){
		rect(bg_lims[i], par("usr")[3], bg_lims[i+1], par("usr")[4], col = bg_color[i], border = bg_color[i])
		}
}
if (!remove_frame){
axis(1,labels=FALSE,cex.axis=txt_size)}
yticks =seq(min(ordered_exp),max(ordered_exp),length.out=5)
yticks_round = round(yticks,digits=(-3))
if(sum(duplicated(yticks_round))>0){
yticks_round = round(yticks,digits=(-2))}
if(sum(duplicated(yticks_round))>0){
yticks_round = round(yticks,digits=(-1))}
if(sum(duplicated(yticks_round))>0){
yticks_round = round(yticks,digits=0)}
yticks_round = seq(min(yticks_round),max(yticks_round),length.out=5)
magaxis(side=2,
		majorn=4,
		minorn="auto",
		frame.plot=!remove_frame,cex.axis=txt_size,cex=txt_size,las=values_las)
points(xs,rollmean(list(colSums(mat_n[gene_to_plot,ordered_meristems]),log2(colSums(mat_n[gene_to_plot,ordered_meristems])+log_reg))[[ifelse(log2_exp,2,1)]],k=k_roll,na.pad=TRUE),pch=19,cex=k_cex,col=trend_col) 
if (plot_grid) {grid()}

if (save_to_file) { dev.off()     }

} else {

if (save_to_file) {
if (!svg) {
	png(paste0(base_dir,gsub("/","_",gene_nm),"_sorted_meristems.png"),width=my_width,height=my_height)
	} else {
	svglite(paste0(base_dir,gsub("/","_",gene_nm),"_sorted_meristems.svg"),width=my_width/pp,height=my_height/pp)
	}
}

par(mar=plot_margins)
plot(xs,
	 rollmean(colSums(mat_n[gene_to_plot,ordered_meristems]),k=k_roll,na.pad=TRUE),type="l",lwd=4,col=trend_col, 
	 xlab="",
	 frame.plot = !remove_frame,
	 ylab="",
	  main=figure_main,
	  cex.main = txt_size,
	 xaxs = plot_axs,
 	 las = values_las,
	 cex.axis=txt_size,
	 cex.lab=txt_size,
	 	 			ylim=list( c(min(colSums(mat_n[gene_to_plot,ordered_meristems])), max(colSums(mat_n[gene_to_plot,ordered_meristems]))),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
	 xaxt="n",
	 yaxt="n")

	
	if(!is.null(bg_vec) | (!is.null(bg_group_colors) & !is.null(bg_group_fracs))){
		bg_lims= cumsum(c(par("usr")[1],cumsum(par("usr")[2]-par("usr")[1])*bg_frac))
				for (i in 1:length(bg_color)){
				rect(bg_lims[i], par("usr")[3], bg_lims[i+1], par("usr")[4], col = bg_color[i], border = bg_color[i])
				}
		}

if (!remove_frame){
axis(side=1,at=1:length(comb_colors),labels=FALSE,cex.axis=txt_size,tck=-0.02)}
yticks =seq(min(ordered_exp),max(ordered_exp),length.out=5)
yticks_round = round(yticks,digits=(-3))
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=(-2))}
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=(-1))}
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=0)}
magaxis(
		side=2,majorn=4,minorn="auto",
		frame.plot=!remove_frame,cex.axis=txt_size,cex=txt_size,las=values_las)


title(xlab=my_xlab,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=xlab_dist)
title(ylab=gene_nm,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=ylab_dist, col.lab = ylab_col)
if (plot_grid) {grid()}
	if (show_colorbar) {
	counter_pos = 1;
		for (i in 1:length(comb_colors)){
		segments(xpd=TRUE,
				x0=counter_pos-1,
				 y0=ifelse(comb_colors[i]=="green",-5*image_mult_cols,
					ifelse(comb_colors[i]=="darkgreen",-8*image_mult_cols,
					ifelse(comb_colors[i]=="yellow",-11*image_mult_cols,
					ifelse(comb_colors[i]=="orange",-14*image_mult_cols,-17*image_mult_cols)))),
				 x1=counter_pos,
				 y1=ifelse(comb_colors[i]=="green",-5*image_mult_cols,
					ifelse(comb_colors[i]=="darkgreen",-8*image_mult_cols,
					ifelse(comb_colors[i]=="yellow",-11*image_mult_cols,
					ifelse(comb_colors[i]=="orange",-14*image_mult_cols,-17*image_mult_cols)))),
				 col=comb_colors[i],
					lwd=8)
		counter_pos=counter_pos+1;
		}
	}

if (save_to_file) {
dev.off()     }
}


} # more than one gene
} # function
####################################################################################################
####################################################################################################





# This function assign a WT-sister meristem to each mutant meristem, based on meristem-meristem correlations
####################################################################################################
smt_find_wt_sisters <- function (chimeric_ordered_correlation_matrix = NULL,
								 n_top_meristems = 3) { 

if (is.null(chimeric_ordered_correlation_matrix)){ stop ("Provide chimeric, ordered meristem-meristem correlation matrix")}
if (n_top_meristems %% 2 != 1) { stop("in order to select 'median meristem', n_top_meristems must be an odd integer")}					 
top_wts = t(apply(chimeric_ordered_correlation_matrix,1,function(x) {
			match(names(tail(sort(x),n_top_meristems)),colnames(chimeric_ordered_correlation_matrix)) } ))
top_wts_median = apply(top_wts,1,median)
names(top_wts_median) = rownames(chimeric_ordered_correlation_matrix)

return (top_wts_median)
}
####################################################################################################
####################################################################################################




# This function plot examples of how WT-sister meristems are assigned to each mutant meristem, based on meristem-meristem correlations
####################################################################################################
smt_plot_wt_sister_assigment <- function (chimeric_ordered_correlation_matrix = NULL,
										  n_top_meristems = 3,
										  assigned_wt_sisters = NULL,
										  meristems_metadata = smer_md,
										  wt_phase_colors = NULL,
										  w = 600,
										  h = 200,
										  wt_phase_assigment = NULL,
										  mutant_genotype_name = "",
										  base_dir = "./",
										  svg=FALSE,
										  ppi=72,
										  example_meristems_to_plot = 5,
										  myDate = substring(Sys.time(),1,10)) {

if (is.null(chimeric_ordered_correlation_matrix)){ stop ("Provide chimeric, ordered meristem-meristem correlation matrix")}
if (is.null(assigned_wt_sisters)){ stop ("Provide a vector of assigned WT sisters")}
								  
top_wts = t(apply(chimeric_ordered_correlation_matrix,1,function(x) {
			match(names(tail(sort(x),n_top_meristems)),colnames(chimeric_ordered_correlation_matrix)) } ))

if (!is.null (wt_phase_assigment) & !is.null (wt_phase_colors)) {
background_limits = c(1,cumsum(table(wt_phase_assigment))) }


			
mutant_meristem_to_plot = sample(1:nrow(chimeric_ordered_correlation_matrix),size=example_meristems_to_plot,replace=FALSE)
sft_sisters_colors = meristems_metadata[colnames(chimeric_ordered_correlation_matrix)[assigned_wt_sisters],"color"]
#browser()

for (j in 1:length(mutant_meristem_to_plot)){

if (svg) {
svglite(paste0(base_dir,myDate,"_example",as.character(j),"_",mutant_genotype_name,"sister_assignment_example_",mutant_meristem_to_plot[j],".svg"),width=w/ppi,height=h/ppi) } else {
png(paste0(base_dir,myDate,"_example",as.character(j),"_",mutant_genotype_name,"sister_assignment_example_",mutant_meristem_to_plot[j],".png"),width=w,height=h) }
par(mar=c(2,5,2,1),bty="L")
plot(chimeric_ordered_correlation_matrix[mutant_meristem_to_plot[j],],
ylim=c(-1,1),
col=meristems_metadata[colnames(chimeric_ordered_correlation_matrix),"color"],
xlab ="",
xaxt="n",
yaxt="n",
ylab = paste0(mutant_genotype_name, " meristem #",mutant_meristem_to_plot[j]),
pch=19,
xaxs="i",
yaxs="i",
cex.lab=1.5,
cex.axis=1.5,
cex=1.3
)
magaxis(side=2,at=c(-1,0,1),lables=c("-1","0","1"),las=2)
points(top_wts[mutant_meristem_to_plot[j],],chimeric_ordered_correlation_matrix[mutant_meristem_to_plot[j],top_wts[mutant_meristem_to_plot[j],]],col="black",pch=1,lty=3,lwd=2.2,cex=1.4)
points(assigned_wt_sisters[mutant_meristem_to_plot[j]],chimeric_ordered_correlation_matrix[mutant_meristem_to_plot[j],assigned_wt_sisters[mutant_meristem_to_plot[j]]],col="black",pch=1,lty=3,lwd=4.5,cex=1.4)
abline(h=0,
  lty=2,col="grey",lwd=0.8)

if (!is.null (wt_phase_assigment) & !is.null (wt_phase_colors)) {
  
	for (i in 1:(length(background_limits)-1)) {
			rect(ybottom = (-1),
			ytop = 1,
			xleft = background_limits[i],
			xright = background_limits[i+1],
			col=wt_phase_colors[i],
			border=NA)
	} 
}

dev.off()

}

}
####################################################################################################
####################################################################################################



# This function show the comparison between mutant pseudo-time and the order of their corresponding WT-sister along the WT trajectory.
# This is the basis for temporal alignment, and if it doesn't look good we might have a problem (technical, or biological) in aligning those two process.
####################################################################################################
smt_compare_pseudotime_ranks <- function (chimeric_ordered_correlation_matrix = NULL,
										  assigned_wt_sisters = NULL,
  										  wt_phase_assigment = NULL,
										  wt_phase_colors = NULL,
										  meristems_metadata = smer_md,
										  svg = FALSE,
										  ppi = 72,
										  base_dir = "./",
										  myDate = substring(Sys.time(),1,10),
										  mutant_genotype_name = "",
										  w = 600,
										  h = 600) {

										  
										  
										  
sister_colors = meristems_metadata[colnames(chimeric_ordered_correlation_matrix)[assigned_wt_sisters],"color"]

if (!is.null (wt_phase_assigment) & !is.null (wt_phase_colors)) {
background_limits = c(1,cumsum(table(wt_phase_assigment))) }

if (svg) { svglite(paste0(base_dir,myDate,"_",mutant_genotype_name,"_rank.svg"),width = w/ppi, height=h/ppi) } else {
		       png(paste0(base_dir,myDate,"_",mutant_genotype_name,"_rank.png"),width = w, height=h) }

par(mar=c(5,5,1,1),bty="L")
plot(assigned_wt_sisters,ylim=c(0,ncol(chimeric_ordered_correlation_matrix)),
		col=sister_colors,
        xlab = paste0("ordered mutant (",mutant_genotype_name,") meristems"),
        ylab = "rank of WT sister",
        pch=19,
        xaxs="i",
        yaxs="i",
        cex.lab=1.5,
        cex.axis=1.5,
        cex=1.5,
		yaxt="n",
		xaxt="n"
        )
magaxis(side=1,cex.axis=1.5)
magaxis(side=2,cex.axis=1.5)

abline(b=ncol(chimeric_ordered_correlation_matrix)/nrow(chimeric_ordered_correlation_matrix),
           a=0,col="grey",lwd=0.8)
		   
if (!is.null (wt_phase_assigment) & !is.null (wt_phase_colors)) {
for (i in 1 : (length(background_limits)-1))
rect(xleft = 1,
         xright = length(assigned_wt_sisters),
         ybottom = background_limits[i],
         ytop = background_limits[i+1],
         col=wt_phase_colors[i],
         border=NA)
}

dev.off()
									  
}
####################################################################################################
####################################################################################################										  




## this function compare trends of expression along aligned order 
## (for example - of mutant and its corresponding WT sisters expression)
####################################################################################################
smt_compare_trends_by_ord_season2 <- function(
							mat1 = NULL,
							mat2 = NULL,
							base_dir = "./",
							gene_to_plot,
							gene_nm,
							k_roll = 11,
							log2_exp = FALSE,
							log_reg = 0.1,
							trend_col1 = "black",
							trend_col2 = "red",
							ignore_meristems = NULL,
							my_xlab = "ordered meristems",
							ylab_col = "black",
							ylab_dist = 1.8,
							txt_size = 1.8,
							plot_margins = c(8,5,1,1),
							bg_group_colors = NULL,
							bg_group_fracs = NULL,
							remove_frame = FALSE,
							trend_lwd = 0.7,
							my_ylim = NULL,
							values_las = 2,
							xlab_dist = 3.5,
							my_height = 300,my_width=500,
							ppi = 72,
							svg = FALSE)
{
require(svglite)
if (is.null(mat1) | is.null(mat2)){ stop ("insert two datasets to compare trends (ordered_meristems parameter")}
		
if(!is.null(ignore_meristems)) { 
mat1 = mat1[,colnames(mat1)[!colnames(mat1) %in% ignore_meristems]]
mat2 = mat2[,colnames(mat2)[!colnames(mat2) %in% ignore_meristems]]
}

if (!svg) {
png(paste0(base_dir,gsub("/","_",gene_nm),"_sorted_meristems_trends.png"),width=my_width,height=my_height)
} else {
svg(paste0(base_dir,gsub("/","_",gene_nm),"_sorted_meristems_trends.svg"),width=my_width/ppi,height=my_height/ppi)}
par(mar=plot_margins)
ordered_exp1 = list(c(mat1[gene_to_plot,]),c(log2(mat1[gene_to_plot,]+log_reg)) )[[ifelse(log2_exp,2,1)]]
ordered_exp2 = list(c(mat2[gene_to_plot,]),c(log2(mat2[gene_to_plot,]+log_reg)) )[[ifelse(log2_exp,2,1)]]
ordered_exp1 = rollmean(ordered_exp1,k = k_roll)
ordered_exp2 = rollmean(ordered_exp2,k = k_roll)
concat_exp = c(ordered_exp1,ordered_exp2)
plot(ordered_exp1,
	type="l",
	lwd = trend_lwd,
	 col=trend_col1,
	 xlab="",
	 xaxs="i",
	 las = values_las,
	 frame.plot = !remove_frame,
	 ylab="",
	 ylim=list( c(min(concat_exp), max(concat_exp)),
				my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
	 cex.axis=txt_size,
	 cex.lab=txt_size,
	 xaxt="n",
	 yaxt="n")
lines(ordered_exp2,lwd = trend_lwd,col=trend_col2)
title(xlab=my_xlab,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=xlab_dist)
title(ylab=gene_nm,cex.lab=txt_size,cex.lab=txt_size,cex=txt_size,line=ylab_dist, col.lab = ylab_col)

if (!remove_frame){	 
axis(1,labels=FALSE,cex.axis=txt_size)}
yticks =seq(min(concat_exp),max(concat_exp),length.out=5)
yticks_round = round(yticks,digits=(-3))
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=(-2))}
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=(-1))}
if(sum(duplicated(yticks_round)>0)){
yticks_round = round(yticks,digits=0)}
magaxis(
		side=2,majorn=4,minorn="auto",
		frame.plot=!remove_frame,cex.axis=txt_size,cex=txt_size,las=values_las)

	if( (!is.null(bg_group_colors) & !is.null(bg_group_fracs))){
	bg_lims= cumsum(c(par("usr")[1],cumsum(par("usr")[2]-par("usr")[1])*bg_group_fracs))
			for (i in 1:length(bg_group_colors)){
			rect(bg_lims[i], par("usr")[3], bg_lims[i+1], par("usr")[4], col = bg_group_colors[i], border =NA)
			}
	}
	
dev.off()  						  
}						  					  
#######################################################################################################################
#######################################################################################################################





## This function limits maximal UMI per meristem
#######################################################################################################################
smt_squeeze_meristems <- function (mat = smer_umi,
								   min_umi = 1e4,
								   max_umi = 2e5)	{
return (cbind(mat[, colSums(mat) >= min_umi & colSums(mat)<max_umi], smt_ds(mat, max_umi)))								 
}
#######################################################################################################################
#######################################################################################################################






#######################################################################################################################
smt_filter_genes_by_two_groups <- function (mat = NULL,
											group1 = NULL,
											group2 = NULL,
											min_gene_umi = 25,
											min_umi_fraction_group2 = 0.9) {
umi1 = rowSums(mat[rowSums(mat)>= min_gene_umi, group1])
umi2 = rowSums(mat[rowSums(mat)>= min_gene_umi, group2])
return(rownames(mat)[rowSums(mat)>= min_gene_umi][(umi2 / (umi1 + umi2) ) > min_umi_fraction_group2])									
}
#######################################################################################################################
#######################################################################################################################





#######################################################################################################################
smt_plot_genes_cummulative <- function (mat = NULL,
 										group1 = NULL,
 										group2 = NULL,
										add_common_names = TRUE,
										base_dir = "./",
										plot_height = 250,
										plot_width = 500,
										show_legend = TRUE,
										file_name = "early_acting_transition_genes.png",
										heatmap_colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)[15:100],
										stage_genes_csv_path = "Supplementary_Tables/sup_table_stage_gene_modules.csv",
 										umi_fraction_threshold = 0.3,
										plot_svg = FALSE,
										font_size = 14) {
smt_load_common_nms()
mat_cum = t(apply(mat,1,cumsum))
mat_cum_n =  mat_cum / apply (mat_cum,1,max)
gene_order_mat = apply(mat_cum_n,1,function(x){which.max(x>umi_fraction_threshold)})
ordered_mat_cum_n = mat_cum_n[order(gene_order_mat),]

if (add_common_names) {
stage_genes_csv = read.csv(stage_genes_csv_path,row.names=1)
new_names = ifelse(is.na(stage_genes_csv[substring(rownames(ordered_mat_cum_n),1,14),"gene_name"]),
				 substring(rownames(ordered_mat_cum_n),1,14),
				 as.character(stage_genes_csv[substring(rownames(ordered_mat_cum_n),1,14),"gene_name"]))
				 
#add names of some of the early acting genes
new_names = gsub(substring(uf,1,14),"UF",new_names)
new_names = gsub(substring(fpf1,1,14),"FPF1",new_names)
new_names = gsub(substring(macrocalyx,1,14),"MC",new_names)
rownames(ordered_mat_cum_n) = new_names
}

ph_color = heatmap_colors

smt_mypheatmap(ordered_mat_cum_n,
				cluster_cols=FALSE,
				cluster_rows=FALSE,
				show_rownames=TRUE,
				show_colnames=FALSE,
				cellheight=12.5,
				fontsize=font_size,
				cellwidth=3,
				gaps_col=length(group1),
				my_color=ph_color,
				center0=FALSE,
				legend = show_legend,
				border = NA,
				annotation_legend = FALSE,
				h = plot_height,
				w = plot_width,
				outdir=base_dir,
				svg = plot_svg,
				fname = file_name)
}
#######################################################################################################################
#######################################################################################################################





#######################################################################################################################
smt_load_snpChip_data <- function (snp_csv_path = "Supplementary_Tables/sup_table_dst_SNPchip.csv",
								   chromosome_length_table_path = "Supplementary_Tables/SL2_4_chromosomeLength.csv",
								   bin_size_SNP_summation = 5e4,
								   myDate = substring(Sys.time(),1,10))

{
# As the SNP-chip array is based on this version, let's obtain tomato chromosome length according to version SL2.4 (from NCBI)
chrs = read.csv(chromosome_length_table_path)
snps = read.csv(snp_csv_path)
snps[,2] = as.numeric(eval(gsub(",","",as.character(snps[,2]))))

# filter the 7,720 markers to 7,370:
#
# filter non mapped markers
print(paste0("filtering out ", sum (snps$chrom=="NM"), " un-mapped markers"))
snps = snps[!snps$chrom=="NM",]
snps[,"chrom"] = factor(snps[,"chrom"],levels=c("CHR01","CHR02","CHR03","CHR04",
												"CHR05","CHR06","CHR07","CHR08",
												"CHR09","CHR10","CHR11","CHR12"))

chrs[,1] = c("CHR01","CHR02","CHR03","CHR04",
			 "CHR05","CHR06","CHR07","CHR08",
			 "CHR09","CHR10","CHR11","CHR12")
			 
# filter markers without calls in all 3 samples
colnames(snps) = c("chrom","pos","e137","modifier","M82")
print(paste0("filtering out ", sum (!grepl("A|C|G|T",snps$e137) | !grepl("A|C|G|T",snps$modifier) | !grepl("A|C|G|T",snps$M82)), " genomic markers without calls in all 3 samples"))
snps = snps[grepl("A|C|G|T",snps$e137) & grepl("A|C|G|T",snps$modifier) & grepl("A|C|G|T",snps$M82),]			 
##

#SNP calling
snps$e137_snp = ifelse(snps$e137==snps$M82,0,1)
snps$modifier_snp = ifelse(snps$modifier==snps$M82,0,1)

snps$chrom_id = as.numeric(substring(as.character(snps$chrom),4,5))
snps$chrom_id = as.integer(snps$chrom_id)
bg_map = data.frame(pos=c(),chrom_id=c())
print(paste0("Summing SNP counts per geomnic bins of ",bin_size_SNP_summation," base-pairs.."))
for (i in 1:nrow(chrs)){
temp_df = data.frame(pos=seq(1,chrs[i,2],bin_size_SNP_summation),chrom_id=i)
bg_map=rbind(temp_df,bg_map)
}

# compute n SNPS per bin
bin_map=bg_map
bin_map$e137_count = -1
bin_map$modifier_count = -1
for (i in 1:(nrow(bg_map)-1)){
relevant_markers = snps[(snps$chrom_id==bg_map[i,"chrom_id"]) & (snps$pos >= bg_map[i,"pos"]) & (snps$pos < bg_map[i+1,"pos"]),]
if (nrow(relevant_markers)==0) {bin_map[i,"e137_count"]=0;bin_map[i,"modifier_count"]=0;next;} else {
bin_map[i,"e137_count"] = sum(relevant_markers$e137_snp)
bin_map[i,"modifier_count"] = sum(relevant_markers$modifier_snp)}
}
bin_map = bin_map[-c(nrow(bin_map)),]

return(list(markers = snps, snp_counts = bin_map, background_counts = bg_map))}
#######################################################################################################################
#######################################################################################################################



#######################################################################################################################
smt_plot_snpChip_data <- function (snp_data = NULL,
								   background_data = NULL,
								   markers_data = NULL,
								   sample_name="",
   								   highlight_dst_modifiers_chrom4 = FALSE,
								   dst_modifiers_chrom4_path = "Supplementary_Tables/chrom4_dstGenesAnnotation_SL2_4_assembly.csv",
								   text_cnv=2,
								   base_dir = "./",
								   snp_color = "blue",
								   myDate = substring(Sys.time(),1,10),
								   factor_img=0.225,
								   fc_factor=0.025,
								   width_factor=2,
								   ppi = 72,
								   svg = FALSE)
{
# Optionally - we can point out at genes that are globally induced when DST gene is impaired, and that reside within a genomic island which is associated with dst fertility
if (highlight_dst_modifiers_chrom4) {
chr4_genes = read.csv("/net/mraid14/export/data/users/zoharme/tomato_Sc/workdir/chrom4_dstGenesAnnotation_SL2_4_assembly.csv") }

if (svg) { svglite(paste0(base_dir,myDate,"_tomatoChromes_chip_SNP_",sample_name,".svg"),width=1300*0.63*(12/23)*1.5/ppi,height=800*0.63/ppi) } else {
			   png(paste0(base_dir,myDate,"_tomatoChromes_chip_SNP_",sample_name,".png"),width=1300*0.63*(12/23)*1.5,height=800*0.63)

}
par(mar=c(4,6.5,3,1.5),mgp=c(3,0.8,0),bty="L")
plot(background_data$pos, 13-(background_data$chrom_id)+factor_img*fc_factor*0,
		cex=1*width_factor,
		xlab="",
		ylab="",
		yaxt="n",
		xaxt="n",
		cex.lab=text_cnv,
		cex.axis=text_cnv,
		col = "#DCDCDC02",
		pch=19)
		
#plot position of markers (no SNP where there are no markers..
points(markers_data[,"pos"],
	   (13-markers_data$chrom_id)+factor_img*fc_factor*0,
		cex=0.45*width_factor,col=alpha("grey53",0.3),pch="|")
		
#plot DST modifiers
if (highlight_dst_modifiers_chrom4) {		
points(chr4_genes[,"mid"],
	   (13-chr4_genes$chrom_ID)-factor_img*fc_factor*45,
		cex=0.2*width_factor,col="red",pch="|") }

		
points(snp_data[,"pos"],
	   (13-snp_data$chrom_id)+factor_img+fc_factor*snp_data[,3],
		cex= ifelse(snp_data[,3]==0,0.1*width_factor,0.3*width_factor),col=snp_color,pch=19)
axis(side=1,at=c(0,2e7,4e7,6e7,8e7),labels=c("0","20M","40M","60M","80M"),
	 cex=text_cnv,cex.axis=text_cnv,cex.lab=text_cnv)
axis(2,at=12:1,las=2,
labels=c(paste0("chr",1:12)),
cex.axis=text_cnv)
title(xlab="coordinate",line=2.7,cex.lab=text_cnv,cex.axis=text_cnv,cex=text_cnv,outer=FALSE)
legend(x="bottomright",
	   bty="n",
	   cex=text_cnv,
	   legend=sample_name,
	   text.col=c(snp_color),
	   x.intersp=0.2,
	   y.intersp=0.8)
dev.off()
}
#######################################################################################################################
#######################################################################################################################


# This function finds genes that are correlated with Number of leaves produced by SAMs
#######################################################################################################################
smt_find_age_genes <- function (veg_mat = smer_n,
								dstsft_mat = smer_n[,smer_md$genotype=="dst sft"],
								corr_method = "pearson",
								corr_thresh = 0.25,
								meristems_metadata = smer_md,
								min_gene_umi = 200) {
								
dstsft_leaf_corr = apply(dstsft_n[rowSums(smer_umi[,colnames(dstsft_n)])>=min_gene_umi,],
						 1,function(x){cor(x,meristems_metadata[colnames(dstsft_n),"Leaf_Num"],method=corr_method)})
vegetative_leaf_corr = apply(veg_n[rowSums(smer_umi[,colnames(veg_n)])>=min_gene_umi,],
						 1,function(x){cor(x,meristems_metadata[colnames(veg_n),"Leaf_Num"],method=corr_method)})
shared_genes = intersect(names(dstsft_leaf_corr), names(vegetative_leaf_corr))
dstsft_leaf_corr_f = dstsft_leaf_corr[shared_genes]
vegetative_leaf_corr_f = vegetative_leaf_corr[shared_genes]
age_genes = shared_genes[dstsft_leaf_corr[shared_genes] >= corr_thresh & vegetative_leaf_corr[shared_genes] >= corr_thresh]
age_df = data.frame(gene = age_genes,
					veg_corr = vegetative_leaf_corr[age_genes],
					dst_sft_corr = dstsft_leaf_corr[age_genes])
return(age_df)															
}
#######################################################################################################################
#######################################################################################################################




# This function just plot a scaled heatmap of normalized expression of selected genes (e.g, age-dependent ones) across genotype
#######################################################################################################################
smt_plot_age_genes <- function(mat = smer_n,
							   age_genes_to_plot = NULL,
							   meristems_metadata = smer_md,
							   scale_expression = TRUE,
							   add_common_names = TRUE,
							   base_dir = paste0("./",substring(Sys.time(),1,10),"_"),
							   file_name = "age_expression_heatmap.png",
							   show_legend = FALSE,
							   font_size = 14,
							   k_clusters_genes = 2,
							   gaps_width = 2,
							   wt_meristems_order = NULL,
							   sft_meristems_order = NULL,
							   dst_meristems_order = NULL,
							   dst_sft_order = colnames(dstsft_n)[order(smer_md[colnames(dstsft_n),"Leaf_Num"])],
							   plot_height=1200,
							   plot_width=2000,
							   plot_cellheight = 3.5,
							   plot_cellwidth = 2,
							   heatmap_colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
							   stage_genes_csv_path = "Supplementary_Tables/sup_table_stage_gene_modules.csv") {
				   
concat = cbind(mat [age_genes_to_plot,wt_meristems_order],
			   mat [age_genes_to_plot,sft_meristems_order],
			   mat [age_genes_to_plot,dst_meristems_order],
			   mat [age_genes_to_plot,dst_sft_order])
			   
if (scale_expression) { concat_n = t(apply(concat,1,scale)) }

if (add_common_names) {
	stage_genes_csv = read.csv(stage_genes_csv_path,row.names=1)
	new_names = ifelse(is.na(stage_genes_csv[substring(rownames(concat_n),1,14),"gene_name"]),
					"",
					as.character(stage_genes_csv[substring(rownames(concat_n),1,14),"gene_name"]))			 
	rownames(concat_n) = new_names
	}

meristem_counts = table(smer_md$genotype)[c("WT","sft","dst","dst sft")]

smt_mypheatmap(pmin(pmax(concat_n,(-3)),3),
		 center0=FALSE,
		 w=plot_width,h=plot_height,
		 cluster_cols = FALSE,
		 show_rownames = TRUE,
		 cellheight = 3.5,
		 cellwidth = 2,
		 fontsize = font_size,
		 border = NA,
		 legend = show_legend,
		 outdir = base_dir,
		 svg = FALSE,
		 my_color = heatmap_colors,
		 treeheight_row = 0,
		 fname = file_name,
		 gaps_col = rep( cumsum(meristem_counts)[1:(length(meristem_counts)-1)], gaps_width),
		 cutree_rows = k_clusters_genes)

}
#######################################################################################################################
#######################################################################################################################


### this function - 
### gets: a vector of solycs (possibly, fusion of more than one Solyc, separated by ";".
### returns: vector of same length, with additional annotations from arabidopsis-tomato table of 6,652 One-to-One Orthologues (derived from ENSEMBL database)
smt_findAnnotations_from_arabidopsis <- function (solycIDs,
												   source_annotation_csvFile = "Supplementary_Tables/Ensembl_6652_tomato_arabidopsis_oneToOne_orthologs.csv")								   
												   {
	anots = read.csv(source_annotation_csvFile)
	rownames(anots) = anots[,"Sl_ID"]
    names_per_gene = sapply(solycIDs,function(x){str_count(x,";")+1})
	new_names=rep("aaa",length(names_per_gene))
	for (i in 1:length(names_per_gene))
	{
		if (names_per_gene[i]==1)
		{
			if (!solycIDs[i] %in% as.character(anots[,"Sl_ID"])) 
			{new_names[i]=solycIDs[i]
			} else{
			new_names[i] = paste0(as.character(anots[ as.character(anots[,"Sl_ID"]) %in% solycIDs[i] ,"At_ID"]),"_", as.character(anots[ as.character(anots[,"Sl_ID"]) %in% solycIDs[i] ,"At_name"]))
			}
		} else if (names_per_gene[i]>1) {
			temp_name="";
			for (j in 1:names_per_gene[i]){
			temp_gene = unlist(strsplit(solycIDs[i],";"))[j]
				if (!temp_gene %in% as.character(anots[,"Sl_ID"])) {temp_gene_name=temp_gene
				} else {
				temp_gene_name = paste0(as.character(anots[ as.character(anots[,"Sl_ID"]) %in% temp_gene ,"At_ID"]),"_", as.character(anots[ as.character(anots[,"Sl_ID"]) %in% temp_gene ,"At_name"]))
				}
			if (j==1){
				temp_name=temp_gene_name
				}else{
				temp_name=paste0(temp_name,"__",temp_gene_name)
				}
			new_names[i] = temp_name
			}
		}
	}
return(new_names)
}
##############################################################################################
##############################################################################################






### these function makes a panel of gene expression in all meristems by their order
##############################################################################################
smt_plot_panel_by_order <- function(genes, label, 
									 trendline_color=NA,
									 trendline_k = 15,
									 trendline_cex = 2,
									 dot_size=3.5,unified_ylim=FALSE,
									 w = 2000,
									 h = 240,
									 base_dir = "./",
									 myDate = substring(Sys.time(),1,10),
									 meristems_metadata = smer_md,
									 wt_tree = NULL,
									 sft_tree = NULL,
									 dst_tree = NULL,
									 wt_ordered_meristems =  NULL,
									 sft_ordered_meristems = NULL,
									 dst_ordered_meristems = NULL,
									 ylab_distance=4.4,
									 show_genotype_title=TRUE,
									 panel_margins = c(5,5,4,1),
									 plots_margins = c(6,6,4,1),
									 text_size = 2.2,
									 log2_expression = FALSE,
									 svg=FALSE,ppi=72) {
require(svglite)

if (unified_ylim) {
if (length(genes)>1) {ylim_panel = range(colSums(smer_n[genes,])); if (log2_expression) {ylim_panel = log2(ylim_panel+1)} } else {
					  ylim_panel = range(smer_n[genes,]); if (log2_expression) {ylim_panel = log2(ylim_panel+1)} }
} else { ylim_panel=NULL }

# reorder dst;sft meristems
dstsfts = rownames(meristems_metadata)[meristems_metadata$genotype=="dst sft"]
dstsft_ord = dstsfts[order(meristems_metadata[dstsfts,"Leaf_Num"])]

bg_colors_sft = c("#00ff0035","#ffff0035","#ffA50035","#ff000035","#8A2BE235")
sft5 = cutree(sft_tree, k=5)
names(sft5) = sft_ordered_meristems
bg_fracs_sft = table(sft5)/length(sft5)

bg_colors_dst = c("#00ff0035","#ffff0035","#ffA50035","#ff000035","#8A2BE235")
dst5 = cutree(dst_tree, k=5)
names(dst5) = dst_ordered_meristems
bg_fracs_dst = table(dst5)/length(dst5)

bg_colors_wt = c("#00ff0035","#ffff0035","#ffA50035","#ff000035","#8A2BE235")
wt5 = cutree(wt_tree, k=5)
names(wt5) = wt_ordered_meristems
bg_fracs_wt = table(wt5)/length(wt5)

bg_colors_dstsft = c("#00ff0035")
bg_fracs_dstsft = 1

if (!svg) {
png(paste0(base_dir,myDate,"expresion_panel_",label,".png"),width=w,height=h)} else {
svg(paste0(base_dir,myDate,"expresion_panel_",label,".svg"),width=w/ppi,height=h/ppi)} 

par(mfcol=c(1,4), mar = panel_margins, mgp=c(3,0.8,0))

smt_plot_by_ord_season2(
	mat_n = smer_n[,names(wt5)],
	svg=FALSE,
    gene_to_plot = genes,
	log2_exp = log2_expression,
	txt_size = text_size,
	my_cex = dot_size,
	remove_frame=TRUE,
	my_width = 800*0.9,
	my_height = 250*0.85,
	mat = smer_umi[,names(wt5)],
	gene_nm = label,
	log_reg = 1,
	values_las = 2,
	bg_vec = NULL,
	plot_axs = "i",
	figure_main = ifelse(show_genotype_title,"WT",""),
	bg_group_colors = bg_colors_wt,
	bg_group_fracs = bg_fracs_wt,
	my_ylim = ylim_panel,
    base_dir = paste0(base_dir,myDate,"_",label,"_genes_wt_"),
	k_roll = trendline_k,
	ylab_dist = ylab_distance,
	xlab_dist=1.2,
	plot_margins = plots_margins,
	my_xlab = "",
	ylab_col = "black",
	trend_col = trendline_color,
	save_to_file=FALSE,
	k_cex = trendline_cex,
	ordered_meristems= names(wt5))

smt_plot_by_ord_season2(
	mat_n = smer_n[,names(sft5)],
	svg=FALSE,
    gene_to_plot = genes,
	log2_exp = log2_expression,
	txt_size = text_size,
	my_cex = dot_size,
	remove_frame=TRUE,
	my_width = 800*0.9,
	my_height = 250*0.85,
	mat = smer_umi[,names(sft5)],
	gene_nm = "",
	log_reg = 1,
	values_las = 2,
	bg_vec = NULL,
	figure_main = ifelse(show_genotype_title,expression(italic("sft")),""),
	plot_axs = "i",
	bg_group_colors = bg_colors_sft,
	bg_group_fracs = bg_fracs_sft,
	my_ylim = ylim_panel,
    base_dir = paste0(base_dir,myDate,"_",label,"_genes_sft_"),
	k_roll = trendline_k,
	ylab_dist=4.4,
	xlab_dist=1.2,
	plot_margins = plots_margins,
	my_xlab = "",
	ylab_col = "black",
	trend_col = trendline_color,
	k_cex = trendline_cex,
	save_to_file=FALSE,
	ordered_meristems= names(sft5))

smt_plot_by_ord_season2(
	mat_n = smer_n[,names(dst5)],
	svg=FALSE,
    gene_to_plot = genes,
	log2_exp = log2_expression,
	txt_size = text_size,
	my_cex = dot_size,
	remove_frame=TRUE,
	my_width = 800*0.9,
	my_height = 250*0.85,
	mat = smer_umi[,names(dst5)],
	gene_nm = "",
	figure_main = ifelse(show_genotype_title,expression(italic("dst")),""),
	log_reg = 1,
	values_las = 2,
	bg_vec = NULL,
	plot_axs = "i",
	bg_group_colors = bg_colors_dst,
	bg_group_fracs = bg_fracs_dst,
	my_ylim = ylim_panel,
    base_dir = paste0(base_dir,myDate,"_",label,"_genes_dst_"),
	k_roll = trendline_k,
	ylab_dist=4.4,
	xlab_dist=1.2,
	plot_margins = plots_margins,
	my_xlab = "",
	ylab_col = "black",
	trend_col = trendline_color,
	k_cex = trendline_cex,
	save_to_file=FALSE,
	ordered_meristems= names(dst5))	
	
smt_plot_by_ord_season2(
	mat_n = smer_n[,dstsft_ord],
	svg=FALSE,
    gene_to_plot = genes,
	log2_exp = log2_expression,
	txt_size = text_size,
	my_cex = dot_size,
	remove_frame=TRUE,
	my_width = 800*0.9,
	my_height = 250*0.85,
	mat = smer_umi[,names(cutree(uf_rows_oclust_reord,3))],
	gene_nm = "",
	log_reg = 1,
	values_las = 2,
	bg_vec = NULL,
	plot_axs = "i",
	bg_group_colors = bg_colors_dstsft,
	bg_group_fracs = bg_fracs_dstsft,
	my_ylim = ylim_panel,
    base_dir = paste0(base_dir,myDate,"_",label,"_genes_dstsft_"),
	k_roll = trendline_k,
	save_to_file=FALSE,
	ylab_dist=4.4,
	xlab_dist=1.2,
	figure_main = ifelse(show_genotype_title,expression(italic("dst;sft")),""),
	plot_margins = plots_margins,
	my_xlab = "",
	ylab_col = "black",
	trend_col = trendline_color,
	k_cex = trendline_cex,
	ordered_meristems= dstsft_ord)
dev.off()	
}
##############################################################################################
##############################################################################################



##############################################################################################
smt_compare_meristems <- function (group1 = NULL,
								   group2 = NULL,
								   mat_umi=smer_umi,
								   max_umi_per_meristem = 2e6,
								   min_umi_per_meristem = 1e5,
								   ds=FALSE)
{

if (ds) {
		mat = smer_ds(mat_umi,min_umi_per_meristem)
		group1 = group1[group1 %in% colnames(mat)]
		group2 = group2[group2 %in% colnames(mat)]
		} else {
		mat=smt_squeeze_meristems (mat_umi,
								   min_umi = min_umi_per_meristem,
								   max_umi = max_umi_per_meristem)}
							
tot1 = rowSums(mat[,group1])
tot2 = rowSums(mat[,group2])
N = min(sum(tot1), sum(tot2))

gene_tot1 <<- tot1*N/sum(tot1)
gene_tot2 <<- tot2*N/sum(tot2)
}
##############################################################################################
##############################################################################################


##############################################################################################
smt_plot_meristems_scatter <- function(tot1 = NULL, 
										tot2 = NULL,
										dots_cex = 0.5, 
										text_cex=2,
										base_dir = "./",
										myDate = substring(Sys.time(),1,10),
										n_reg = 5, 
										fig_nm = NULL,
										top_genes = 10,
										bad_genes = NULL,
										highlight_genes = NULL,
										highlight_genes_nms = NULL,
										highlight_genes_x = NULL,
										highlight_genes_y = NULL,
										highlight_genes_colors = "cyan",
										dots_color = "black",
										gene_characters_to_show = 28,
										highlight_cex_bonus = 0.5,
										fig_h=900, fig_w=1200,
										segments_calibration_right = 1.1,
										segments_calibration_left = 1.1,										
										figure_margins = c(25,25,25,25),
										lab1="grp1", lab2="grp2",
										main="compare bulk",
										plot_grid = FALSE,
										pt_col1 = "darkred", pt_col2 = "darkblue",
										show_gene_ticks = T)

{
   if(is.null(fig_nm)) { stop("insert fig_nm to save it to")}       
   if(!is.null(bad_genes)) { tot1 = tot1[!names(tot1) %in% bad_genes]
							 tot2 = tot2[!names(tot2) %in% bad_genes]
							}       
   
		png(paste0(base_dir,myDate,"_",fig_nm),w=fig_w,h=fig_h)
		par(mgp=c(3,0.7,0),mar = figure_margins, bty="L")
        lr = log2(n_reg+tot2)-log2(n_reg+tot1)

        xs = log2(n_reg+tot1)
        ys = log2(n_reg+tot2)

		max_value = max(c(xs,ys))
		min_value = min(c(xs,ys))
		
        plot(xs, ys, cex = dots_cex, pch=19, col = dots_color, 
			 xlim=c(min_value, max_value),
			 ylim=c(min_value, max_value),
			 cex.axis= text_cex,
			 cex.lab = text_cex,
             xlab=lab1,
			 ylab=lab2)
        top5 = names(head(sort(lr),n=top_genes))
        bottom5 = names(tail(sort(lr),top_genes))
        top5 = top5[order(ys[top5])]
        bottom5 = bottom5[order(ys[bottom5])]
        points(xs[bottom5], ys[bottom5], cex=dots_cex+0.5, pch=19, col=pt_col1)
        points(xs[top5], ys[top5], cex=dots_cex+0.5, pch=19, col=pt_col2)

        xmin = min(xs)
        xmax = max(xs)
        lab_xs = seq(xmin, xmax, length.out=length(highlight_genes))
        ymin = min(ys)
        ymax = max(ys)
        lab_ys = seq(ymin, ymax, length.out=top_genes)
		
		
        if(!is.null(highlight_genes)) {
                points(xs[highlight_genes], ys[highlight_genes],
                                                                cex=dots_cex+highlight_cex_bonus, pch=21, lwd=2,col = highlight_genes_colors)
				text(x=xs[highlight_genes] + highlight_genes_x,
					 y=ys[highlight_genes] + highlight_genes_y,
					 highlight_genes_nms,
					 cex=text_cex,
					 col=highlight_genes_colors)
                mtext(main, side=1)
        } else {
                mtext(main, side=3)
        }

        if(show_gene_ticks) {
                mtext(substring(top5,1,gene_characters_to_show), at=lab_ys, side=4, line=2, las=2, cex=text_cex)
                mtext(substring(bottom5,1,gene_characters_to_show), at=lab_ys, side=2, line=3, las=2, cex=text_cex)
                segments(x0 = min_value*segments_calibration_left, y0=lab_ys, x1=xs[bottom5], y1=ys[bottom5], xpd=T)
                segments(x0 = max_value*segments_calibration_right, y0=lab_ys, x1=xs[top5], y1=ys[top5], xpd=T)
        }
        if (plot_grid) {grid()}
                dev.off()

}
##############################################################################################
##############################################################################################	



## this function just load some common gene names into the workspace
##############################################################################################
smt_load_common_nms <- function() {
wox13 <<- "Solyc02g082670_YE_WOX_WOX13"
ham <<- "Solyc01g090950_YE_HAM"
wus <<- "Solyc02g083950_YE_WOX_WUS"
tfam2 <<- "Solyc05g055020_YE_TFAM_2"
tfam3 <<- "Solyc09g025280_YE_TFAM_3"
aspargine <<- "Solyc06g007180_YE_Fals_direct_target_asparagine_synthase_1"
asparaginase <<- "Solyc04g078460_YE_Florigen_asparaginase"
fa <<- "Solyc03g118160_YE_Falsiflora"
hdzip <<- "Solyc02g024070_YE_Fals_direct_target_Class_III_homeodomain_leucine_zipper"
agl66_a <<- "Solyc07g052700_YE_Novel"
agl66_b <<- "Solyc07g052707_YE_Novel__Solyc07g052730"
stylish <<- "Solyc02g084680_YE_stylish"
aux_reg <<- "Solyc04g010330_YE_AuxR_Auxin_regulated_protein"
Solyc07g040890 <<- "Solyc07g040890_SG_Gibberellin_receptor"
Solyc05g012380 <<- "Solyc05g012380_SG_Glucan_endo-1,3-beta-glucosidase-like"
Solyc04g054800 <<- "Solyc04g054800_SG_UPSTREAM_OF"
wrky_chr1 <<- "Solyc01g095100_SG_WRKY_transcription"
jointless2<<- "Solyc12g038510_YE_SEP_J2"
procera <<-"Solyc11g011260_YE_GA_PROCERA"
mazg<<- "Solyc09g097930_YE_MazG_(promoter_flowering)"
pny <<- "Solyc09g011380_YE_PNY"
sp <<- "Solyc06g074350_YE_CETS_SP"
ssp <<- "Solyc02g083520_YE_Fals_direct_target_SSP"
ap1 <<- "Solyc02g089210_YE_FUL/AP2_AP1"
jip56 <<-"Solyc00g179240_YE_SOC1_Jip56"
tkn1 <<- "Solyc04g077210_YE_TKN1"
wiry <<- "Solyc03g118770_YE_WOX_WOX1_wiry"
stm <<- "Solyc02g081120_YE_STM"
compound_infloresence <<- "Solyc02g077390_YE_WOX_WOX9_S"
yabby <<- "Solyc01g091010_YE_Fals_direct_target_YABBY_like_transcription_factor"
substilin <<- "Solyc01g006660_SG_Subtilisin-like_protease"
log7_gene <<- "Solyc10g082020_YE_CK_LOG_(LOGB)_(OX:LOG7)"
bop2 <<- "Solyc10g079460_YE_BOP2"
bop3 <<- "Solyc10g079750_YE_BOP3"
slySBP3 <<-"Solyc10g009080_YE_Fals_direct_target_SPL"
lob30 <<- "Solyc01g091420_YE_LOB_domain_30_CR_LBD30"
puchi <<- "Solyc01g067540_YE_ERF_12_PUCHI"
at_hook_chr1 <<- "Solyc01g080960_YE_AT_hook_CHR1"
at_hook_nsf1 <<-  "Solyc01g091470_YE_AT_HOOK_NSF_1"
at_hook_nsf37 <<- "Solyc08g077010_YE_AT_HOOK_NSF_37"
ck_crf <<- "Solyc08g081960_YE_CK_CRF"
pectinesterase <<- "Solyc02g080200_SG_Pectinesterase__Solyc02g080210"
macrocalyx <<- "Solyc05g056620_YE_FUL/AP1_MC__Solyc05g012020_YE_SEP_rin"
jointless2<<- "Solyc12g038510_YE_SEP_J2"
stylish3 <<- "Solyc11g064800_YE_stylish_CR_sty3"
tfam1 <<- "Solyc02g069510_YE_TFAM_1"
dst <<- "Solyc03g006900_YE_DST"
ap2c <<- "Solyc02g093150_YE_Fals_direct_target_APETALA2c"
lmi2 <<- "Solyc05g048830_YE_Fals_direct_target_MYB_TF_(LMI2)"
tmf <<- "Solyc09g090180_YE_TMF_(TFAM_founder)"
ful1 <<- "Solyc06g069430_YE_FUL/AP1_FUL1"
ful2 <<- "Solyc03g114830_YE_FUL/AP1_FUL2"
pin1 <<- "Solyc10g078370_YE_pin1_e9310"
bop1 <<- "Solyc04g040220_YE_BOP1"
arr_b <<- "Solyc01g065540_YE_CK_ARR_B"
stylish1 <<-"Solyc04g080970_YE_stylish_CR_sty1_"
wrky18 <<- "Solyc07g065260_YE_wrky18"
phantastica <<- "Solyc09g010840_YE_Fals_direct_target_phantastica"
spgb2 <<- "Solyc02g061990_YE_Gbox2_CR_SPGB2"
rossman <<- "Solyc01g091660_SG_NAD(P)-binding_Rossmann-fold"
svp04 <<- "Solyc04g076280_YE_AGL24/SVP_SVP04"
gata19 <<- "Solyc02g062750_SG_sugar_facilitator__Solyc02g062760"
sep3 <<- "Solyc05g015750_YE_SEP_TM5/SEP3"
agl6 <<- "Solyc01g093960_YE_AGL6" 
uf <<- "Solyc09g005070_YE_bHLH_UF"
fpf1 <<- "Solyc01g066970__Solyc01g066957_SG_Flowering_promoting__Solyc01g066980"
clv3 <<- "Solyc11g071380_YE_CLV3"
sft <<- "Solyc03g063100_YE_CETS_SFT"
bark <<- "Solyc08g066880_SG_Bark_storage"
fat1 <<- "Solyc02g020920_YE_FAT1"
spl15 <<- "Solyc10g078700_SG_Squamosa_promoter"
lax10 <<- "Solyc10g076790_YE_lax-10"
lypoxigenase3 <<-"Solyc03g122340_SG_lipoxygenase_D"
mlp550 <<- "Solyc09g014550_YE_Florigen_Major_latex_like/CR-MLP550"
ej2 <<- "Solyc03g114840_YE_SEP_ej2"
myb17 <<- "Solyc05g048830_YE_myb17"
tsjt1 <<-"Solyc01g079880_SG_Stem-specific_TSJT1"
tsjt1_like <<-"Solyc03g006490_SG_glutaminase_domain-containing"
tsjt1_like2 <<-"Solyc02g078500_SG_Stem-specific_protein"
anantha <<-  "Solyc02g081670_YE_Fals_partner_Anantha"
lin <<- "Solyc04g005320_YE_SEP_lin"
han <<- "Solyc02g062750_SG_sugar_facilitator__Solyc02g062760"



#batchy gene modules
mod_lateral_dissection_1 <<- c(
"Solyc00g009020_SG_MALE_GAMETOPHYTE",
"Solyc01g005810_SG_MAK16_protein-like",
"Solyc01g006450_SG_Enoyl_reductase",
"Solyc01g007060_SG_BTB/POZ_domain-containing",
"Solyc01g008820_SG_signal_peptide",
"Solyc01g022750_SG_NADH_dehydrogenase",
"Solyc01g058260_SG_Poly(A)_polymerase",
"Solyc01g059880_SG_ATP-citrate_synthase",
"Solyc01g059930_SG_Adenine_nucleotide",
"Solyc01g080750_SG_ARM_repeat",
"Solyc01g087120_SG_F1-ATP_synthase",
"Solyc01g091890_SG_Ubiquitin-related_modifier",
"Solyc01g094420_SG_Ulp1_protease",
"Solyc01g094470_SG_ketose-bisphosphate_aldolase",
"Solyc01g094550_SG_Acyl-CoA_thioesterase",
"Solyc01g094560_SG_60S_ribosomal",
"Solyc01g095410_SG_Eukaryotic_translation",
"Solyc01g096270_SG_Cytochrome_b5",
"Solyc01g097310_SG_Sec-independent_protein",
"Solyc01g097500_SG_LSTK-1-like_kinase",
"Solyc01g098340_SG_DnaJ_domain-containing",
"Solyc01g098550_SG_Indole-3-glycerol_phosphate",
"Solyc01g100320_SG_Disulfide-isomerase-like_protein",
"Solyc01g100530_SG_Protein_EFR3",
"Solyc01g101030_SG_FACT_complex",
"Solyc01g103020_SG_Survival_of",
"Solyc01g103250__Solyc01g103245_SG_Zinc_finger",
"Solyc01g104000_SG_Serine_hydroxymethyltransferase__Solyc01g104015_SG_Embryo-specific_3",
"Solyc01g106540_SG_neurofilament_heavy",
"Solyc01g106750_SG_Transmembrane_protein",
"Solyc01g109910_SG_RING/U-box_superfamily",
"Solyc01g110410_SG_ubiquitin_fusion",
"Solyc01g111520_SG_Calcium-dependent_lipid-binding",
"Solyc01g111560_SG_Rac-like_GTP-binding",
"Solyc01g111730_SG_Protein_phosphatase",
"Solyc02g014870_SG_Thioredoxin-like_protein",
"Solyc02g021440_SG_Non-specific_serine/threonine",
"Solyc02g050287_SG_SNF2_domain-containing",
"Solyc02g065300_SG_Leucyl-tRNA_synthetase",
"Solyc02g069680_SG_Vacuolar_protein",
"Solyc02g070510_SG_Proteasome_subunit",
"Solyc02g072130_SG_Protein_transport",
"Solyc02g077680_SG_Alpha-1,4_glucan",
"Solyc02g078360_SG_Glutaredoxin",
"Solyc02g079970_SG_bHLH_transcription",
"Solyc02g080400_SG_Abnormal_spindle-like",
"Solyc02g080430_SG_Mediator_of",
"Solyc02g081680_SG_Nucleolar_complex",
"Solyc02g081700_SG_Proteasome_subunit",
"Solyc02g082760_SG_ethylene-responsive_catalase",
"Solyc02g083890_SG_Chalcone-flavanone_isomerase",
"Solyc02g086880_SG_formate_dehydrogenase",
"Solyc02g089820_SG_proteasome_maturation",
"Solyc02g091970",
"Solyc02g093300_SG_DNA_polymerase",
"Solyc02g094160_SG_Ribonucleoside-diphosphate_reductase",
"Solyc03g005300_SG_ADP-ribosylation_factor",
"Solyc03g007440_SG_Plastid_division",
"Solyc03g007670_SG_defense_signal",
"Solyc03g025420_SG_SKP1_family",
"Solyc03g026040_SG_LRR_receptor-like",
"Solyc03g026340_SG_Calcium-dependent_protein",
"Solyc03g062680_SG_Serine/threonine-protein_phosphatase",
"Solyc03g063600_SG_Guanylate_kinase",
"Solyc03g082440_SG_Lung_seven",
"Solyc03g083530_SG_40S_ribosomal",
"Solyc03g093220_SG_transmembrane_protein",
"Solyc03g095600_SG_Kinase_interacting",
"Solyc03g096460_SG_wound/stress_protein",
"Solyc03g097970_SG_B-cell_receptor-associated-like",
"Solyc03g098050_SG_Calmodulin_6",
"Solyc03g098450_SG_Kinase_interacting",
"Solyc03g113030_SG_Aldose_1-epimerase",
"Solyc03g113405_SG_Plasma_membrane",
"Solyc03g113450_SG_Receptor-like_protein",
"Solyc03g113730_SG_B12D_protein",
"Solyc03g114070_SG_Rac-like_small",
"Solyc03g114960_SG_ABC_transporter",
"Solyc03g115820_SG_Ribulose-phosphate_3-epimerase",
"Solyc03g117030_SG_Polynucleotidyl_transferase",
"Solyc03g117040_SG_XH/XS_domain-containing",
"Solyc03g117730_YE_Fals_direct_target_TLP",
"Solyc03g117780",
"Solyc03g118410_SG_Acyl_carrier",
"Solyc03g118480_SG_Zinc-finger_domain",
"Solyc03g120280_SG_Ran-binding_protein",
"Solyc03g120495_SG_RNA-binding_(RRM/RBD/RNP__Solyc03g120500_SG_auxin-regulated_IAA27",
"Solyc03g120630_SG_ribosomal_protein",
"Solyc03g121310_SG_RWD_domain-containing",
"Solyc03g123530_SG_CCAAT/enhancer-binding_zeta",
"Solyc04g007370_SG_Chaperone_protein",
"Solyc04g008820_SG_High_mobility",
"Solyc04g011360_SG_GTP-binding_protein",
"Solyc04g011370_SG_Photosystem_I",
"Solyc04g011400_SG_UDP-glucuronate_decarboxylase",
"Solyc04g051370_SG_26S_proteasome",
"Solyc04g051830_SG_Developmentally-regulated_GTP-binding",
"Solyc04g055170_SG_annexin_p35",
"Solyc04g058090_SG_Phosphoribosylformylglycinamidine_cyclo-ligase",
"Solyc04g074240_SG_5-adenylylsulfate_reductase-like",
"Solyc04g074910_SG_40S_ribosomal",
"Solyc04g080160_SG_COP9_signalosome",
"Solyc04g080610_SG_ornithine_carbamoyltransferase",
"Solyc04g082260_SG_Kinase_family",
"Solyc04g082450__Solyc04g082460_SG_Catalase_3",
"Solyc04g082640_SG_Chloroplast_inner",
"Solyc05g005160_SG_ATP-citrate_synthase",
"Solyc05g005873_SG_Ribosomal_protein",
"Solyc05g010300_SG_DNA-directed_RNA",
"Solyc05g010420_SG_S-adenosylmethionine_decarboxylase",
"Solyc05g010423_SG_S-adenosylmethionine_decarboxylase",
"Solyc05g013920_SG_Cysteine_protease",
"Solyc05g023680_SG_RING-finger_E3",
"Solyc05g025600_SG_Photosystem_II",
"Solyc05g044480_SG_RNA_helicase",
"Solyc05g050200_SG_Eukaryotic_translation",
"Solyc05g052470_SG_Ferritin",
"Solyc05g052800_SG_60S_ribosomal",
"Solyc05g053070_SG_Dehydrin_family",
"Solyc05g054100_SG_Transmembrane_protein",
"Solyc05g054640_SG_2-oxoglutarate_dehydrogenase",
"Solyc05g054800_SG_LEM3_(Ligand-effect",
"Solyc05g055310_SG_Heavy_metal",
"Solyc05g055450_SG_RNA-binding_Nova-1",
"Solyc05g056310_SG_T-complex_protein",
"Solyc06g005940_SG_Protein_disulfide-isomerase",
"Solyc06g007800_SG_Mitochondrial_import",
"Solyc06g008010_SG_Mediator_of",
"Solyc06g034220_SG_Vesicle-associated_1-1-like",
"Solyc06g043150_SG_cell_cycle",
"Solyc06g050770_SG_Alpha-soluble_NSF",
"Solyc06g050980_SG_Ferritin",
"Solyc06g051810_SG_XH/XS_domain-containing",
"Solyc06g053460_SG_Prefoldin_subunit",
"Solyc06g053470_SG_Proteasome_assembly",
"Solyc06g053490_SG_FAM136A-like_protein",
"Solyc06g062760_SG_50S_ribosomal",
"Solyc06g062940_SG_RING-box_protein",
"Solyc06g064810_SG_Zinc_finger",
"Solyc06g069230_SG_DNA_mismatch",
"Solyc06g069790_SG_Gibberellin-regulated_protein",
"Solyc06g071310_SG_Pollen-specific_protein",
"Solyc06g071960_SG_Nucleoside_diphosphate",
"Solyc06g073330_SG_Lysine--tRNA_ligase",
"Solyc06g073430_SG_Ribosomal_protein",
"Solyc06g074670_SG_UDP-apiose/xylose_synthase",
"Solyc06g074730_SG_Argonaute_5",
"Solyc06g075580_SG_Kinesin-like_protein",
"Solyc06g082090_SG_Methionine_aminopeptidase,related",
"Solyc06g082580_SG_Eukaryotic_translation",
"Solyc06g082650_SG_60S_ribosomal",
"Solyc06g082870_SG_60S_ribosomal",
"Solyc06g082930_SG_Protein_FRIGIDA-like",
"Solyc06g083110_SG_YGGT_family",
"Solyc06g083530_SG_Vesicle-associated_membrane__Solyc06g083535_SG_multidrug_resistance-associated",
"Solyc06g083870_SG_structural_maintenance",
"Solyc06g083880",
"Solyc06g084230_SG_40S_ribosomal",
"Solyc06g084430_SG_Cupredoxin_superfamily",
"Solyc07g005560_SG_Eukaryotic_translation",
"Solyc07g006280_SG_Tetraspanin_family",
"Solyc07g026770_SG_Mitochondrial_ATP",
"Solyc07g026960_SG_Metallopeptidase_M24",
"Solyc07g039330_SG_Transducin_family",
"Solyc07g041490_SG_Photosystem_II",
"Solyc07g049480_SG_Cleavage_and",
"Solyc07g061930_SG_response_regulator",
"Solyc07g062310_SG_UNE1-like_protein",
"Solyc07g063960",
"Solyc07g064040_SG_bHLH_transcription",
"Solyc07g064680_SG_Autophagy-related_protein",
"Solyc07g065030_SG_Syntaxin-51",
"Solyc07g065310_SG_DUF21_domain-containing",
"Solyc07g066580_SG_Sulfurtransferase",
"Solyc08g008210_SG_vacuolar_proton",
"Solyc08g015860_SG_Rab_GDP",
"Solyc08g016510_SG_Proteasome_subunit",
"Solyc08g016670_SG_Calcyclin-binding_protein",
"Solyc08g041890_SG_Importin_subunit",
"Solyc08g062630_SG_aminopeptidase_M1",
"Solyc08g067030_SG_transmembrane_protein",
"Solyc08g075090_SG_bHLH_transcription",
"Solyc08g075160_SG_Bifunctional_purine",
"Solyc08g075950_SG_Growth-regulating_factor",
"Solyc08g076120_SG_MIP18_family",
"Solyc08g076160_SG_Guanine_nucleotide-binding",
"Solyc08g078390_SG_peroxisomal_acyl-CoA",
"Solyc08g080830_SG_Receptor_kinase",
"Solyc08g082740_SG_Signal_recognition",
"Solyc08g082850_SG_ABC_transporter",
"Solyc08g083300_SG_Ubiquitin_carboxyl-terminal",
"Solyc09g007540_SG_valyl-tRNA_synthetase",
"Solyc09g007920",
"Solyc09g008340_SG_transmembrane_protein",
"Solyc09g009210_SG_transcriptional_activator",
"Solyc09g009300_SG_Protein_BUD31",
"Solyc09g009640_SG_Small_nuclear",
"Solyc09g010210_SG_endo-1,4-beta-glucanase_precursor",
"Solyc09g011170_SG_Prf_interactor",
"Solyc09g011840_SG_Lectin_protein",
"Solyc09g018630_SG_Bis(5'-adenosyl)-triphosphatase",
"Solyc09g018750_SG_CBS_domain-containing",
"Solyc09g020130_SG_60S_ribosomal",
"Solyc09g042710_SG_weak_chloroplast",
"Solyc09g042750_SG_Acyl-CoA_thioesterase",
"Solyc09g074470_SG_Kinase_interacting",
"Solyc09g082320_SG_Proteasome_subunit",
"Solyc09g082520_SG_40S_ribosomal",
"Solyc09g083080_SG_MAR-binding_protein",
"Solyc09g090820_SG_phospholipase_D",
"Solyc09g091280_SG_Retinoblastoma-related_protein",
"Solyc09g092510_SG_Pyrroline-5-carboxylate_reductase",
"Solyc09g092550_SG_30S_ribosomal",
"Solyc10g006490_SG_Trafficking_protein",
"Solyc10g007400_SG_DNA_polymerase",
"Solyc10g054060_SG_Ribosomal_protein",
"Solyc10g054660_SG_Carbohydrate-binding-like_fold",
"Solyc10g055630_SG_plasma_membrane",
"Solyc10g074860_SG_Transmembrane_protein",
"Solyc10g076350_SG_Macrophage_migration",
"Solyc10g077120_SG_Photosystem_II",
"Solyc10g078300_SG_Single-stranded_nucleic",
"Solyc10g078325_SG_MSCS-like_2__Solyc10g078330",
"Solyc10g078660_SG_60S_ribosomal",
"Solyc10g079090_SG_Chaperone_protein",
"Solyc10g080110_SG_Cytoplasmic_tRNA",
"Solyc10g080390_SG_Kinase_family",
"Solyc10g080940_SG_Tubulin_beta",
"Solyc10g081310_SG_Histidine_triad",
"Solyc10g083600_SG_Brix_domain-containing",
"Solyc10g084160_SG_Small_RNA",
"Solyc10g086730_SG_Fructose-1,6-bisphosphatase",
"Solyc11g006660_SG_Eukaryotic_peptide",
"Solyc11g008270__Solyc11g008280_SG_Carboxypeptidase",
"Solyc11g008430_SG_Ras-related_small",
"Solyc11g009030_SG_glycine-rich_protein",
"Solyc11g012070_SG_Acyl-protein_thioesterase",
"Solyc11g012870_SG_lysM_domain-containing",
"Solyc11g012920_SG_translation_initiation",
"Solyc11g062300_SG_ARM_repeat__Solyc11g062310",
"Solyc11g068400_SG_Cytochrome_b-c1",
"Solyc11g068480_SG_Nuclear_factor",
"Solyc11g068490",
"Solyc11g068810_SG_Aspartyl/glutamyl-tRNA_(Asn/Gln)",
"Solyc11g071260",
"Solyc11g072260_SG_40S_ribosomal",
"Solyc11g072450_SG_ATP_synthase",
"Solyc11g072740_SG_Protein_BPS1",
"Solyc11g072880_SG_Calcium-transporting_ATPase",
"Solyc12g010060_SG_Eukaryotic_translation",
"Solyc12g011290_SG_Kinesin-like_protein",
"Solyc12g014400__Solyc12g014397_SG_Cell_differentiation",
"Solyc12g019750_SG_Polypyrimidine_tract-binding-like",
"Solyc12g038980_SG_H/ACA_ribonucleoprotein",
"Solyc12g040680_SG_mitogen-activated_protein",
"Solyc12g042950_SG_Plastidic_ATP/ADP-transporter",
"Solyc12g055800_SG_vacuolar_H+-ATPase",
"Solyc12g055830_SG_Soluble_inorganic",
"Solyc12g056100_SG_Ubiquitin_conjugating",
"Solyc12g056800_SG_Oxidoreductase_family",
"Solyc12g077590_SG_Peptidyl-tRNA_hydrolase",
"Solyc12g089030_SG_Ubiquitin-conjugating_enzyme",
"Solyc12g095910_SG_cysteine_protease",
"Solyc12g095940_SG_DUF1764_domain",
"Solyc12g096040_SG_PHD_finger",
"Solyc12g099660_SG_Glucosidase_2",
"Solyc12g099820_SG_Signal_recognition")

mod_lateral_dissection_2 <<- c(
"Solyc01g005620_SG_oxoglutarate/malate_translocator",
"Solyc01g068530_SG_Ribosomal_protein",
"Solyc01g080600_SG_Histone_H3",
"Solyc01g086970_SG_A20/AN1_zinc",
"Solyc01g088505_SG_Dynamin,_putative__Solyc01g088520",
"Solyc01g099670_SG_meloidogyne-induced_giant",
"Solyc01g099830_SG_60S_ribosomal",
"Solyc01g111280__Solyc01g111275_SG_cold_shock__Solyc01g111300",
"Solyc02g014310_SG_glycine_rich",
"Solyc03g005470_SG_HR-like_lesion-inducing",
"Solyc03g071620_SG_Histone_H2B",
"Solyc03g080160_SG_Nascent_polypeptide",
"Solyc04g063290_SG_40S_ribosomal",
"Solyc04g074300_SG_40S_ribosomal",
"Solyc05g005690_SG_Ribosomal_protein",
"Solyc05g051000_SG_40S_ribosomal",
"Solyc05g051290_SG_HMG-Y-related_A",
"Solyc05g053780_SG_RNA_binding",
"Solyc05g055440_SG_Histone_H2B",
"Solyc06g051310_SG_Clathrin_heavy",
"Solyc06g065475_SG_RNA-binding_family",
"Solyc06g069120_SG_translocase_subunit",
"Solyc06g071720_SG_60S_ribosomal",
"Solyc06g071870_SG_Ribosomal_protein",
"Solyc06g073790_SG_40s_ribosomal",
"Solyc06g075800_SG_Histone_H2B",
"Solyc07g053690_SG_Nucleic_acid-binding",
"Solyc08g080720_SG_Selenoprotein_H",
"Solyc09g007350",
"Solyc09g062970_SG_GDSL-like_Lipase/Acylhydrolase",
"Solyc09g075160_SG_60S_ribosomal",
"Solyc09g092710_SG_Glycine-rich_protein",
"Solyc09g092720__Solyc09g092715_SG_RGG_repeats",
"Solyc10g051390_SG_RNA-binding_glycine-rich",
"Solyc10g054560_SG_V-type_proton",
"Solyc10g078540_SG_H/ACA_ribonucleoprotein",
"Solyc10g078620_SG_Ribosomal_protein",
"Solyc10g084310_SG_Ribosomal_protein",
"Solyc10g086020_SG_interactor_of",
"Solyc11g007920",
"Solyc11g007930_SG_Histone_H2B",
"Solyc11g012130_SG_Blue_copper",
"Solyc11g072030_SG_Lipid_transfer",
"Solyc12g056540_SG_Histone_H3",
"Solyc12g096540_SG_40S_ribosomal",
"Solyc12g100060_SG_Zinc_finger")


mod_apoptosis <<- c(
"Solyc00g006830_SG_no_exine",
"Solyc00g013140",
"Solyc00g013155_SG_Cytochrome_c",
"Solyc00g013180_SG_NADH-ubiquinone_oxidoreductase",
"Solyc00g014830_SG_NADH_dehydrogenase",
"Solyc00g019750_SG_transmembrane_protein",
"Solyc00g021630_SG_NADH-ubiquinone_oxidoreductase",
"Solyc00g025290_SG_Alpha-1,4_glucan",
"Solyc00g171810_SG_ATP_synthase",
"Solyc01g021640_SG_Katanin_p80",
"Solyc01g058500_SG_TBP-associated_factor",
"Solyc01g108440_SG_calmodulin-binding_protein",
"Solyc02g083200_SG_F-box_and",
"Solyc03g063070_SG_vacuolar_protein",
"Solyc03g063480_SG_high_chlorophyll",
"Solyc05g025580_SG_N-glycosylase/DNA_lyase",
"Solyc06g011490_SG_Unknown_protein",
"Solyc06g024207_SG_alpha-1,2-Mannosidase",
"Solyc07g039295_SG_glutamate_decarboxylase",
"Solyc08g036520_SG_Diphosphomevalonate_decarboxylase",
"Solyc08g066990_SG_Potassium_channel",
"Solyc09g010880_SG_Rhomboid-like_protein",
"Solyc09g050020_SG_cytochrome_b",
"Solyc10g045760_SG_NADH-ubiquinone_oxidoreductase",
"Solyc10g049590_SG_DNA-directed_RNA",
"Solyc11g028160_SG_Cytochrome_c",
"Solyc11g056250_SG_glutamine-dependent_asparagine",
"Solyc11g056320_SG_Photosystem_II",
"Solyc11g056360_SG_Ribosomal_protein__Solyc11g056370",
"Solyc11g063500_SG_nitrite_reductase",
"Solyc11g063510_SG_Ypt/Rab-GAP_domain",
"Solyc11g063520_SG_sequence-specific_DNA",
"Solyc12g036350_SG_nudix_hydrolase",
"Solyc12g036415_SG_SAUR-like_auxin-responsive")



wt_age_mrks <<- c(
"Solyc01g008150_SG_cytochrome_B5-like",
"Solyc01g014180_SG_A20/AN1_zinc",
"Solyc01g066880_SG_Heavy_metal",
"Solyc01g080070_SG_Heavy_metal",
"Solyc01g081050_SG_D-alanine--poly(phosphoribitol)_ligase",
"Solyc01g091590_SG_BON1-associated_protein",
"Solyc01g094870_SG_PLAC8_family",
"Solyc01g094970_SG_Pectin_lyase-like",
"Solyc01g095530_SG_cytoplasmic_tRNA",
"Solyc01g104740_SG_Multiprotein-bridging_factor",
"Solyc01g106940_SG_Myb_transcription__Solyc01g106950__Solyc01g106943_SG_Myosin_heavy",
"Solyc01g107490_SG_basic_helix-loop-helix",
"Solyc02g063520_SG_Homeobox-leucine_zipper",
"Solyc02g069180_SG_SBP_(S-ribonuclease",
"Solyc02g077970_SG_Late_embryogenesis",
"Solyc02g080610_SG_CDP-diacylglycerol--glycerol-3-phosphate_3-phosphatidyltransferase__Solyc02g080615_SG_Ankyrin_repeat",
"Solyc02g084900",
"Solyc03g006920_SG_Kinase_family",
"Solyc03g082390_SG_RNA-binding_protein__Solyc03g082395_SG_RNA-binding_protein__Solyc03g082400",
"Solyc03g083290_SG_40S_ribosomal",
"Solyc03g083350_SG_phosphatidylinositol_4-kinase",
"Solyc03g095300_SG_Calcium-dependent_lipid-binding",
"Solyc03g097380_SG_Heavy_metal",
"Solyc03g113460_SG_NAD(P)-binding_Rossmann-fold",
"Solyc03g118160_YE_Falsiflora",
"Solyc03g119010_SG_Pre-mRNA-splicing_factor",
"Solyc04g010330_YE_AuxR_Auxin_regulated_protein",
"Solyc04g011545_SG_receptor_like__Solyc04g011550",
"Solyc04g018110_SG_Hop-interacting_protein",
"Solyc04g072110_SG_Eukaryotic_translation",
"Solyc04g078070_SG_Zinc_finger",
"Solyc04g079930_SG_Histone_deacetylase",
"Solyc05g011980_SG_Photosystem_II",
"Solyc05g015500_SG_Ran_BP2/NZF",
"Solyc05g053010_SG_Lectin_receptor",
"Solyc05g055160_SG_DnaJ-like_protein",
"Solyc05g055650_SG_Quinohemoprotein_ethanol",
"Solyc06g005530_SG_flowering_promoting",
"Solyc06g010170_SG_lysine_ketoglutarate",
"Solyc06g011370_SG_Chaperone_protein",
"Solyc06g036170",
"Solyc06g049020_SG_Sorghum_bicolor",
"Solyc06g069430_YE_FUL/AP1_FUL1",
"Solyc06g070980_SG_Ubiquitin-conjugating_enzyme__Solyc06g070985_SG_Tubulin_beta-1",
"Solyc06g072040_SG_CONSTANS_interacting",
"Solyc07g006380_SG_Defensin-like_protein",
"Solyc07g040890_SG_Gibberellin_receptor",
"Solyc07g052700_YE_Novel",
"Solyc07g052707_YE_Novel__Solyc07g052730",
"Solyc07g053140_SG_Zinc_finger",
"Solyc07g054850_SG_transmembrane_protein",
"Solyc07g063180_SG_Dynein_light",
"Solyc07g064130_SG_Polyubiquitin",
"Solyc07g066050_SG_RNA_polymerase",
"Solyc07g066260_SG_Protein_phosphatase",
"Solyc08g005090_SG_DUF740_family",
"Solyc08g006320_SG_WRKY_transcription",
"Solyc08g008350_SG_ribosome_maturation",
"Solyc08g062960",
"Solyc08g065160_SG_Mediator_of",
"Solyc08g066880_SG_Bark_storage",
"Solyc08g079700_SG_Zinc_finger",
"Solyc09g005670_SG_Type_I__Solyc09g005690",
"Solyc09g010640_SG_HSP20-like_chaperones",
"Solyc09g014280_SG_HXXXD-type_acyl-transferase",
"Solyc09g014475_SG_Histone-lysine_N-methyltransferase",
"Solyc09g030450_SG_Receptor-like_kinase",
"Solyc09g063140_SG_Protein_LOW",
"Solyc09g064660_SG_Small_nuclear",
"Solyc10g006480_SG_Polyubiquitin",
"Solyc10g054910_SG_Peptidyl-prolyl_cis-trans",
"Solyc10g075070_SG_Non-specific_lipid-transfer",
"Solyc10g078700_SG_Squamosa_promoter",
"Solyc11g005640_SG_Polyubiquitin",
"Solyc11g006460_SG_DnaJ-like_protein",
"Solyc11g068770_SG_Polynucleotidyl_transferase",
"Solyc11g069850_SG_Telomere_repeat-binding",
"Solyc12g010310_SG_Zinc_finger",
"Solyc12g016190_SG_High_mobility",
"Solyc12g044640_SG_Mediator_of",
"Solyc12g095920_SG_BRISC_and",
"Solyc12g096350_SG_WRKY_transcription",
"Solyc12g098850_SG_ethylene-forming_enzyme",
"Solyc03g114830_YE_FUL/AP1_FUL2",
"Solyc03g006900_YE_DST",
"Solyc05g015510_YE_Fals_direct_target_Squamosa_promoter_binding_protein_10",
"Solyc05g015840_SG_Squamosa_promoter__Solyc05g015843_SG_Ribulose_bisphosphate")
}
##############################################################################################
##############################################################################################