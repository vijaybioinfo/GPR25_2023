
#!/usr/bin/R

######################
# Figures compendium #
######################

# ---
# Author: Francisco Emmanuel Castaneda-Castro
# Date: 2023-01-29
# ---

### ================== Figures Lung and Liver day 30 ================== ###

# Locking the clustering objects (lung)
# lung = pct25_pc20_0.4


# source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/lung08/scripts/figures/object_lock_signature.R")
# source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1_day30/scripts/figures/object_lock_signature.R")


source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1/scripts/figures/global_signature.R") # change to directory /home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/gpr25/exp1_day30/
source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')

resources = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  # filters_columns is.file.finished show_commas (cluster_reports)
  "/home/ciro/scripts/clustering/R/plotting.R", # cluster_reports
  "/home/ciro/scripts/clustering/R/utilities.R", # get_top_n
  "/home/ciro/scripts/handy_functions/devel/filters.R", # sample_even
  "/home/ciro/scripts/handy_functions/devel/plots.R", # plot_pct getlegend mytheme
  "/home/ciro/scripts/handy_functions/R/stats_summary_table.R" #stats_summary_table
)
for(i in resources){ source(i) }

{ cat(redb("### Secondary global variables ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
    # As you progress they appear in many tasks and so become more or less global
    source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
 
    gpr25_liver_day30@meta.data = joindf(gpr25_liver_day30@meta.data,
      as.data.frame(gpr25_liver_day30@reductions$umap@cell.embeddings))

    gpr25_exp3@meta.data = joindf(gpr25_exp3@meta.data,
      as.data.frame(gpr25_exp3@reductions$umap@cell.embeddings))

}

###########
{ cat(red("### Final figures ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
  setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/final_figures")
  dir.create("liver_day30")
  setwd("liver_day30")

#### LIVER DAY30 
  setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/final_figures")
  dir.create("liver_day30_corrected")
  setwd("liver_day30_corrected")

  { cat(redb("### Figure 4a. Blanks volcano ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))

    dir.create("volcano")
    results<- read.csv("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/gpr25/exp1_day30/results/dgea/WTvsKO_liver_day30/WTvsKO/JuPa03_Mo_KO_CD8_3HvsJuPa03_Mo_WT_CD8_3H/mastlog2cpm_results.csv")

      group1 = "JuPa03_Mo_KO_CD8_3H"
      group2 = "JuPa03_Mo_WT_CD8_3H"
      padjthr = 0.05
      fcthr = 0.25
      cols_filt = "minExp"
      output = "./"
      verbose = TRUE
      return_report = FALSE
      ##############
      pseuc = 1
      dtype = "CST"
      # volcano parameters
      showgenes = NULL
      # heatmap parameters
      cat("\n%%%%%%%%%%%%%%%%%%%%%% DGEA report %%%%%%%%%%%%%%%%%%%%%\n")
      means = NULL
      if(is.null(cols_filt)) cols_filt = "pattern123"
      the_report = list()
      if(is.null(names(group1))) names(group1) <- group1
      if(is.null(names(group2))) names(group2) <- group2

      if(verbose) cat("---------------------- Filtering results ---------------\n")
      resdf <- data.frame(
        results[which(results$padj <= 1), ],
        stringsAsFactors = FALSE, check.names = FALSE
      ) # order should ideally be by padj all the way through
      resdf <- resdf[order(resdf$padj), ]
      if(!"gene" %in% colnames(resdf)) resdf$genes = rownames(resdf)
      resdf[, "gene"] <- features_parse_ensembl(resdf[, "gene"])
      rownames(resdf) <- features_parse_ensembl(resdf[, "gene"])
      tvar <- list(
        which(resdf[, "log2FoldChange"] <= -fcthr),
        which(resdf[, "log2FoldChange"] >= fcthr)
      )
      if(!"group" %in% colnames(resdf)){
        if(verbose) cat("Addding 'group' column\n")
        resdf$group = NA; resdf$group[tvar[[1]]] <- group1
        resdf$group[tvar[[2]]] <- group2
      }else{
        tvar <- list(which(resdf$group == group1), which(resdf$group == group2))
      }
      group_m <- sapply(c(names(group1), names(group2)), function(x){
        grep(pattern = paste0("^", x, "_mean"),x = colnames(resdf), value = TRUE)
      })
      if(length(group_m) == 2){
        if(verbose){
          cat("Using means as colour\n"); str(tvar, vec.len = 10, no.list = 4)
          str(group_m, vec.len = 10, no.list = 4)
        }
        resdf$Mean <- 0; means = "Mean"
        resdf$Mean[tvar[[1]]] <- round(log2(resdf[tvar[[1]], group_m[1]] + pseuc), 1)
        resdf$Mean[tvar[[2]]] <- round(log2(resdf[tvar[[2]], group_m[2]] + pseuc), 1)
      }
      genes2plot <- mysignames <- getDEGenes(
        resdf, pv = padjthr, fc = fcthr,
        gene_name = "gene", further = NULL, verbose = verbose
      )
      if(any(grep(cols_filt, colnames(resdf)))){
        tvar <- grep(cols_filt, colnames(resdf), value = TRUE)
        genes2plot <- genes2plot[genes2plot %in% resdf$gene[resdf[, tvar]]]
        if(verbose){ cat("Filtered by", tvar, "\n"); str(genes2plot) }
      }
      resdf$degs <- "Not_significant"
      resdf$degs[resdf[, "gene"] %in% genes2plot] <- "DEG"
      if(length(genes2plot) == 0)
        genes2plot<-getDEGenes(resdf,pv=0.2,fc=fcthr,gene_name="gene",v=TRUE)
      if(length(genes2plot) == 0)
        genes2plot<-resdf[bordering(resdf,cnames="log2FoldChange",n=50),"gene"]
      if(length(group_m) == 2) resdf$Mean[!resdf[, "gene"] %in% genes2plot] <- NA
      resdf$group[!resdf[, "gene"] %in% mysignames] <- NA
      if("pct_diff" %in% colnames(resdf))
        resdf$pct_diff <- abs(resdf$pct_diff)
        resdf[!resdf[, "gene"] %in% genes2plot, ]$pct_diff <- 0
      if(is.null(showgenes)){
        showgenes <- unique(c(
          bordering(resdf[genes2plot, ], cnames = "log2FoldChange", n = 10),
          head(resdf[genes2plot, "gene"], 10)))
      }

      cat("---------------------- Volcano -------------------------\n")

      source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1/scripts/volplot.R") # volplot
      # source("/home/ciro/scripts/handy_functions/devel/utilities.R") # getDEGenes
      source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
      source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_breaks
      source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes

      # source("/home/ciro/scripts/handy_functions/devel/volcano.R")

      # source("/mnt/BioAdHoc/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/asthma_biopsy/redo_analysis/figures/ciro/volplot.R")
      #ephitelials
      showgenes<-c("Apoe", "Cxcr4", "Tcf7", "Il7r", "Plac8", "Psmd4", "Cxcr6", "Gzmb", "Prf1", "S1pr1", "Ccr7", "Lef1", "Fasl", "Gzmk")

    pdf("volcano/DAY30_WTvsKO_blank_v2.pdf") #Change name for CD8
      volplot(
        resdf,
        pvalth = padjthr,
        lfcth = fcthr,
        pvaltype = "padj",
        lfctype = "log2FoldChange",
        col_feature = means,
        size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
        gene_name = "gene",
        group = "degs",
        check_genes = list(text = features_parse_ensembl(showgenes)),
        return_plot = TRUE,
        clipp = 4,
        xlims=c(-2,2),
        ylims=c(0,20),
        verbose = verbose
      ) + labs(
        size = "Delta %", color = paste0("Mean (", dtype, ")"),
        title = paste(group2, "(-) vs ", group1, "(+)")
      )
    dev.off()


    pdf("volcano/DAY30_WTvsKO_v2.pdf") #Change name for CD8
      volplot_wnames(
        resdf,
        pvalth = padjthr,
        lfcth = fcthr,
        pvaltype = "padj",
        lfctype = "log2FoldChange",
        col_feature = means,
        size_feature = if("pct_diff" %in% colnames(resdf)) "pct_diff",
        gene_name = "gene",
        group = "degs",
        check_genes = list(text = features_parse_ensembl(showgenes)),
        return_plot = TRUE,
        clipp = 4,
        xlims=c(-2,2),
        ylims=c(0,20),
        verbose = verbose
      ) + labs(
        size = "Delta %", color = paste0("Mean (", dtype, ")"),
        title = paste(group2, "(-) vs ", group1, "(+)")
      )
    dev.off()


    datavis <- GetAssayData(gpr25_exp1_day30)
    annot <- gpr25_exp1_day30@meta.data
    datatype <- "SeuratNormalized"
    verbose <- TRUE

    aver<-stats_summary_table(
      mat = datavis,
      groups = make_list(x = annot, colname = 'origlib', grouping = TRUE),
      rnames = rownames(datavis),
      datatype = datatype,
      verbose = verbose
    )

    write.csv(aver, "./volcano/mean_and_percetage_per_gene.csv")
  }

  { cat (redb("### Figure 4b. and Extended Data Figure 4a violins ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
      dir.create("violins") 

      gpr25_liver_day30$condition<-ifelse(grepl("KO",gpr25_liver_day30$origlib), "KO", "WT")
      table(gpr25_liver_day30$origlib, gpr25_liver_day30$condition)
      # gpr25_exp3$condition_cluster

      fconfigs = list(
            list(result_id = "violins/gpr25_exp1_liver_day30_corrected",
              edata = "gpr25_liver_day30@assays$RNA@data",
              metadata = "gpr25_liver_day30@meta.data",
              axis_x = list(col = 'condition', order=c("WT", "KO")),
              features = c("Klrg1", "Gzmk", "Tcf7", "Lef1", "Ccr7", "Prf1", "Gzmb", "Gzmk", "Cxcr3", "Cd226", "Cd160", "Ccr7", "Gzmk", "Nkg7", "Ctsw", "Fasl", "Cxcr4", "Cxcr6", "Ly6e", "Xcl1", "Cd28")
            )
          )


          source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1/scripts/plots.R")
          source("/home/ciro/scripts/figease/figease.R")
          source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1/scripts/source.R")
          pp_violins = fig_plot_violins(
            fconfigs, verbose = 1,
            theme_extra = labs(x = "Corrected plot", y = "Seurat Normalized"),
            colour_by = "pct", couls = couls_opt$red_gradient$white,
            box_median = TRUE, dots_mean = TRUE, chilli=FALSE)

    }

  saveRDS(gpr25_liver_day30@meta.data, "metadata_gpr25_liver_day30.rds")
  ######

  ###### Lung Day 30
  setwd("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/final_figures")
  dir.create("lung_day30")
  setwd("lung_day30")

  { cat(redb("### Figure 5d UMAP WT and KO and pie  ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
    dir.create("dim_reduction_features")

    mycols<-c('0'="#FF0000", '1'="#356BB4")
    #downsampling KO cells to have an equal number of cells between KO and WT, just for figure purposes
    WT_cells <-  which(gpr25_exp3$orig.lib == 'JuPa08_Mo_WT_CD8_1H')
    KO_cells <- which(gpr25_exp3$orig.lib == 'JuPa08_Mo_KO_CD8_1H')
    downsampled_KO_cells <- sample(KO_cells, 706)
    gpr25_exp3_WT_KO_downsampled <- gpr25_exp3[,c(WT_cells, downsampled_KO_cells)]


    fconfigs = list(
      list(result_id = "./dim_reduction_features/WTvsKO_umap",
        edata = "gpr25_exp3@assays$RNA@data", metadata = "gpr25_exp3@meta.data",
        axis_x = redu[[1]][1], axis_y = redu[[1]][2],
        vars2col = 'RNA_snn_res.0.4', col = mycols, facet = "orig.lib"
      ),
      list(result_id = "./dim_reduction_features/den_WTvsKO_umap_downsampled",
        edata = "gpr25_exp3_WT_KO_downsampled@assays$RNA@data", metadata = "gpr25_exp3_WT_KO_downsampled@meta.data",
        axis_x = redu[[1]][1], axis_y = redu[[1]][2],
        vars2col = 'RNA_snn_res.0.4', col = mycols, facet = "orig.lib"
      )
    )

    pp_clones = fig_plot_base(
      fconfigs[2], return_plot = TRUE, verbose = TRUE,
      theme_extra = function(x){
        x + geom_point(size=3, alpha = 0.8, stroke = 0.2) +  theme(
        plot.margin = unit(c(3, 3, 3, 8), "cm"))
      }
    )

    ##pie

      KO <- data.frame(
        group = c(1, 0),
        value = c(0.62, 0.38)
      )

      WT <- data.frame(
        group = c(1, 0),
        value = c(0.38, 0.62)
      )

    pie_wt<-
    ggplot(WT, aes(x="", y=value, fill=as.character(group)))+
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + scale_fill_manual(values=rev(c("#FF0000", "#356BB4")))

    pie_ko<-ggplot(KO, aes(x="", y=value, fill=as.character(group)))+
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) + scale_fill_manual(values=rev(c("#FF0000", "#356BB4")))

    pdf("./dim_reduction_features/pie_WT_lung_day30.pdf")
    print(pie_wt)
    dev.off()

    pdf("./dim_reduction_features/pie_ko_lung_day30.pdf")
    print(pie_ko)
    dev.off()


  }

  { cat(redb("### Figure 5e. Dotplot ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
    dir.create("dotplots")
    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/exp1/scripts/source.R")
    source("/mnt/bioadhoc-temp/Groups/vd-vijay/fcastaneda/rnaseq-sc-standar/gpr25/lung08/scripts/plots_dotplot.R")

    gpr25_exp3$RNA_snn_res.0.4_order_paper<-ifelse(gpr25_exp3$RNA_snn_res.0.4 == "0", "1", "0")

    fconfigs = list(
      list(result_id = "./dotplots/", sufix = "NewVersion_lung_gpr25_exp3_RNA_snn_res.0.4_16Jan2023",
      object = "gpr25_exp3",
        edata = "gpr25_exp3@assays$RNA@data", metadata = "gpr25_exp3@meta.data",
        size = c(5,4.6),
        axis_x = list(
          col = "RNA_snn_res.0.4_order_paper"),
        features = c("Tcf7", "Il7r", "Cd69", "Xcl1", "Gpr183", "Zeb2", "S1pr1", "S1pr5", "Gzma", "Gzmb", "Prf1", "Cx3cr1", "Fasl", "Klrg1") ,
        col=c('#f0e76e', '#ff0000'),
        scale_mean=TRUE,
        size_scale=10,
        return_plot=TRUE
      )
    )

    pp_curtains = fig_plot_curtain(fconfigs[1], verbose = 2)

    pp_curtains = fig_plot_curtain(fconfigs[1], verbose = 2)

  { cat(red("### Figure 5f and Extended Data Figure 5e Markers: UMAPs expression genes ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
      dir.create("dim_reduction_features")
      fconfigs = list(
        list(result_id = "dim_reduction_features/vPURPLE_sc_gpr25_exp3_", sufix = "/",
          edata = "gpr25_exp3@assays$RNA@data",
          metadata = "gpr25_exp3@meta.data", sample_even = TRUE,
          axis_x = redu[[1]][1], axis_y = redu[[1]][2],
          features = c("Ccr7", "Cd27", "Ctla2a", "Cxcr3", "Gpr183", "Il7r", "Ly6e", "Xcl1", "Cx3cr1", "Cxcr6", "Fasl", "Gzma", "Gzmb", "Itga1", "Itgal", "Klrg1", "Prf1", "S1pr5", "Zeb2" ),
          col = c("#d5d4d6", "#450aab"),
          plot_main = TRUE
          ))

          pp_markers = fig_plot_scatters(fconfigs, return_plot = TRUE, verbose = TRUE, theme_extra= function(x){
            x + geom_point(aes(size=1) +
              geom_point(aes(size = 0.6), shape = 1, color = "black", alpha = 0.1, stroke = 0.4))
          } )
      }


      fconfigs = list(
        list(result_id = "dim_reduction_features/vORANGE_sc_gpr25_exp3_", sufix = "_color_changed2_2/",
      edata = "gpr25_exp3@assays$RNA@data",
      metadata = "gpr25_exp3@meta.data", sample_even = TRUE,
      axis_x = redu[[1]][1], axis_y = redu[[1]][2],
      features = c("Ccr7", "Cd27", "Ctla2a", "Cxcr3", "Gpr183", "Il7r", "Ly6e", "Xcl1", "Cx3cr1", "Cxcr6", "Fasl", "Gzma", "Gzmb", "Itga1", "Itgal"," Klrg1", "Prf1", "S1pr5", "Zeb2"),
      plot_main = TRUE
      ))

      pp_markers = fig_plot_scatters(fconfigs, return_plot = TRUE, verbose = TRUE, theme_extra= function(x){
          x + geom_point(aes(size=1) +
          geom_point(aes(size = 0.6), shape = 1, color = "black", alpha = 0.1, stroke = 0.4)) +
          scale_colour_gradientn(colours =  c("#fff8a6","#f5e74e","#f0dd0c","#FFA500","#FF3500","#670000"), values=c(0,0.1,0.2,0.3,0.5,1))
        } )

    saveRDS(gpr25_exp3@meta.data, "metadata_gpr25_lung_day30.rds")
