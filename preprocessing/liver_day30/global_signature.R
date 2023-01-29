#!/usr/bin/R

###############
# Global data #
###############

# ---
# Author: Ciro Ramirez-Suastegui
# Date: 2021-02-22
# ---

# This script loads all the files used in a global way

if(requireNamespace("crayon", quietly = TRUE)){
  cyan = crayon::cyan; redb = crayon::red$bold; green = crayon::green
}else{ cyan = redb = green = c }

{ cat(redb("------------------ Setting global parameters -----------------\n"))
  cat(cyan("------------------ File names and initial objects\n"))
  fig_dir = "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/gpr25/"
  dir.create(fig_dir, showWarnings = FALSE); setwd(fig_dir)
  cat("Working at (fig_dir):", fig_dir, "\n")
  redu = list(umap = c("UMAP_1", "UMAP_2"), tsne = c("tSNE_1", "tSNE_2"))
  gpr25_exp1_clust = "RNA_snn_res.0.3"; gpr25_exp3_clust = "RNA_snn_res.0.8"
  global_objects_f = c(
    # gcolours = "/home/ciro/pbtumor/info/colours.csv",
    # gpr25_exp1 = "/mnt/beegfs/fcastaneda/gpr25_liver_xdoublets_clean1_clustersCD8n_cellcycleREGRESSION_object_lock_mean0.01_pct25_pc20.rds",
    gpr25_exp3 = "/mnt/beegfs/fcastaneda/gpr25_lung_08_xdoublets_clean1_clusters_noCD8n_cellcycleREGRESSION_5_object_lock_mean0.01_pct25_pc20.rds",
    gpr25_exp1_day30 = "/mnt/beegfs/fcastaneda/DAY30_gpr25_liver_xdoublets_clean1_clustersCD8n_cellcycleREGRESSION_object_lock_mean0.01_pct25_pc20.rds"
  )
  if(!exists("include")) include = names(global_objects_f)


  cat(cyan("------------------ Packages and functions\n")) # -------------------
  packages_funcs = c(
    "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
    "/home/ciro/scripts/handy_functions/devel/filters.R",
    "/home/ciro/scripts/handy_functions/devel/utilities.R",
    "/home/ciro/scripts/handy_functions/devel/plots.R",
    "/home/ciro/scripts/figease/source.R", # fig_global_objects
    "ggplot2", "cowplot", "patchwork", "Seurat", "dplyr"
  )
  loaded <- lapply(X = packages_funcs, FUN = function(x){
    cat("*", x, "\n")
    if(!file.exists(x)){
      suppressMessages(require(package = x, quietly = TRUE, character.only = TRUE))
    }else{ source(x) }
  }); theme_set(theme_cowplot())


  cat(cyan("------------------ Objects\n")) # ----------------------------------
  object_names = fig_global_objects(
    global_objects_files = global_objects_f,
    object_names = ls(pattern = "_f$|^redu|_clust$|_ident$"),
    reading_fun = 'readfile', stringsAsFactors = FALSE, check.name = FALSE, row.names = 1
  )
  source("/home/ciro/scripts/clustering/R/utilities.R")
  }

  # cat(green("------------------ Identities\n")) # ------------------------------
  #
  #   celltype_subset = c(
  #     "0"  =       "TAct",
  #     "1"  =       "TACTLow",
  #     "2"  =       "TIFNR_eff_reg",
  #     "3"  =       "TH1poly",
  #     "4"  =    "TFH",
  #     "5"  =        "TREG",
  #     "6"  = "TRegzma",
  #     "7"  =       "TH2",
  #     "8"  =          "TH17like",
  #     "9"  =          "CTXhi",
  #     "10"  =          "CTXlow",
  #     "11" =          "TH17",
  #     "12" =       "THDR",
  #     "13" =     "TFH10",
  #     "14" =    "c14",
  #     "15" =        "c15"
  #   )
  #   fgal_ident = list(
  #     celltype = gsub("~.*", "", celltype_subset),
  #     celltype_subset = gsub("Tcell~", "", celltype_subset),
  #     order = names(celltype_subset)
  #   )
  #   fgal = ident_set_names(
  #     object = fgal,
  #     ident = fgal_ident,
  #     cluster_column = fgal_clust
  #   )
  #   fgal_ident$colours = c(
  #     TAct = c("#feffba"),
  #     TACTLow = c("#e1e5eb"),
  #     TIFNR_eff_reg = c("#00FF00"),
  #     TH1_poly = c("#000dff"),
  #     TFH = c("#38d9f2"),
  #     TREG = c("#009900"),
  #     TRegzma = c("#96FFB4"),
  #     TH2 = c("#ff0000"),
  #     TH17like = c("#CB7BF4"),
  #     CTXhi = c("#A30032"),
  #     CTXlow = c("#ff33c9"),
  #     TH17 = c("#AF33FF"),
  #     THDR = c("#FFBD33"),
  #     TFH10 = c("#20999F"),
  #     c14 = c("#34211D"),
  #     c15= c("#52342E")
  #   )
  #   fgal_ident <- ident_colours(fgal_ident, mdata = fgal@meta.data)

  print(format(
    object_names,
    justify = "centre", quote = FALSE)
  ); rm(include)
  system("ls -loh")
