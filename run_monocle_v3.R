# user parameters
#------------------------------------
finrds <- "/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/test_monocle3/K_B6_NSDS5_Mins_Stromal.rds"
earlyMarkerGene <- "CD38"
#------------------------------------



# library
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(reticulate)
library(monocle3)
#library(tidyverse)
library(SeuratWrappers)
library(ggridges)



# check requirement packages to install monocle3
if(F){
  #monocle3 requirement:
  library('BiocGenerics')
  library('DelayedArray')
  library('DelayedMatrixStats')
  library('limma')
  library('lme4')
  library('S4Vectors')
  library('SingleCellExperiment')
  library('SummarizedExperiment')
  library('batchelor')
  library('HDF5Array')
  library('terra') #
  library('ggrastr')
  devtools::install_github('cole-trapnell-lab/monocle3')
}



# set src dir
source("/storage/chen/data_share_folder/22_10x_hs_AnteriorSegment_data/scAtacQC/R_functions_scATAC.R")


# get parameters
if(F){
  args = commandArgs(trailingOnly=TRUE)
  print(args)
  
  if (length(args)> 1) {
    stop("please type 1 parameters: finrds(r object)", call.=FALSE)
  } else {
    # default output file
    finrds   <- args[1] # input folder names
    cat("finrds_r.object: ",finrds,"\n")
  }
}

# get dir/file names
bName  <- basename(finrds)
dirIn  <- dirname(finrds)
prefix <- tools::file_path_sans_ext(bName)

# create out_dir
setwd(dirIn)
dirOut <- paste0("out_monocle3_",prefix)
dirOut <- setDir(dirIn, dirOut)

# load data
setwd(dirIn)
objA <- readRDS(finrds)
head(objA); objA

# plot
p1 <- DimPlot(objA, label = T)
#p2 <- DimPlot(data, split.by="sample", label = F,ncol=2)

setwd(dirOut)
png(paste0("01_umap_",prefix,".png"), width=300, height = 300, units="px")
print(p1)
dev.off()


# clustering
#cellgrp <- data.frame(table(Idents(data),data$sample))
#write.csv(cellgrp,"summary_table_cellcounts_bysample_clusters.csv")


# Converting seuratobject to celldataset object for Monocle3
cds <- as.cell_data_set(objA)
# Get cell metadata
head(colData(cds))
# Get feature/gene metadata
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))

#Get counts
#head(counts(cds))

# Retrieve clustering information from Surat object
# Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
summary(recreate.partitions)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions


#Assign cluster information
list.cluster <- objA@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

#Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- objA@reductions$umap@cell.embeddings
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")


png(paste0("02_1_umap_cluster_",prefix,"_before_trajectory.png"), width=500, height = 400, units="px")
cluster.before.traj
dev.off()


#Learn Trajectory
cds <- learn_graph(cds, use_partition = F)

png(paste0("02_2_umap_monocle3_",prefix,"_byCluster.png"), width=500, height = 400, units="px")
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,group_label_size = 5)
dev.off()

# check celltypes 
table(colData(cds)$celltype)


# order cells
#earlyMarkerGene <- "RRM2"
max.avp <- which.max(unlist(FetchData(objA, earlyMarkerGene)))
max.avp <- colnames(objA)[max.avp]

cds <- order_cells(cds, root_cells = max.avp)
png(paste0("03_1_monocle3_pseudotime_",prefix,".png"),width = 400,height = 300)
  plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
             label_branch_points = FALSE)
dev.off()

# save monocle objects

save_monocle_objects(cds=cds, directory_path=paste0('monocle3_cds_objects_',prefix), comment=paste0('This is cds. Stored on', Sys.Date()))




