
#######################################################
#  Setup
#######################################################
rm(list = ls())
set.seed(1)

saveName <- "PreprocessedSeurat_v6_250110"

annotationDir  <- "~/MultiomeProject/reference/Ppatens_v6/genes/"
tenXOutputPath <- "~/MultiomeProject/runCellRanger-arc_Count_2024-02-01/"
SeuratObjPaths <- "~/MultiomeProject/R_Analysis/SeuratObjects/"
figuresPath    <- "~/MultiomeProject/R_Analysis/Figures/"

# Load essential packages

pacman::p_load(Seurat,Signac,tidyverse, BiocGenerics, SeuratDisk, harmony, BSgenome.Ppatens.v6, JASPAR2020, JASPAR2024, TFBSTools, patchwork, motifmatchr)


#
#pacman::p_load('tidyverse','Seurat','SeuratDisk','SeuratData','Signac','SeuratObject','BSgenome.Ppatens.v6','rtracklayer','hdf5r','ggplot2','AnnotationFilter','harmony','BiocGenerics','GenomeInfoDb','GenomicFeatures','GenomicRanges','IRanges','Rsamtools','S4Vectors','TFBSTools','ggbio','motifmatchr','AnnotationDbi','R.utils','qlcMatrix','devtools','data.table','presto','remotes','gridExtra')

#######################################################
#  Save Progress TryCatch <~ head
#######################################################

tryCatch({

#######################################################
#  Load Annotation
#######################################################
# Check if the GTF file for Physcomitrium patens (ppsn) exists

if(!file.exists(paste0(annotationDir,"genes.gtf"))){
# If not, unzip the gzipped GTF file
gunzip(paste0(annotationDir,"genes.gtf.gz"), remove = FALSE)
}

# Import the GTF file
gtf<-rtracklayer::import(paste0(annotationDir,"genes.gtf"))

# Prepare gene coordinates
gene.coords <- gtf
gene.coords$gene_biotype <- "protein_coding"
gene.coords$gene_id <- str_replace_all(gene.coords$gene_id, "_", "-")
gene.coords$gene_name <- str_replace_all(gene.coords$gene_id, "_", "-")
gene.coords$tx_id <- gene.coords$transcript_id
gene.coords <- setNames(gene.coords, gene.coords$gene_id)
genome(gene.coords) <- "Assembly v6"

# Extract coding sequences (CDS) and untranslated regions (UTRs)
cds <- gene.coords[gene.coords$type == 'CDS']
cds$type <- "cds"
utr <- gene.coords[gene.coords$type %>% str_detect("utr")]
utr$type <- "utr"

# Create an annotation object combining genes, exons, CDS, and UTRs
annotation <- c(gene.coords[gene.coords$type == 'gene'], gene.coords[gene.coords$type == 'exon'], cds, utr)

#Clean up workspace
rm(annotationDir, gtf, gene.coords, cds, utr)

#######################################################
#  Create Seurat
#######################################################

sample_list = c("Wt_T00", "Wt_T06", "Wt_T12", "Wt_T24", "ste_T00", "ste_T06", "ste_T12", "ste_T24")#, "Wt_T00_2021")
inputdata <- Read10X_h5
ppsn <- CreateSeuratObject

# Read in 10x Genomics data for each sample
for (sample in sample_list) {
  inputdata <-
    Read10X_h5(
      paste0(
        tenXOutputPath,
        sample,
        "/outs/filtered_feature_bc_matrix.h5"
      )
    )
  
  # extract RNA data
  rna_counts <- inputdata$`Gene Expression`
  
  
  # Create Seurat object (P. patens single nuclei)

  options(Seurat.object.assay.version = 'v4')
  sample_seurat <-
    CreateSeuratObject(counts = rna_counts,
                       project = sample,
                       min.cells = 0)
  sample_seurat[["percent.mt"]] <-
    PercentageFeatureSet(sample_seurat, pattern = "^MT-")
  
  # Add the ATAC-seq data
  atac_counts <- inputdata$Peaks
  
  # we'll only use peaks in standard chromosomes 
  # standardChromosomes() function works on human genome,
  # so we generate the standard P. patensV6 chromosomes with sprintf())

  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% sprintf("Chr%02d", 1:26)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  
  
  frag.file <- paste0(tenXOutputPath,sample,"/outs/atac_fragments.tsv.gz")
  
  sample_seurat[["ATAC"]] <- CreateChromatinAssay(
     counts = atac_counts,
     sep = c(":", "-"),
     fragments = frag.file,
     min.cells = 0,
     annotation = annotation
   )
  
  assign(paste0(sample, "_seurat"), sample_seurat)
  
  rm(
    rna_counts,
    atac_counts,
    grange.counts,
    grange.use,
    frag.file,
    inputdata,
    sample_seurat
  )
}

#saveRDS(as.vector(mget(ls(pattern = "_seurat"))), file= paste0(SeuratObjPaths,"separate_samples_seurat.RDS"))

# merge individual samples into a single seurat object
ppsn <-
  merge(
    Wt_T00_seurat,
    y = c(
      Wt_T06_seurat,
      Wt_T12_seurat,
      Wt_T24_seurat,
      ste_T00_seurat,
      ste_T06_seurat,
      ste_T12_seurat,
      ste_T24_seurat   #,Wt_T00_2021_seurat
    ),
    add.cell.ids = sample_list,
    project = "ppsn"
  )

rm(list = c(ls(pattern = "_seurat"),'sample','sample_list','tenXOutputPath'))
ls()

#######################################################
#  Metadata Part I
#######################################################

ppsn@meta.data <-  ppsn@meta.data %>%

## Unabbreviated original identity
  mutate(Sample = orig.ident %>%
           str_replace_all(setNames(
             c("Wild type\n","\u394stemin\n","T0","T3;4.5;6","T10;12;14","T24;36;prot."),
             c("Wt_","ste_","T00","T06","T12","T24")))
         ) %>%
  {
    .$Sample <- factor(.$Sample, 
                       levels = paste0(c(rep("Wild type\n",4), rep("\u394stemin\n",4)),
                                       rep(c("T0", "T3;4.5;6","T10;12;14", "T24;36;prot."), 2))) 
    ;.} %>%

## Genotype
  mutate(Genotype = Sample %>% str_split("\n", simplify = T) %>% .[, 1]) %>%
  {.$Genotype <- factor(.$Genotype, levels = c("Wild type", "\u394stemin")) ;.} %>%

## Time shorthand
  mutate(Time = orig.ident %>% str_split("_", simplify = T) %>% .[, 2] %>% 
           str_remove("T") %>% as.numeric()) %>%

## Unabbreviated time
  mutate(`Time II` = Sample %>% str_split("\n", simplify = T) %>% .[, 2]) %>%
           {.$`Time II` <- factor(.$`Time II`, levels = c("T0","T3;4.5;6","T10;12;14","T24;36;protonema")) ;.}

#######################################################
#  TSS Enrichment
#######################################################

DefaultAssay(ppsn) <- "ATAC"
ppsn <- NucleosomeSignal(ppsn, n = ncol(ppsn) * 5000, verbose = TRUE)
ppsn <- TSSEnrichment(ppsn, fast = F)
#ggplot2::ggsave(paste0(figuresPath,"TSSPlot.png"),TSSPlot(ppsn))

save.image(file = paste0(SeuratObjPaths, saveName, "_preFilter.RData"))

#######################################################
#  Filtering data based on count:
#  low count or high count cells are discarded
#######################################################

plotVln <- function() {
  return(
    VlnPlot(
      object  = ppsn,
      ncol = 4,
      features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "TSS.enrichment"),
      log = TRUE
    ) + theme(legend.position = 'right') &
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
))}

ggplot2::ggsave(paste0(figuresPath, "unfilteredPlot.png"),plotVln(), width = 20)

ppsn <- subset(
  x = ppsn,
  subset = nCount_ATAC > 5000 &
    nCount_RNA < 20000 &
    nCount_RNA > 1000 &
    TSS.enrichment > 1
)
ppsn@misc$subsetInfo <-
  "nCount_ATAC > 5000 & nCount_RNA < 20000 & nCount_RNA > 1000 & TSS.enrichment > 1"

ggplot2::ggsave(paste0(figuresPath, "filteredPlot.png"),plotVln(), width = 20)

rm(plotVln)
#######################################################
#  Peak Calling
#######################################################
DefaultAssay(ppsn) <- "ATAC"
peaks <- CallPeaks(ppsn, macs2.path ="~/miniconda3/envs/snR-env/bin/macs3", combine.peaks = FALSE)

#discard peaks too close to edge of scaffold

peaks <- peaks[as.data.frame(ranges(peaks))$start >= 2, ]

# quantify counts in each peak

macs2_counts <- FeatureMatrix(
  fragments = Fragments(ppsn),
  features = peaks,
  cells = colnames(ppsn)
)

ppsn[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,
                                        fragments = Fragments(ppsn),
                                        annotation = annotation)


rm(peaks, macs2_counts)

#######################################################
#  Normalise Gene expression
#######################################################

DefaultAssay(ppsn) <- "RNA"
ppsn <- SCTransform(ppsn, verbose = FALSE)
ppsn <-  RunPCA(ppsn)


#######################################################
#  DNA accessibility data processing
#######################################################

DefaultAssay(ppsn) <- "peaks"
ppsn <- RunTFIDF(ppsn)
ppsn <-
  FindTopFeatures(ppsn, min.cutoff = 'q0') #q0 ~ just cut off outliers basically
ppsn <-
  RunSVD(ppsn) #Run singular value decomposition. Outputs a latent semantic indexing (lsi)

#SaveH5Seurat(ppsn, filename = paste0(SeuratObjPaths,"ppsn_preHarmony.h5Seurat"), overwrite = T)


#######################################################
#  Remove Batch Effects
#######################################################

# For gene expression

ppsn <-
  RunHarmony(
  ppsn,
  group.by.vars = "orig.ident",
  assay.use = "SCT"
)
ppsn@reductions$sct_harmony <- ppsn@reductions$harmony
ppsn[['harmony']] <- NULL

ppsn <- FindNeighbors(ppsn,reduction = 'sct_harmony',
dims = 1:11,
return.neighbor = T)%>%
RunUMAP(n.neighbors = 20,
min.dist = 0.01,
nn.name = 'SCT.nn',
reduction.name = 'gex_umap',
reduction.key = 'gexUMAP_'
)

# For chromatin accessibility

ggplot2::ggsave(paste0(figuresPath,"ElbowPlot.png"), ElbowPlot(ppsn, reduction = "lsi"))

ppsn <- RunHarmony(
  ppsn,
  group.by.vars = "orig.ident",
  assay.use = "peaks",
  reduction = "lsi",
  project.dim = FALSE,
  reduction.save = "lsi_harmony"
)

ggplot2::ggsave(paste0(figuresPath,"DepthCor.png"),DepthCor(ppsn, reduction = "lsi_harmony"))

ppsn <- FindNeighbors(
 object = ppsn,
 reduction = 'lsi_harmony',
 dims = 2:13,
 return.neighbor = T) %>%
RunUMAP(
 n.neighbors = 20,
 min.dist = 0.01,
 nn.name = 'peaks.nn',
 reduction.name = 'atac_umap')

#######################################################
# UMAP Visualisations
#######################################################

# build a joint neighbor graph using both assays
ppsn <- FindMultiModalNeighbors(
  object = ppsn,
  reduction.list = list("sct_harmony", "lsi_harmony"),
  modality.weight.name = c("SCT.weight", "peaks.weight"),
  dims.list = list(1:11, 2:13)
)

# build a joint UMAP visualization
ppsn <- RunUMAP(
  object = ppsn,
  n.neighbors = 20,
  min.dist = 0.01,
  nn.name = "weighted.nn",
  reduction.name = "umap",
  reduction.key = "wnnUMAP_"
)

#######################################################
# Clustering
#######################################################

ppsn <- FindClusters(
  ppsn,
  graph.name = "wsnn",
  #"weighted shared nearest-neighbour
  algorithm = 3,
  resolution = 0.4,
  verbose = FALSE
)
ppsn$wsnn_res.0.4<- as.numeric(as.character(ppsn$wsnn_res.0.4))%>%factor(0:max(.))

#save.image(file = paste0(SeuratObjPaths, saveName, ".RData"))
print("save1 done")

#######################################################
# Differentially Expressed Genes
#######################################################

DefaultAssay(ppsn) <- "SCT"
Idents(ppsn) <- "wsnn_res.0.4"

all_da_genes <- FindAllMarkers(
 ppsn,
 recorrect_umi=FALSE,
 assay = "SCT",
 logfc.threshold = 0,
 test.use = "wilcox",
 min.pct = 0,
 return.thresh = 1)

#save.image(file = paste0(saveName, ".RData"))
print("save2 done")

#######################################################
#  Metadata Part II
#######################################################
cluster_amount = levels(ppsn$wsnn_res.0.4) %>% length
genotype_amount = levels(ppsn$Genotype) %>% length

ppsn@meta.data <-  ppsn@meta.data %>%
  mutate(Clusters = wsnn_res.0.4) %>%
  { .$wsnn_res.0.4 <- NULL;.} %>%
  
## Reprogramming Cluster
  mutate(Reprogramming = ifelse(Clusters %in% c(1), TRUE, FALSE)) %>%

## Cell clusters divided by genotype
  mutate(`Cluster and Genotype` = paste0(Clusters,". ", Genotype)) %>%
  {
    .$`Cluster and Genotype` <- factor(.$`Cluster and Genotype`,
                                   levels = paste0(rep(x = levels(.$Clusters), 
                                                       times = rep(genotype_amount, times = cluster_amount)),
                                                   ". ", rep(levels(.$Genotype), cluster_amount))
    ) ;.} %>%
## Reprogramming, protonema apical or other, by genotype
  mutate(`Reprogramming by Genotype` = case_when(
    .$Clusters %in% c(1) & .$Genotype == "Wild type" ~ "Cluster 1,  Wild type",
    .$Clusters %in% c(1) & .$Genotype == "Δstemin" ~ "Cluster 1, Δstemin",
    .$Clusters %in% c(7) & .$Genotype == "Wild type" ~ "Cluster 7, Wild type",
    .$Clusters %in% c(7) & .$Genotype == "Δstemin" ~ "Cluster7, Δstemin",
    !.$Clusters %in% c(1,7) & .$Genotype == "Wild type" ~ "Other, Wild type",
    TRUE ~ "Other, ∆stemin"
  )) %>% {.$`Reprogramming by Genotype` <- factor(.$`Reprogramming by Genotype`,
                                                levels = c("Cluster 1,  Wild type",
                                                           "Cluster 7,  Wild type",
                                                           "Other, Wild type", 
                                                           "Cluster 1, Δstemin",
                                                           "Cluster 7, Δstemin",
                                                           "Other, ∆stemin"));.}

#######################################################
# Link Peaks
#######################################################

DefaultAssay(ppsn) <- "peaks"

# first compute the GC content for each peak
ppsn <- RegionStats(ppsn, genome = BSgenome.Ppatens.v6)

# link peaks to genes
ppsn <- LinkPeaks(
  object = ppsn,
  peak.assay = "peaks",
  expression.assay = "SCT",
  verbose = TRUE
)

save.image(file = paste0(SeuratObjPaths, saveName, ".RData"))
print("save3 done")


#######################################################
# Gene Activity
#######################################################

DefaultAssay(ppsn) <- "peaks"
geneActivity <- Signac::GeneActivity(ppsn,process_n = 8000)
ppsn[['activity']] <- CreateAssayObject(counts = geneActivity)
#Seurat::CreateSeuratObject(ppsn@assays$activity) %>% SeuratDisk::as.h5Seurat('ppsn_activity')
#SeuratDisk::Convert('ppsn_activity.h5seurat', dest = 'h5ad')

save.image(file = paste0(saveName, ".RData"))
print("save4 done")

#######################################################
#  TryCatch <~ tail
#######################################################
}, 
error = function(){},
finally = {
save.image(file = paste0(saveName, ".RData"))
}
)

rm(list = ls())
q()




#######################################################
# Differentially Accessible Peaks
#######################################################

DefaultAssay(ppsn) <- "peaks"
Idents(ppsn) <- "Clusters"

all_wilcox_da_peaks <- FindAllMarkers(ppsn, only.pos = TRUE, min.pct = 0.01, test.use = 'wilcox')
clust1_wilcox_da_peaks <- all_wilcox_da_peaks%>%filter(cluster == 1 & p_val < 0.05 & avg_log2FC > 1 & p_val_adj < 0.05)


all_LR_da_peaks <- FindAllMarkers(ppsn, only.pos = TRUE, min.pct = 0.01, test.use = 'LR')
clust1_LR_da_peaks <- all_LR_da_peaks%>%filter(cluster == 1 & p_val < 0.05 & avg_log2FC > 1 & p_val_adj < 0.05)

save.image(file = paste0(saveName, ".RData"))
print("save5 done")
#######################################################
# Motif Analysis
#######################################################

JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

pfm <- getMatrixSet(
x = JASPAR,
opts = list(collection = "CORE", tax_group = 'plants', all_versions = FALSE)
)

ppsn <- AddMotifs(
ppsn,
genome = BSgenome.Ppatens.v6,
pfm = pfm
)

clust1_LR_enriched.motifs <- FindMotifs(
object = ppsn,
features = rownames(clust1_LR_da_peaks)
)

clust1_wilcox_enriched.motifs <- FindMotifs(
object = ppsn,
features = rownames(clust1_wilcox_da_peaks)
)

ppsn <- RunChromVAR(
object = ppsn,
genome = BSgenome.Ppatens.v6)




save.image(file = paste0(saveName, ".RData"))


#SaveH5Seurat(ppsn, filename = paste0(SeuratObjPaths,"ppsn_Final.h5Seurat"), overwrite = T)


#######################################################
#  TryCatch <~ tail
#######################################################
}, 
error = function(){},
finally = {
save.image(file = paste0(saveName, ".RData"))
}
)

rm(list = ls())
q()