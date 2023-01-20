# 0. Download packages ########################################################
bio_pkgs = c("BayesSpace","SingleCellExperiment","scran","scater","BiocNeighbors","ComplexHeatmap")
BiocManager::install(bio_pkgs, update = T, force =TRUE)
install.packages("Seurat")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("patchwork")
install.packages("cowplot")
devtools::install_github("sqjin/CellChat")

# 1. Clustering Analysis - BayesSpace #########################################

## 0. Load packages ####
bio_pkgs = c("BayesSpace","SingleCellExperiment","scran","scater","BiocNeighbors","ComplexHeatmap")
invisible(lapply(c(bio_pkgs,"dplyr","Seurat","ggplot2","patchwork","cowplot"), function(x) library(x, character.only=TRUE))) 
set.seed(1234)

## 1. Load Data ####

### Load pre-built dataset ######
melanoma = getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")

# For R < 4.1 # 
library(assertthat)
library(RCurl)
getRDS <- function(dataset, sample, cache=TRUE) {
  
  url <- "https://fh-pi-gottardo-r-eco-public.s3.amazonaws.com/SpatialTranscriptomes/%s/%s.rds"
  url <- sprintf(url, dataset, sample)
  assert_that(url.exists(url), msg="Dataset/sample not available")
  
  if (cache) {
    bfc <- BiocFileCache()
    local.path <- bfcrpath(bfc, url)
  } else {
    local.path <- tempfile(fileext=".rds")
    download.file(url, local.path, quiet=TRUE, mode="wb")
  }
  
  readRDS(local.path)
}
melanoma = getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")


### Explore overall data ######
melanoma

### Explore row data ######
rowData(melanoma)

### Explore column data ######
colData(melanoma)

## 2. Preprocessing  ####

### Add metadata  ######
metadata(melanoma)$BayesSpace.data = list()
metadata(melanoma)$BayesSpace.data$platform = "ST"
metadata(melanoma)$BayesSpace.data$is.enhanced = FALSE

### Log normalize ######
melanoma = logNormCounts(melanoma)
melanoma@assays@data$logcounts[1:10,1:10]

## 3. Feature selection / Dimension reduction ####
### Feature selection ######
dec = scran::modelGeneVar(melanoma, assay.type="logcounts")
top = scran::getTopHVGs(melanoma, n=2000)

### Dimension reduction - PCA ######
melanoma  = scater::runPCA(melanoma, subset_row = top,ncomponents = 50)

### Preprocessing, Feature selection, Dimension reduction at once ######
# melanoma = spatialPreprocess(melanoma, platform="ST", n.PCs=50, n.HVGs=2000, log.normalize=TRUE)


## 4. Clustering  ####

### Select number of PCs to use ######
percent.var = attr(reducedDim(melanoma), "percentVar")
plot(percent.var[1:20], log="y", xlab="PC", ylab="Variance explained (%)", main = "Elbow Plot")
d = 7

### Select number of clusters to generate ######
melanoma = qTune(melanoma, qs=seq(2, 10), platform="ST",d=d)
qPlot(melanoma)
q = 4

### Clustering with BayesSpace ######
melanoma = spatialCluster(melanoma,q=q, d=d, gamma=2,platform="ST",nrep=10000,
                          save.chain=TRUE)
palette = c("purple", "red", "blue", "yellow", "darkblue")
clusterPlot(melanoma, palette=palette,color="black", size=0.1) +labs(title="BayesSpace")

# saveRDS(melanoma,'melanoma.rds')

### Enhancing resolution######
# melanoma.enhanced = spatialEnhance(melanoma, q=q, d=d, platform="ST", gamma=2,nrep=200000, verbose=TRUE, save.chain=TRUE,jitter_scale=3.5,jitter_prior=0.3)
melanoma.enhanced = 	readRDS(url("https://parkjooyoung99.github.io/KOGO_2023/data/KOGO_2023_spatial_transcriptomics_melanoma_enhanced_object.rds"))

melanoma.enhanced

clusterPlot(melanoma.enhanced,palette=palette,color="black", size=0.1) +labs(title="BayesSpace")

## 5. Annotation ####

### Define marker genes for plotting ######
markers = list()
markers[["Tumor"]] = c("PMEL")
markers[["Fibroblast"]] = c("COL1A1")
markers[["Macrophage"]] = c("CD14", "FCGR1A", "FCGR1B")
markers[["B-cell"]] = c("CD19", "MS4A1")
markers[["T-cell"]] = c("CD2", "CD3D", "CD3E", "CD3G", "CD7")

sum_counts = function(sce, features) {
  if (length(features) > 1) {
    colSums(logcounts(sce)[features, ])
  } else {
    logcounts(sce)[features, ]
  }
}
spot_expr = purrr::map(markers, function(xs) sum_counts(melanoma, xs))
enhanced_expr = purrr::map(markers, function(xs) sum_counts(melanoma.enhanced, xs))
plot_expression = function(sce, expr, title){ 
  featurePlot(sce, expr, color=NA) + viridis::scale_fill_viridis(option="A")+ labs(title=title, fill="Log-normalized\nexpression")
}
plot_expression_comparison = function(cell_type){ 
  spot.plot = plot_expression(melanoma,spot_expr[[cell_type]],"Spot")
  enhanced.plot = plot_expression(melanoma.enhanced, enhanced_expr[[cell_type]],"Enhanced")
  (spot.plot + enhanced.plot) + plot_annotation(title=cell_type, theme=theme(plot.title=element_text(size=18)))
}

### Rough cell type check ######
# hvgs = top[grep("^RP[LS]", top, invert=TRUE)]
# melanoma.enhanced = enhanceFeatures(melanoma.enhanced, melanoma,  model="xgboost",feature_names=hvgs, nrounds=0)
melanoma.enhanced = readRDS(url("https://parkjooyoung99.github.io/KOGO_2023/data/KOGO_2023_spatial_transcriptomics_melanoma_enhancedfeature_object.rds"))

melanoma.enhanced

p1 = plot_expression_comparison("Tumor")
p2 = plot_expression_comparison("Fibroblast")
p3 = plot_expression_comparison("Macrophage")
p4 = plot_expression_comparison("B-cell")
p5 = plot_expression_comparison("T-cell")

cowplot::plot_grid(p1,p2,p3,p4,p5, ncol=2)

### Find Marker genes for each cluster ######
sobj = Seurat::CreateSeuratObject(counts=logcounts(melanoma.enhanced), assay="Spatial", meta.data=as.data.frame(colData(melanoma.enhanced)))
sobj = Seurat::SetIdent(sobj, value = "spatial.cluster")
sobj@assays$Spatial@scale.data = sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
top_markers = Seurat::FindAllMarkers(sobj, assay="Spatial", slot="data", group.by="spatial.cluster", only.pos=TRUE) %>% group_by(cluster) %>% top_n(5, avg_log2FC)

head(top_markers,2)

### Visualize marker genes expression with heatmap ######
Seurat::DoHeatmap(sobj,  features = top_markers$gene, 	slot="scale.data", group.by = "spatial.cluster", 	group.colors=palette, angle=0, size=4, label = FALSE, 	raster=FALSE) + guides(col = FALSE)

### Annotate spots for enhaced object ######
coldata = colData(melanoma.enhanced) 
coldata$spatial.cluster.annotation = ifelse(coldata$spatial.cluster == 1, "Macrophage",ifelse(coldata$spatial.cluster == 2,"Fibroblasts",ifelse(coldata$spatial.cluster == 3, "Tumor","B cell")))
coldata$spatial.cluster.annotation = factor(coldata$spatial.cluster.annotation , levels = c('Macrophage','Fibroblasts','Tumor','B cell'))
colData(melanoma.enhanced)= coldata

clusterPlot(melanoma.enhanced, 	palette=palette, color="black", size=0.1,
            label = "spatial.cluster.annotation") +labs(title="Annotation")

# saveRDS(melanoma.enhanced, 'melanoma_enhanced.rds')

### Annotate spots for melanoma object ######
coldata = colData(melanoma) 
coldata$spatial.cluster.annotation = ifelse(coldata$spatial.cluster == 1, "Macrophage",ifelse(coldata$spatial.cluster == 2,"Fibroblasts",ifelse(coldata$spatial.cluster == 3, "Tumor","B cell")))
coldata$spatial.cluster.annotation = factor(coldata$spatial.cluster.annotation , levels = c('Macrophage','Fibroblasts','Tumor','B cell'))
colData(melanoma)= coldata

clusterPlot(melanoma, 	palette=palette, color="black", size=0.1,
            label = "spatial.cluster.annotation") +labs(title="Annotation")

# saveRDS(melanoma, 'melanoma.rds')

# Clear R environment for further analysis #############################################
rm(list = ls())
gc()
#############################################

# 2. Cell-cell interaction - CellChat ####

## 0. Load packages ####
invisible(lapply(c("CellChat", "patchwork"), function(x) library(x,character.only=TRUE))) 
set.seed(1234)

## 1. Load dataset ####

### Load data ######
visium.brain = readRDS(url("https://parkjooyoung99.github.io/KOGO_2023/data/KOGO_2023_spatial_transcriptomics_cortex_object.rds"))
visium.brain$manual_annotation = factor(visium.brain$manual_annotation, levels=c("Astro", "L2/3 IT", "L4", "L5 IT","L6 IT", "L6 CT", "L6b", "Oligo"))
Idents(visium.brain) = visium.brain$manual_annotation
colors = scPalette(nlevels(visium.brain))
names(colors) = c("Astro", "L2/3 IT", "L4", "L5 IT","L6 IT", "L6 CT", "L6b", "Oligo")
SpatialDimPlot(visium.brain, label = T, label.size = 3, cols = colors)

### Prepare input data for CellChat analysis ######
data.input = GetAssayData(visium.brain, slot = "data", assay = "SCT")
meta = data.frame(labels = Idents(visium.brain), row.names = names(Idents(visium.brain)))

unique(meta$labels) # check the cell labels

### Load spatial imaging information ######
spatial.locs = GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol")) 
scale.factors = jsonlite::fromJSON(txt = url("https://parkjooyoung99.github.io/KOGO_2023/data/scalefactors_json.json"))
scale.factors = list(spot.diameter=65, spot=scale.factors$spot_diameter_fullres, fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef
)

### Create a CellChat object ######
cellchat = createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)

cellchat

### Set the ligand-receptor interaction database ######
CellChatDB = CellChatDB.mouse
cellchat@DB = CellChatDB

## 2. Preprocessing ####

### Preprocessing the expression data for cell-cell communication analysis ######
cellchat = subsetData(cellchat)
future::plan("multiprocess", workers = 4) # do parallel
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)

## 3. Infer cell-cell communication network ####

### Compute the communication probability and infer cellular communication network
# cellchat = computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
#                              distance.use = TRUE, interaction.length = 200, 
#                              scale.distance = 0.01)
cellchat = readRDS(url("https://parkjooyoung99.github.io/KOGO_2023/data/KOGO_2023_spatial_transcriptomics_cortex_cellchat_computeCommunProb_object.rds"))

cellchat = filterCommunication(cellchat, min.cells = 10)

### Infer the cell-cell communication at a signaling pathway level ######
cellchat = computeCommunProbPathway(cellchat)

### Calculate the aggregated cell-cell communication network ######
cellchat = aggregateNet(cellchat)

### Visualization of the aggregated cell-cell communication network ######
groupSize = as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

## 4. Visualization ####

### Visualization of cell-cell communication network ######
pathways.show = c("CXCL") 

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

### Spatial Plot ######
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

### Compute the network centrality scores ######
cellchat = netAnalysis_computeCentrality(cellchat, slot.name = "netP")

par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

## 5. Find ligand-receptor pairs ####

### Identify ligand-receptor pairs associated with CXCL pathway ######
levels(cellchat@idents)
CellChat::netVisual_bubble(cellchat, sources.use =  c(3), targets.use =c(4),remove.isolate = FALSE, angle.x = 90,thresh = 0.05) + coord_flip()

# If you are encountering NA in row, modified function below might work #

netVisual_bubble_windows <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, color.heatmap = c("Spectral","viridis"), n.colors = 10, direction = -1, thresh = 0.05,
                                    comparison = NULL, group = NULL, remove.isolate = FALSE, max.dataset = NULL, min.dataset = NULL,
                                    min.quantile = 0, max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = F, color.text = NULL,
                                    title.name = NULL, font.size = 10, font.size.title = 10,show.legend = TRUE,
                                    grid.on = TRUE, color.grid = "grey90", angle.x = 90, vjust.x = NULL, hjust.x = NULL,
                                    return.data = FALSE){
  
  angle=c(0, 45, 90)
  hjust=c(0, 1, 1)
  vjust=c(0, 1, 0.5)
  vjust.x = vjust[angle == angle.x]
  hjust.x = hjust[angle == angle.x]
  color.use <- color.heatmap
  color.use <- rev(color.use)
  
  cells.level <- levels(object@idents)
  sources.use <- cells.level[sources.use]
  targets.use <- cells.level[targets.use]
  df.net <- subsetCommunication(object, slot.name = "net",
                                sources.use = sources.use, targets.use = targets.use,
                                signaling = signaling,
                                pairLR.use = pairLR.use,
                                thresh = thresh)
  
  df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")
  source.target <- paste(rep(sources.use, each = length(targets.use)), targets.use, sep = " -> ")
  source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
  
  df.net$pval[df.net$pval > 0.05] = 1
  df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
  df.net$pval[df.net$pval <= 0.01] = 3
  df.net$prob[df.net$prob == 0] <- NA
  df.net$prob.original <- df.net$prob
  df.net$prob <- -1/log(df.net$prob)
  
  idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 0)
  # rownames(df.net) <- df.net$interaction_name_2
  
  df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% unique(df.net$source)])
  df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% unique(df.net$target)])
  group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), levels(df.net$target), sep = " -> ")
  
  df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
  df.net <- with(df.net, df.net[order(interaction_name_2),])
  df.net$interaction_name_2 <- factor(df.net$interaction_name_2, levels = unique(df.net$interaction_name_2))
  cells.order <- group.names
  df.net$source.target <- factor(df.net$source.target, levels = cells.order)
  df <- df.net
  
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target),unique(df$source.target)))
  
  # g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
  #   geom_point(pch = 16) +
  #   theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x, color = 'black'),axis.text.y = element_text(color = 'black'),
  #         axis.title.x = element_blank(),
  #         axis.title.y = element_blank()) +
  #   scale_x_discrete(position = "bottom")
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() + theme(panel.grid.major = element_blank()) +
    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_x_discrete(position = "bottom")
  
  values <- c(1,2,3); names(values) <- c("p > 0.05", "0.01 < p < 0.05","p < 0.01")
  g <- g + scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)),labels = names(values)[values %in% sort(unique(df$pval))], name = "p-value")
  g <- g + scale_colour_gradientn(colors = rev(brewer.pal(11,"Spectral")), na.value = "white", limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                                    breaks = c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)), labels = c("min","max")) +
      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  g <- g + theme(text = element_text(size = font.size),plot.title = element_text(size=font.size.title)) +
    theme(legend.title = element_text(size = 8), legend.text = element_text(size = 6))
  
  g <- g + geom_hline(yintercept=seq(1.5, length(unique(df$interaction_name_2))-0.5, 1),lwd=0.1,colour=color.grid)
  g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  return(g)
}


netVisual_bubble_windows(cellchat, sources.use =  3, targets.use =c(4), angle.x = 90,thresh = 0.05, color.heatmap = 'viridis') + coord_flip()
