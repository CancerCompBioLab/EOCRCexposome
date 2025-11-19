### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Lea Lemler, Maria Butjosa-Espín, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)


build_gene_sets_from_metadata <- function(
    metadata_df,
    annotation_file,
    output_file = "Dosis_all.csv",
    base_path = "Results/MRS_Validation_IPC/",
    verbose = TRUE
) {
  suppressPackageStartupMessages({
    library(DESeq2)
    library(readr)
    library(dplyr)
    library(stringr)
  })
  
 
  counts_file_path <- function(plate_id) {
     paste0(base_path, "GSE262419_", plate_id, ".gz")
    
  }
  
  all_gene_sets <- list()
  
  for (trait_name in setdiff(unique(metadata_df$trait), "VEH")) {
    sub_meta <- metadata_df[metadata_df$trait == trait_name, ]
    plate_id <- unique(sub_meta$plate)
    cutoff_str <- unique(sub_meta$cutoff)
    directions_raw <- unique(sub_meta$directions)
    
    if (length(plate_id) != 1 || length(cutoff_str) != 1 || length(directions_raw) != 1) {
      warning(paste("⚠ Skipping trait", trait_name, ": ambiguous plate/cutoff/directions."))
      next
    }
    
    file_path <- counts_file_path(plate_id)
    if (!file.exists(file_path)) {
      warning(paste("❌ File not found:", file_path))
      next
    }
    
    if (verbose) cat("\n▶ Trait:", trait_name, "| Plate:", plate_id, "\n")
    
    
    pasCts = read.csv(gzfile(file_path), header = TRUE, row.names = "Genes")
    colnames(pasCts) <- gsub("\\.", "_", colnames(pasCts))
   
    QCct <- pasCts
    QCct$genes <- sub("_.*", "", rownames(QCct))
    
    gene_counts <-as.data.frame(QCct %>%
                                  group_by(genes) %>%
                                  summarise(across(where(is.numeric), sum)))
    rownames(gene_counts) <- gene_counts$genes
    gene_counts$genes <- NULL
    cts <- as.matrix(gene_counts)
    
    sub_df <- metadata_df[metadata_df$plate == plate_id &
                            metadata_df$trait %in% c("VEH", trait_name), ]
    
    rownames(sub_df) <- sub_df$sample_id
    sub_df <- sub_df[rownames(sub_df) %in% colnames(cts), ]
    cts <- cts[, colnames(cts) %in% rownames(sub_df)]
    sub_df <- sub_df[match(colnames(cts), rownames(sub_df)), ]
    
    sub_df$treatment <- sub("_.*", "", sub_df$treatment.ch1)
    sub_df$dose <- as.numeric(gsub(".*_(\\d+(\\.\\d+)?)_uM", "\\1", sub_df$treatment.ch1))
    sub_df$dose[grepl("VEH", sub_df$treatment)] <- 0
    sub_df$dose <- as.numeric(sub_df$dose)
    
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = sub_df, design = ~ dose)
    dds <- dds[rowMeans(counts(dds)) >= 10, ]
    dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
    
    dds <- DESeq(dds)
    res <- results(dds, name = "dose")
    res <- res[order(res$pvalue), ]
    
    annot <- read_tsv(gzfile(annotation_file), show_col_types = FALSE)
    
    # parse directions and cutoffs
    directions_list <- unlist(strsplit(directions_raw, ","))
    if ("all" %in% directions_list) {
      directions_list <- c("positive", "negative", "combined")
    }
    cutoff_list <- as.numeric(unlist(strsplit(cutoff_str, ",")))
    
    for (cutoff in cutoff_list) {
      for (dir in directions_list) {
        if (verbose) cat("  ⤷ Direction:", dir, "| Cutoff:", cutoff, "\n")
        
        pos_genes <- rownames(res)[res$pvalue < cutoff & res$log2FoldChange > 0]   
        neg_genes <- rownames(res)[res$pvalue < cutoff & res$log2FoldChange < 0]
        
        sig <- switch(dir,
                      positive = pos_genes,
                      negative = neg_genes,
                      combined = union(pos_genes, neg_genes),
                      rownames(res)[res$pvalue < cutoff]
        )
        
        symbols <- sub("_.*", "", sig)
        annots <- annot[annot$Symbol %in% symbols, ]
        annots$trait <- trait_name
        annots$direction <- dir
        annots$plate <- plate_id
        annots$cutoff <- cutoff
        annots$n_genes <- nrow(annots)
        
        if (nrow(annots) > 0) {
          key <- paste(trait_name, dir, "p", cutoff, sep = "_")
          all_gene_sets[[key]] <- annots[, c("plate", "trait", "direction", "cutoff", "n_genes", "Symbol", "EnsemblGeneID")]
        }
      }
    }
  }
  
  output <- paste0(base_path, output_file)
  
  gene_sets_df <- bind_rows(all_gene_sets)
  write.csv(gene_sets_df, output, row.names = FALSE)
  
  if (verbose) cat("\n✔ Gene sets saved to:", output_file, "\n")
  return(gene_sets_df)
}


library(GSVA)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)

## Download data from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE262419
# Data is not provided in this github

coldata <- read.csv("pDATA.csv", row.names= 1)
coldata$treatment <- sub("_.*", "", coldata$treatment.ch1)
coldata$sample_id <- rownames(coldata)

# prepare input file:

#   | sample\_id | treatment.ch1              | plate   | trait    | cutoff     | directions | 
#   | ---------- | -------------------------- | ------- | -------- | ---------- | ---------- | 
#   | S1         | treatment: Picloram\_1\_uM | Plate10 | Picloram | 0.01,0.005 | all        | 
#   | S2         | treatment: VEH\_0\_uM      | Plate10 | VEH      | NA         | NA         | 

meta <- coldata[,c("sample_id", "treatment.ch1", "description", "treatment")]  
names(meta) <- c("sample_id", "treatment.ch1", "plate", "trait")  

# Select the exposures 
trait <- c("Acetochlor", "Atrazine", "Dicamba", "Heptachlor", "Lindane", "Metolachlor", "Picloram")

clean_trait_name <- function(trait_vec) {
  cleaned <- gsub("[ ,\\-]", ".", trait_vec)                # Replace spaces, commas, hyphens with dots
  cleaned <- gsub("\\.+", ".", cleaned)                     # Replace multiple dots with single dot
  cleaned <- sub("^\\.+|\\.+$", "", cleaned)                # Remove leading/trailing dots
  
# If trait starts with a digit, prepend "X"
cleaned <- ifelse(grepl("^[0-9]", cleaned), paste0("X", cleaned), cleaned)

  return(cleaned)
}

meta$trait_clean <- clean_trait_name(meta$trait)
meta$directions <- "combined"                  # negative, positive, combined, all
meta$cutoff <- "0.005" # Uses non-adjusted p-value, can be changed in rows 99/100


# annotation_file = https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz

gene_sets <- build_gene_sets_from_metadata(
  metadata_df = meta,
  annotation_file = "Human.GRCh38.p13.annot.tsv.gz",
  base_path = "EOCRCexposome/Results/",
  output_file = "Dosis_all.csv"
)


# download TCGA COAD gene expression data using https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html

COADexp<-readRDS("COAD_TCGA_gene_expression.RDS")

TCGA_expr<-as.data.frame(assays(COADexp)$tpm_unstrand)
TCGA_expr$genes <- rownames(TCGA_expr)
TCGA_expr$genes <- gsub("\\..*","",TCGA_expr$genes)

tcga_expr = TCGA_expr
id_col_expr = "genes"
tcga_expr <- tcga_expr[complete.cases(tcga_expr[[id_col_expr]]), ]


tcga_expr <- tcga_expr[! rownames(tcga_expr) %like% "PAR_", ]
tcga_expr <- tcga_expr[!duplicated(tcga_expr[[id_col_expr]]), ]
rownames(tcga_expr) <- tcga_expr[[id_col_expr]]
tcga_expr[[id_col_expr]] <- NULL
expr_mat <- as.matrix(tcga_expr)
mode(expr_mat) <- "numeric"
expr_mat2 <- as.data.frame(expr_mat)
expr_mat3 <- as.matrix(expr_mat2)

#filter expression, and log2 transform

expr_mat4 <- expr_mat3[which(rowSums(expr_mat3>=10) >50 & rowVars(expr_mat3)>0.01 ), ]
expr_mat4 = t(apply(expr_mat4,1,function(x) log2(x+1)))

gene_sets_df <- read.csv("Results/MRS_Validation_IPC/Dosis_all.csv", header=T)


##### PICLORAM # other pesticides included in the file: Acetochlor, Atrazine, Dicamba, Heptachlor, Lindane, and Metolachlor
gene_sets_df_piclo = gene_sets_df[which(gene_sets_df$trait=="Picloram"),] 

gene_sets_df_piclo_COAD <- gene_sets_df_piclo[gene_sets_df_piclo$EnsemblGeneID %in% rownames(expr_mat4), ]

gene_sets_df_piclo_COAD$gene_set <- paste(gene_sets_df_piclo_COAD$trait, gene_sets_df_piclo_COAD$direction, "p", gene_sets_df_piclo_COAD$cutoff, sep = "_")
gene_list_piclo_COAD <- split(gene_sets_df_piclo_COAD$EnsemblGeneID, gene_sets_df_piclo_COAD$gene_set)


k <-37  #change this to reflect the number of gene overlapping between the gene set and TCGA-COAD = Differs per *pesticide*!

set.seed(569)
t10k <- lapply(1:1000, function(x) sample(rownames(expr_mat4), k, replace = F))
names(t10k)=paste0("random_",1:1000)

ssgsea_mat <- gsva(expr_mat4, c(gene_list_piclo_COAD,t10k), method = "ssgsea",BPPARAM=BiocParallel::MulticoreParam(20))

id1 = substr(colnames(ssgsea_mat),14,16)
ssgsea_mat1 = ssgsea_mat[,id1=="01A"]
colnames(ssgsea_mat1)=substr(colnames(ssgsea_mat1),1,12)

clin <- read.csv("Results/MRS/GW_COAD.csv", header=T) # Other pesticides use P1E5_COAD.csv due to low significant CpGs in the GW MRS

head(clin)
clintr <- clin[, colnames(clin) %in%  c("ID", "PICLORAM", "age") ]  # Acetochlor    Atrazine     Dicamba  Heptachlor     Lindane Metolachlor    Picloram 

ssgsea_df = as.data.frame(t(ssgsea_mat1))
ssgsea_df$SampleID <- rownames(ssgsea_df)

gc <- merge(clintr, ssgsea_df, by = "SampleID") 
#write.csv(gc, "Results/MRS_Validation_IPC/permutation/Picloram_dosis_P0_005_SCORES.csv", row.names = F)

ggplot(gc, aes(x=Picloram_combined_p_0.005))+
  geom_histogram()+
  theme_bw()+
  labs(x="Picloram ssGSEA score in TCGA-COAD")


num_cols <- setdiff(colnames(gc), c("SampleID", "Picloram","age"))

# Compute correlations + p-values
get_corr <- function(var) {
  x <- gc[[var]]
  spear <- suppressWarnings(cor.test(gc$Picloram, x, method = "spearman",use="complete.obs"))
  data.frame(
    variable   = var,
    spearman_r = spear$estimate,
    spearman_p = spear$p.value
  )
}

# Apply across all variables
results_df <- do.call(rbind, lapply(num_cols, get_corr))

df <- results_df %>%
  mutate(
    rank_spearman_p = rank(spearman_p, ties.method = "min"),
    rank_spearman_r = rank(-abs(spearman_r), ties.method = "min")
  )

#pdf("Figure3F.pdf",width = 2, height = 1.75)
ggplot(gc, aes(x=scale(Picloram),y=Picloram_combined_p_0.005))+
  geom_point(color = "#98A7D4", size=0.75)  + 
  geom_smooth(method='lm', formula= y~x)+
  stat_cor(method = "spearman", size = 2)+
   labs(
    x = "Picloram MRS-GW",
    y = "Picloram ssGSEA (37 genes)"
  ) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

#dev.off()


#pdf("Figure3H.pdf",width = 2, height = 1.45)
ggplot(df, aes(x = "", y = abs_spearman)) +
  geom_violin(fill = "#4C72B7", color = "gray30", alpha = 0.6, width = 0.6) +
  geom_boxplot(fill = "#4C72B7", width=0.1, alpha = 0.6, outlier.size =0.4)+
  annotate("point", x = 1, y = highlight_val, color = "#ae3960", size = 2, shape = 8, stroke=1) +
  annotate("text", x = 1.2, y = highlight_val, label = "Picloram", 
           color = "#ae3960",  size = 1.5) +
  labs(
    x = NULL,
    y = expression("|Spearman's rho|")) +
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()
  )


#dev.off()

