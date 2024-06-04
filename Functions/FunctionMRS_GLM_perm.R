### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)


GLM_MRS_perm <-  function(datameth, 
                     dataclinic, 
                     histodiag = NULL,
                     PhenoCpG, 
                     datasetname, 
                     adjust = NULL
){
  
  
  met_mat_M  <- datameth
  clinic <- dataclinic 
  cov <- adjust
  
  
  
  
  if (!is.null(histodiag)){
    
    clinic <- clinic[clinic$histologic_diagnosis == histodiag, ]
    
  } 
  
  
  clinic <- clinic[complete.cases(clinic$EOCRC), ]
  met_mat_M <- met_mat_M[, colnames(met_mat_M) %in% clinic$ID]
  
  clinic <- clinic[clinic$ID %in% colnames(met_mat_M), ]
  
  CpG <- as.data.frame(rownames(met_mat_M))
  names(CpG)[1] <- "CpG"
  
  datasetname2 <- paste(PhenoCpG, datasetname, cov, sep = "_")
  
  # Open the output file from the GLM_MRS function
  saveglm <- paste("GLM_MRS/", datasetname2, ".csv", sep="")
  GLM_disc <-read.csv(saveglm, header=T)
  
  # Permutation is only done for triats that were significant in the discovery data set
  sign <- GLM_disc[GLM_disc$Pr...z.. <= 0.05, ]
  
  signs <- sign$name
  
  
  
  if(PhenoCpG=="GW"){
    Pheno = read.csv(file="InputData/CombGW.csv", header=TRUE)
    
  }
  
  if(PhenoCpG=="P1E5"){
    Pheno = read.csv(file="InputData/CombP1E5.csv", header=TRUE)
    
  }
  
  if(PhenoCpG=="F001"){
    Pheno = read.csv(file="InputData/CombF001.csv", header=TRUE)
    
  }
  
  if(PhenoCpG=="F005"){
    Pheno = read.csv(file="InputData/CombF005.csv", header=TRUE)
    
  }
  
  if(PhenoCpG=="F01"){
    Pheno = read.csv(file="InputData/CombF01.csv", header=TRUE)
    
  }
  
  pheno <- Pheno[,c("CpG", "Beta", "Trait")]
  
  pheno <- pheno[complete.cases(pheno$Trait), ]
  
  
  for (a in signs){

    trait <- as.data.frame(pheno[pheno$Trait == a,])
    
    # We exclude per trait the CpGs that were included in the original MRS
    CpG2 <- as.data.frame(CpG[!CpG$CpG %in% trait$CpG,])
    names(CpG2)[1] <- "CpG"
    
    k <- nrow(trait)
    t10k <- lapply(1:10000, function(x) sample(CpG2$CpG, k, replace = F))
    
    newlist <- append(list(trait$CpG), t10k)
    names(newlist) <- c("MRS", 1:10000)
    
    dflist <- as.data.frame(newlist)
    
    met_mat2 <- met_mat_M

    
  Scores <- matrix(ncol= ncol(dflist), nrow= ncol(met_mat2))
  colnames(Scores) <- colnames(dflist)
  row.names(Scores) <- colnames(met_mat2)
  
  labels <- colnames(dflist)
  
  for (j in labels){

    
    pheno1 <- trait
    # use only CpGs present in both
    DNAm <- as.data.frame(met_mat2[rownames(met_mat2) %in% dflist[ ,j],])
    pheno1$CpG2 <- dflist[ ,j]
    
    df2 <- DNAm[order(rownames(DNAm),decreasing=TRUE),]
    df3 <- pheno1[order(pheno1$CpG2,decreasing=TRUE),]
    
    scores <- as.numeric(rep(NA,ncol(df2)))
    
    for (m in 1:ncol(df2)){
      scores[m]<-sum(df2[,m]*df3$Beta)}  
    
    Scores[ ,j] <- scores
    
  }
  
  Scores <- as.data.frame(Scores)
  library(dplyr)
  
  Scores <- Scores %>% mutate_all(~(scale(.) %>% as.vector))

  Scorescl  <- merge(Scores, clinic, by.x= 0, by.y="ID")
  rownames(Scorescl) <- Scorescl$Row.names
  Scorescl$Row.names <- NULL
  

  GSE <- Scorescl
  GSE <- na.omit(GSE)
  GSVA <- GSE[, names(GSE) %in% labels] 
  
  
  clinic2 <- GSE[, names(GSE) %in% names(clinic)]
  clinic2$ID <- rownames(clinic2)
  clinic2$EOCRC <- as.factor(clinic2$EOCRC) 
  
    clinic2$gender[clinic2$gender == "FEMALE"] <- 1
    clinic2$gender[clinic2$gender == "MALE"] <- 0
    clinic2$gender <- as.factor(clinic2$gender)
    
      model <- paste("EOCRC ~ value + gender")
    
  GSVA <- as.matrix(t(GSVA))
  COggpl <- reshape2::melt(GSVA)
  COggpl$Var2 <- as.character(COggpl$Var2)
  clinic2$ID <- as.character(clinic2$ID)
  
  COggpl <- dplyr::left_join(COggpl, clinic2, by=c("Var2" = "ID"))
  
  var <- unique(COggpl$Var1)
  
  GLM <- matrix(ncol=4, nrow=length(var))
  colnames(GLM) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  row.names(GLM) <- var

  for (k in var){
    newggpl <- COggpl
    newggpl <- newggpl[which(newggpl$Var1 == k), ]
    test <- glm(model, family = binomial, data = newggpl) 
      results_df <-summary.glm(test)$coefficients
      results_df <- as.data.frame(t(results_df[2,]))
      GLM[k,1] <- results_df[1,1]
      GLM[k,2] <- results_df[1,2]
      GLM[k,3] <- results_df[1,3]
      GLM[k,4] <- results_df[1,4]
     
}
  
  
  GLM2 <- as.data.frame(GLM)
  GLM2$study <- datasetname
  GLM2$name <- rownames(GLM2)
  
  saveglm <- paste0("GLM_MRS_perm/", datasetname2, a,".csv")
  
  write.csv(GLM2, saveglm, row.names = F)
  
  
}
}
