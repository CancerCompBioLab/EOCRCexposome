### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)



GLM_MRS <-  function(datameth, 
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
  
  # we only include traits that have at least 2 CpGs
  pheno <- pheno[as.numeric(ave(pheno$Trait, pheno$Trait, FUN=length)) >= 2, ]
  
  phenos <- unique(pheno$Trait)
  
  met_mat2 <- met_mat_M
  
  Scores <- matrix(ncol= length(phenos), nrow= ncol(met_mat2))
  colnames(Scores) <- phenos
  row.names(Scores) <- colnames(met_mat2)
  
  
  for (j in phenos){
      pheno1 <- pheno[pheno$Trait == j, ]
    
    # use only CpGs present in both
    DNAm <- as.data.frame(met_mat2[rownames(met_mat2) %in% pheno1$CpG, ])
    pheno1 <- as.data.frame(pheno1[pheno1$CpG %in% rownames(met_mat2), ])
    
    df2 <- DNAm[order(rownames(DNAm),decreasing=TRUE),]
    df3 <- pheno1[order(pheno1$CpG,decreasing=TRUE),]
    
    scores <- as.numeric(rep(NA,ncol(df2)))
    
    for (i in 1:ncol(df2)){
      scores[i]<-sum(df2[,i]*df3$Beta)}  
    
    Scores[ ,j] <- scores
    
  }
  
  Scores <- as.data.frame(Scores)
  
  library(dplyr)
  
  # The obtained scores are scaled
  Scores <- Scores %>% mutate_all(~(scale(.) %>% as.vector))
  
  Scorescl  <- merge(Scores, clinic, by.x= 0, by.y="ID")
  rownames(Scorescl) <- Scorescl$Row.names
  Scorescl$Row.names <- NULL
  
  
  # Save the scaled MRS 
  saveScores <- paste0("MRS/", PhenoCpG, datasetname, ".csv")
  write.csv(Scorescl, saveScores)

  GSE <- Scorescl
  GSE <- na.omit(GSE)
  GSVA <- GSE[, names(GSE) %in% pheno$Trait] 
  clinic2 <- GSE[, names(GSE) %in% names(clinic)]
  clinic2$ID <- rownames(clinic2)
  clinic2$EOCRC <- as.factor(clinic2$EOCRC) 
  
  if ("gender" %in% adjust){
    
    clinic2$gender[clinic2$gender == "FEMALE"] <- 1
    clinic2$gender[clinic2$gender == "MALE"] <- 0
    clinic2$gender <- as.factor(clinic2$gender)
    
    gen <- as.data.frame(table(clinic$gender, clinic$EOCRC))
    
    if (0 %in% gen$Freq) {
      adjust <- adjust[adjust != "gender"]  
    }  
    
  }
  
  if (identical(adjust, character(0))) {
    
    adjust <- NULL
    
  }
  
  
  if (is.null(adjust)){
    
        model <- "EOCRC ~ value"
    
  } 
  
  
  if(length(adjust)>0){
    
  adjust <- paste("+", adjust) 
    
    
      if(length(adjust)==1){
      
      b <- paste(adjust[1])}
    
    if(length(adjust)==2){
      
      b <- paste(adjust[1], adjust[2])}
      
      model <- paste("EOCRC ~ value",  b, sep ="")
    
  }    
  
  GSVA <- as.matrix(t(GSVA))
  COggpl <- reshape2::melt(GSVA)
  COggpl$Var2 <- as.character(COggpl$Var2)
  clinic2$ID <- as.character(clinic2$ID)
  
  COggpl <- dplyr::left_join(COggpl, clinic2, by=c("Var2" = "ID"))
  
  var <- unique(COggpl$Var1)
  
  GLM <- matrix(ncol=6, nrow=length(var))
  colnames(GLM) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "CI_l", "CI_H")
  row.names(GLM) <- var
  
  
  
  #
  # If a model gives a warning or does not generate the confidence intervals, we set the results to NA
  #
  
  myTryCatch <- function(expr) {
    warningMessage <- errorMessage <- NULL
    result <- tryCatch(
      expr,
      warning = function(w) {
        warningMessage <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      },
      error = function(e) {
        errorMessage <<- conditionMessage(e)
        NULL
      }
    )
    list(value = result, warning = warningMessage, error = errorMessage)
  }
  
  
  
  
  
  
  for (i in var){
    newggpl <- COggpl
    newggpl <- newggpl[which(newggpl$Var1 == i), ]
    
    warn1 <- myTryCatch(test <- glm(model, family = binomial, data = newggpl))
    
    if(is.null(warn1$warning)){ 
      results_df <-summary.glm(test)$coefficients
      results_df <- as.data.frame(t(results_df[2,]))
      GLM[i,1] <- results_df[1,1]
      GLM[i,2] <- results_df[1,2]
      GLM[i,3] <- results_df[1,3]
      GLM[i,4] <- results_df[1,4]
      
    } else {
      GLM[i,] <- NA
    }
    
    warn <- myTryCatch(CIl <-  confint(test)["value","2.5 %"])
    warn2 <- myTryCatch(CIh <- confint(test)["value","97.5 %"])
    
    if(is.null(warn$warning)){ 
      GLM[i,5] <- as.numeric(CIl)
    } else {
      GLM[i,5] <- NA
    }
    
    
    if(is.null(warn$warning)){ 
      GLM[i,6] <- as.numeric(CIh)
    } else {
      GLM[i,6] <- NA
    }
  }
  
  
  GLM2 <- as.data.frame(GLM)
  GLM2$study <- datasetname
  GLM2$name <- rownames(GLM2)
  
  datasetname2 <- paste("_", datasetname, sep="")
 
  
  if (!is.null(adjust)){
    
    covadj <- paste("_", cov, sep="")
    
  } else {
    
    covadj <- "NoAdjust"
  }  
  
  
  
  if(length(covadj)==1){
    b <- covadj[1]}
  
  if(length(covadj)==2){
    b <- paste(covadj[1], covadj[2], sep="")
  
 }
  
  covmod <- b
 
  saveglm <- paste0("GLM_MRS/", PhenoCpG, datasetname2, covmod, ".csv")
  
  write.csv(GLM2, saveglm, row.names = F)
  
  
}

