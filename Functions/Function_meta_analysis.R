### Script to reproduce results of:
### Epigenetic Fingerprints Link Early-Onset Colon and Rectal Cancer to Pesticide Exposure  
### Silvana C.E. Maas, Iosune Baraibar, Odei Blanco Irazuegui, Elena Elez, Jose A. Seoane 
### author: Silvana C.E. Maas (silvanamaas at vhio.net)



meta<-  function(method, 
                       PhenoCpG, 
                       datasets,
                       type
){
  
  
  
  if(method=="MRS"){
    setwd("Results/GLM_MRS")  # Change this if you want to use your own data
  }
  
  
  if(PhenoCpG=="GW"){
    files <- list.files(pattern = "GW")
  }
  
  if(PhenoCpG=="P1E5"){
    files <- list.files(pattern = "P1E5")
    files
  }
  
  if(PhenoCpG=="F001"){
    files <- list.files(pattern = "F001")
    files
    
  }
  
  if(PhenoCpG=="F005"){
    files <- list.files(pattern = "F005")
    files
    
  }
  
  if(PhenoCpG=="F01"){
    files <- list.files(pattern = "F01")
    files
  }
  
  data <- paste(datasets, collapse="|")
  files2 <- files[grep(data, files)]

  tables <- lapply(files2, read.csv, header = TRUE)
  comb <- do.call(rbind , tables)
  names(comb) <- c("beta", "se", "z.value", "pvalue", "CI_L","CI_H", "study", "name") 
  comb <- comb[complete.cases(comb),]
  
  name <- unique(comb$name)
  meta <- matrix(ncol=6, nrow=length(name))
  colnames(meta) <- c("Beta", "CI_L","CI_H", "pval", "study", "name") 
  row.names(meta) <- name
  
  library(dmetar)
  library(meta)
  
  
  for (i in name){
    newggpl <- comb
    newggpl <- newggpl[which(newggpl$name == i), ]
    
    title <- paste(i, "OR",sep="_")
    
    m.gen_bin <- metagen(TE = beta,
                         seTE = se,
                         lower = CI_L,
                         upper = CI_H,
                         studlab = study,
                         data = newggpl,
                         sm = "OR",
                         method.tau = "PM",
                         fixed = FALSE,
                         random = TRUE,
                         title = title)
    
    
    meta[i,1] <- as.numeric(m.gen_bin$TE.common)
    meta[i,2] <- as.numeric(m.gen_bin$lower.common)
    meta[i,3] <- as.numeric(m.gen_bin$upper.common)
    meta[i,4] <- as.numeric(m.gen_bin$pval.common)
    st <- unique(newggpl$study)
    meta[i,5] <- paste(st, collapse=",")
    meta[i,6] <- newggpl[1,8]
    
  
  }  
  
  meta <- as.data.frame(meta)
  outpM <- paste("meta_analysis/", method, type, PhenoCpG, ".csv", sep="_") 
  write.csv(meta, outpM, row.names = F)
  
}