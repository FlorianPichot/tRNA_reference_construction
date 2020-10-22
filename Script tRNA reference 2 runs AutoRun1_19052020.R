######## tRNA reference creation by sequence clustering ##########
# version script YM 19-05-2020  run 1 in automatic mode
# THREE STEPS Step 1 run REMOVE DUPLICATES ONLY
#SECOND COLLAPS SEQUENCES
#STEP THREE CREATE SINGLE FILE
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)                            }
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)                     }

library(seqinr)
library(msa)  
rm(list=ls(all=TRUE))

pdf.width=15
pdf.height=15

#BLOCK READ FILE  COPY-PASTE from WIN KEYBOARD
extDataDir <- gsub ( "\\\\", "/", readClipboard () ) #READ FROM CLIPBOARD AND send to extDataDir
setwd(extDataDir)

#GIVE NAME of Fasta FILE
Fasta_file_RNApattern <-"sacCer3-mature-tRNAs_16052020"

Fasta_file_RNA <- dir(paste0(extDataDir), pattern = paste0(Fasta_file_RNApattern), full.names = T)
Fasta_file_RNA #VERIFICATION of the file name

ResultsDir <- paste0(extDataDir,"/counting_distances_lim8")
dir.create(ResultsDir, showWarnings = FALSE)  #CREATE FOLDER FOR tRNA analysis
setwd(ResultsDir)

set.seed(123)

AAlist <- list("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "SeC")

#READING FASTA file and transformation to DATAFRAME
# Maybe use "Biostrings" package? <-- Already done by library(msa)

## Functions ##
readinteger <- function() { 
  n <- readline(prompt="Enter nc, the number of clusters : ")
  if(!grepl("^[0-9]+$",n))
  {
    return(readinteger())
  }
  return(as.integer(n))
                          }

#EXTRACT AND CHECK ORGANISM NAME
Fasta_RNA <-read.fasta(file = Fasta_file_RNA)
head(Fasta_RNA)
RNAstring <- getName(Fasta_RNA)[1]
print(RNAstring)
organism_name <- unlist(strsplit(RNAstring,"tRNA"))[1]     # something like"Saccharomyces_cerevisiae_"  #TAKE FROM FASTA file
print(organism_name)


## Run 1 ONLY TO REMOVE DUPLICATES##
for (aa in AAlist) {
 
  Fasta_RNA_AA <- Fasta_RNA[grep(pattern = aa, getName(Fasta_RNA))]
  print (aa)
  if (length (Fasta_RNA_AA) > 1) 
    {
    pattern <- character()
    for (i in 1:length(Fasta_RNA_AA)){pattern <- c(pattern, c2s(Fasta_RNA_AA[[i]]))}
      
   Dist_matrix <- matrix(nrow = length(Fasta_RNA_AA), ncol = length(Fasta_RNA_AA)) #EMPTY DISTANCE MATRIX
    #FILL DISTANCE MATRIX by calculation of distances type = "overlap", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE), gapOpening = 1, gapExtension = 0
    for (i in 1:length(Fasta_RNA_AA)){
      alignment <- pairwiseAlignment(DNAStringSet(pattern),DNAStringSet(c2s(getSequence(Fasta_RNA_AA[[i]]))), scoreOnly = TRUE, type = "overlap", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE), gapOpening = 1, gapExtension = 0) 
      Dist_matrix[,i] <- alignment
      Dist_matrix[,i] <- getLength(Fasta_RNA_AA[[i]]) - Dist_matrix[,i]
    }
    #Dist_matrix may be assymetric!
    print(Dist_matrix)
    print (unique(Dist_matrix[1,]))
    print(unique(Dist_matrix[nrow(Dist_matrix),]))
    
    #Create heatmap to visualize and to save                                
   # heatmap.2(Dist_matrix, key = FALSE, cellnote = Dist_matrix, Rowv=T, Colv="Rowv", notecex=0.8, notecol="black", na.color=par("bg"), dendrogram="row", density.info="none",  labRow=getName(Fasta_RNA_AA), labCol=getName(Fasta_RNA_AA), margins = c(8, 8), trace=c("none"))

    #SAVE file in pdf   
    pdf(paste0("Step1_Heatmap2_", aa,".pdf"),width=pdf.width,height=pdf.height, paper='special') 
    heatmap.2(Dist_matrix, key = FALSE, cellnote = Dist_matrix, Rowv=T, Colv="Rowv", notecex=0.8, notecol="black", na.color=par("bg"), dendrogram="row", density.info="none",  labRow=getName(Fasta_RNA_AA), labCol=getName(Fasta_RNA_AA), margins = c(15, 15), trace=c("none"))
    dev.off()
    #DEFINE number of clusters nc as max number of unique distances in the 1st and last lines of the heatmap
    nc <-max(length(unique(Dist_matrix[1,])), length(unique(Dist_matrix[nrow(Dist_matrix),])))
    
    #select the right number of clusters define nc value
    #print (paste0("tRNAs : ", aa,"  Suggested number of clusters : ", nc ,"  confirm or enter another number"))
    #nc= readinteger()     # Ne prendre que les clusters avec 0 IF OK FOR AUTOMATIC nc definition, not necessary anymore
    
    #graphics.off()
    
    d<-as.dist(Dist_matrix) #converrt to distances
    #attribute cluster numbers
    if (nrow(Dist_matrix) > nc) {
      fit <- hclust(d, method="ward.D2")
      fit.cluster <- cutree(fit, k=nc) # cut tree into nc clusters
      Dist_matrix <- data.frame(Dist_matrix, fit.cluster) 
                                }
    else {
      Dist_matrix <- data.frame(Dist_matrix, fit.cluster = 1:nc)
         }
    CluSeq <- data.frame(names = getName(Fasta_RNA_AA), ClusterGroup = Dist_matrix$fit.cluster)
    CluSeq <- CluSeq[order(CluSeq$ClusterGroup),] #ORDER CLUSTERS BY Cluster number
    
    if (nrow(CluSeq) != nc){
      CluSeqList <- unstack(CluSeq)  #SPLIT into groups
                           }
    else {
      CluSeqList <- as.list(as.character(CluSeq$names))
         }
    
    ## Cluster Treatment FUSION OF SEQUENCES##
    for (c in 1:length(CluSeqList)) {
      count <- integer()
      for (r in CluSeqList[[c]]){		
        count <- c(count, grep(paste0("^",r,"$"), getName(Fasta_RNA_AA)))
      }
      Fasta_RNA_col <- Fasta_RNA_AA[count]
      #WRITE LIST OF MERGED SEQUENCES
      write.fasta(Fasta_RNA_col, names = getName(Fasta_RNA_col), file.out = paste0("Step1_Collapsed_sequences_", aa, "_cluster",c,".fa"))
      
      #COLLAPS SEQUENCES in ONE entry
      alignment <- readDNAStringSet(paste0("Step1_Collapsed_sequences_", aa, "_cluster",c,".fa"), seek.first.rec=TRUE)
      
      if (length (alignment)>1) {
        alignment2 <- msaClustalW(alignment, 
                                  type = 'DNA',
                                  cluster = "upgma",              # matrice de distance --> arbre enraciné
                                  substitutionMatrix = "iub")     # prend en compte la notation IUPAC. + 1.9 si match, 0 si non.
        
        consensusSeq <- consensusString(consensusMatrix(alignment2))
        #consensusSeq
        #EVENTUALLY REMOVE SEQUENCE? SHOULD NOT HAPPEN HERE, TO REMOVE?
        if (countPattern(pattern = "?", subject = c2s(getSequence(consensusSeq))) >= 10) {
          write.fasta(Fasta_RNA_col, names = getName(Fasta_RNA_col), file.out = paste0("Step1_Removed_sequences_", aa, "_cluster",c,".fa"))
        }
        else {
          consensusSeq <- gsub(pattern = "[?]", "N", consensusSeq)
          consensusSeq <- gsub(pattern = "[-]", "", consensusSeq)
          consensusSeq <- as.SeqFastadna(s2c(consensusSeq), Annot = paste(gsub(organism_name, "", getName(Fasta_RNA_col)), collapse = "="))
          
          write.fasta(consensusSeq, names = getAnnot(consensusSeq), file.out = paste0("Step1_Sequence_NoDup_", aa, "_collapsed_cluster",c,".fa"))
        }}
      else {
        write.fasta(Fasta_RNA_col, names = gsub(organism_name, "", getName(Fasta_RNA_col)), file.out = paste0("Step1_Sequence_NoDup_", aa, "_untouched_cluster",c,".fa"))
      }
    }
  
  }
   else {if (length(Fasta_RNA_AA) != 0) {write.fasta(Fasta_RNA_AA, names = paste(gsub(organism_name, "", getName(Fasta_RNA_AA)), collapse = "|"), file.out = paste0("xStep2_UniqueSequence_NoDup_", aa, "_untouched__cluster",c,".fa"))}
  }
}

#####################################################################################################################
#READ AND CREATE NON-RUDUNDANT FASTA from all Step1_Sequence_NoDup_ files
#MERGE OF SEQUENCE in one Fasta file take by aminoacid

out.file<-NULL
file.names <- dir(ResultsDir, pattern = paste0("^Step1_Sequence_NoDup"), full.names = T)

for(i in 1:length(file.names)){
  
  file <- read.table(file.names[i],header=F, sep=";", stringsAsFactors=FALSE)
  out.file <- rbind(out.file, file)
}

print (i)

setwd(extDataDir)
write.table(out.file, file = paste0("Merge_Step1_NoDup",organism_name, Sys.Date(),".fa"),sep=";", col.names = F,
            row.names = FALSE, quote= F,fileEncoding="windows-1252")


#####################################################################################################################
## 2nd Run ## MERCGE NOW TO CLUSTERS
#WILL USE NON-repetitive sequence obtained during the first step "Step1_Sequence_NoDup_"
#####################################################################################################################

setwd(ResultsDir)
Fasta_RNA <- list()

for (fi in list.files(pattern = "^Step1_Sequence_NoDup_")){
  Fasta_RNA <-c(Fasta_RNA, read.fasta(file = fi))
                                                          }
head(Fasta_RNA) #CREATE AND VERIFY fasta file 

for (aa in AAlist) {
print (aa)
  Fasta_RNA_AA <- Fasta_RNA[grep(pattern = aa, getName(Fasta_RNA))] 
  
  if (length (Fasta_RNA_AA) > 1) {
    pattern <- character()
    for (ii in 1:length(Fasta_RNA_AA)){
      pattern <- c(pattern, c2s(Fasta_RNA_AA[[ii]]))
    }
    
    Dist_matrix <- matrix(nrow = length(Fasta_RNA_AA), ncol = length(Fasta_RNA_AA))  #MAKE EMPTY MATRIX
    #FILL DISTANCE MATRIX by calculation of distances type = "overlap", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE), gapOpening = 1, gapExtension = 0
    for (i in 1:length(Fasta_RNA_AA)){
      alignment <- pairwiseAlignment(DNAStringSet(pattern),DNAStringSet(c2s(getSequence(Fasta_RNA_AA[[i]]))), scoreOnly = TRUE, type = "overlap", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE), gapOpening = 1, gapExtension = 0) 
      Dist_matrix[,i] <- alignment
      Dist_matrix[,i] <- getLength(Fasta_RNA_AA[[i]]) - Dist_matrix[,i]
    }
    #Dist_matrix may be assymetric!
    print(Dist_matrix)
    print(rowSums(Dist_matrix>8))
    #Order matrix by heatmap
    heatmap.2(Dist_matrix, key = FALSE, cellnote = Dist_matrix, Rowv=T, Colv="Rowv", notecex=0.8, notecol="black", na.color=par("bg"), dendrogram="row", density.info="none",  labRow=getName(Fasta_RNA_AA), labCol=getName(Fasta_RNA_AA), margins = c(12, 12), trace=c("none"))
    
    #SAVE file in pdf   
    pdf(paste0("Step2_Heatmap2_", aa,".pdf"),width=pdf.width,height=pdf.height, paper='special') 
    heatmap.2(Dist_matrix, key = FALSE, cellnote = Dist_matrix, Rowv=T, Colv="Rowv", notecex=0.8, notecol="black", na.color=par("bg"), dendrogram="row", density.info="none",  labRow=getName(Fasta_RNA_AA), labCol=getName(Fasta_RNA_AA), margins = c(12, 12), trace=c("none"))
    dev.off()
    
    #DEFINE number of clusters nc as number distances >8 in all lines +1 (for Zero cluster)
    nc <-1+ min(rowSums(Dist_matrix>8))
  
    #select the right number of clusters define nc value
    
    print (paste0("tRNAs : ", aa,"  Suggested number of clusters : ", nc ,"  confirm or enter another number"))
    nc= readinteger()     # Ne prendre que les clusters avec 0 IF OK FOR AUTOMATIC nc definition, not necessary anymore
    
    graphics.off()
   
    ## Clustering ##
    d<-as.dist(Dist_matrix) #converrt to distances
    print (d)
    #attribute cluster numbers   
    if (nrow(Dist_matrix) > nc) {
      fit <- hclust(d, method="ward.D2")
      fit.cluster <- cutree(fit, k=nc) # cut tree into nc clusters
      Dist_matrix <- data.frame(Dist_matrix, fit.cluster) 
    }
    else {
      Dist_matrix <- data.frame(Dist_matrix, fit.cluster = 1:nc)
    }
    CluSeq <- data.frame(names = getName(Fasta_RNA_AA), ClusterGroup = Dist_matrix$fit.cluster)
    CluSeq <- CluSeq[order(CluSeq$ClusterGroup),]
 
    if (nrow(CluSeq) != nc){
      CluSeqList <- unstack(CluSeq)
    }
    else {
      CluSeqList <- as.list(as.character(CluSeq$names))
    }
    
    ## Cluster Treatment ##
    for (c in 1:length(CluSeqList)) {
      count <- integer()
      for (r in CluSeqList[[c]]){		
        count <- c(count, grep(paste0("^",r,"$"), getName(Fasta_RNA_AA)))
      }
      Fasta_RNA_col <- Fasta_RNA_AA[count]
      
      
      
      write.fasta(Fasta_RNA_col, names = getName(Fasta_RNA_col), file.out = paste0("Step2_Collapsed_sequences_", aa, "_cluster",c,".fa"))
      
      alignment <- readDNAStringSet(paste0("Step2_Collapsed_sequences_", aa, "_cluster",c,".fa"), seek.first.rec=TRUE)
      
      if (length (alignment)>1) {
        alignment2 <- msaClustalW(alignment, 
                                  type = 'DNA',
                                  cluster = "upgma",              # matrice de distance --> arbre enraciné
                                  substitutionMatrix = "iub")     # prend en compte la notation IUPAC. + 1.9 si match, 0 si non.
        
        consensusSeq <- consensusString(consensusMatrix(alignment2))
        consensusSeq
        
        if (countPattern(pattern = "?", subject = c2s(getSequence(consensusSeq))) >= 10) {
          write.fasta(Fasta_RNA_col, names = getName(Fasta_RNA_col), file.out = paste0("Step2_Removed_sequences_", aa, "_cluster",c,".fa"))
        }
        else {
          consensusSeq <- gsub(pattern = "[?]", "N", consensusSeq)
          consensusSeq <- gsub(pattern = "[-]", "", consensusSeq)
          consensusSeq <- as.SeqFastadna(s2c(consensusSeq), Annot = paste0(aa, "_cluster",c," consensus sequence ", paste(gsub(organism_name, "", getName(Fasta_RNA_col)), collapse = "|")))
          
          write.fasta(consensusSeq, names = getAnnot(consensusSeq), file.out = paste0("Step2_Sequence_", aa, "Consensus_cluster",c,".fa"))
        }}
      else {
        write.fasta(Fasta_RNA_col, names = paste0(aa, "_cluster",c," consensus sequence ", paste(gsub(organism_name, "", getName(Fasta_RNA_col)), collapse = "|")), file.out = paste0("Step2_Sequence_", aa, "_Unique_cluster",c,".fa"))
      }
    }
  } 
  else {if(length(Fasta_RNA_AA) != 0) {write.fasta(Fasta_RNA_AA, names = paste0(aa, "_cluster",c," consensus sequence ", paste(gsub(organism_name, "", getName(Fasta_RNA_AA)), collapse = "|")), file.out = paste0("Step2_Sequence_", aa, "_Unique_cluster",c,".fa"))}
}
}

#############################################################################################################
#PART 3
#TAKE xStep2_Uniques_sequences and xStep2_Consensus_sequences to merge in unique FASTA file

#MERGE OF SEQUENCE in one Fasta file take by aminoacid

out.file<-NULL
file.names <- dir(ResultsDir, pattern = paste0("^Step2_Sequence_"), full.names = T)
  
  for(i in 1:length(file.names)){
  
    file <- read.table(file.names[i],header=F, sep=";", stringsAsFactors=FALSE)
    out.file <- rbind(out.file, file)
  }

print (i)

setwd(extDataDir)
write.table(out.file, file = paste0("Merge_Step2_Collapsed_",organism_name, Sys.Date(), ".fa"),sep=";", col.names = F,
              row.names = FALSE, quote= F,fileEncoding="windows-1252")


