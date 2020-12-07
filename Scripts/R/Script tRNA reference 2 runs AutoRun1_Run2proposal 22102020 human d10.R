######## tRNA reference creation by sequence clustering ##########
# version script YM 20-10-2020  run 1 in automatic mode
# THREE STEPS Step 1 run REMOVE DUPLICATES ONLY
#SECOND COLLAPS SEQUENCES
#STEP THREE CREATE SINGLE FILE
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)                            }
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)                     }
if (!require("seqinr")) {
  install.packages("seqinr", dependencies = TRUE)
  library(seqinr)
}
if (!require("msa")) {
  if (!require("BiocManager")){
    install.packages("BiocManager", dependencies = TRUE)
    library(BiocManager)
  }
  install(msa)
  library(msa)
}

rm(list=ls(all=TRUE))
#SET pdf size for heatmaps
pdf.width=15
pdf.height=15

#BLOCK READ FILE  COPY-PASTE from WIN KEYBOARD
###############################################
#COPY WINDOWS PATH TO CLIPBOARD
###############################################
extDataDir  <- gsub  ( "\\\\",  "/",  readClipboard ()  ) #READ FROM CLIPBOARD

#OR INDICATE extDataDir directly
#extDataDir <-  "/run/user/1001/gvfs/smb-share:server=193.54.30.197,share=runs/PROJECTS/2primeO/FlorianBioinfo/cytoplasmic_tRNAs/counting_distances_Bsubtilis" #READ AND send to extDataDir
setwd(extDataDir)

#GIVE NAME of Fasta FILE
#Fasta_file_RNApattern <-"stapAure_AUREUS_NCTC_832-mature-tRNAs_03062020"
Fasta_file_RNApattern <-"*..fa" #TAKE ANY *.fa file in folder

Fasta_file_RNA <- dir(paste0(extDataDir), pattern = paste0(Fasta_file_RNApattern), full.names = T)
Fasta_file_RNA #VERIFICATION of the file name

ResultsDir <- paste0(extDataDir,"/counting_distances_lim10") #10 for human
dir.create(ResultsDir, showWarnings = FALSE)  #CREATE FOLDER FOR tRNA analysis
setwd(ResultsDir)

set.seed(123)

AAlist <- list("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "SeC")

#READING FASTA file and transformation to DATAFRAME
# Maybe use "Biostrings" package? <-- Already done by library(msa)

## Functions ##
readinteger <- function(default) {
  n <- readline(prompt="Enter nc, the number of clusters : ")
  if(grepl("^[0-9]+$",n))
    return(as.integer(n))
  else
    return(default)
}

#readinteger <- function() { 
#  n <- readline(prompt="Enter nc, the number of clusters : ")
#  if(!grepl("^[0-9]+$",n))
#  {
#    return(readinteger())
#  }
#  return(as.integer(n))
#                          }

#EXTRACT AND CHECK ORGANISM NAME
Fasta_RNA <-read.fasta(file = Fasta_file_RNA)
head(Fasta_RNA)
RNAstring <- getName(Fasta_RNA)[1]
print(RNAstring)
organism_name <- unlist(strsplit(RNAstring,"tRNA"))[1]     # something like"Saccharomyces_cerevisiae_"  #TAKE FROM FASTA file
print(organism_name)

###########################################################################################################################
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
    
    #SAVE heatmap file in pdf   
    pdf(paste0("Step1_Heatmap2_", aa,".pdf"),width=pdf.width,height=pdf.height, paper='special') 
    heatmap.2(Dist_matrix, key = FALSE, cellnote = Dist_matrix, Rowv=T, Colv="Rowv", notecex=0.8, notecol="black", na.color=par("bg"), dendrogram="row", density.info="none",  labRow=getName(Fasta_RNA_AA), labCol=getName(Fasta_RNA_AA), margins = c(15, 15), trace=c("none"))
    dev.off()
    #DEFINE number of clusters nc as max number of unique distances in the 1st and last lines of the heatmap
    nc <-max(length(unique(Dist_matrix[1,])), length(unique(Dist_matrix[nrow(Dist_matrix),])))
    
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
                                  cluster = "upgma",              # matrice de distance --> arbre enracin?
                                  substitutionMatrix = "iub")     # prend en compte la notation IUPAC. + 1.9 si match, 0 si non.
        
        consensusSeq <- consensusString(consensusMatrix(alignment2))
        
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
   else {if (length(Fasta_RNA_AA) != 0) {
     write.fasta(Fasta_RNA_AA, names = paste(gsub(organism_name, "", getName(Fasta_RNA_AA)), collapse = "|"), file.out = paste0("Step2_Sequence_",aa, "_Unique_untouched__cluster",c,".fa"))
                                         }
  }
}

#####################################################################################################################
#READ AND CREATE NON-RUDUNDANT FASTA from all Step1_Sequence_NoDup_ files
#MERGE OF SEQUENCE in one Fasta file take by aminoacid

out.file<-NULL
file.names <- dir(ResultsDir, pattern = paste0("^Step1_Sequence_NoDup"), full.names = T) #TAKE all collapsed Step1 sequences
file.namesU <- dir(ResultsDir, pattern = paste0("Unique_untouched"), full.names = T) #TAKE Uniques sequences Step1
file.names <- append(file.names,file.namesU)

for(i in 1:length(file.names)){
  
  file <- read.table(file.names[i],header=F, sep=";", stringsAsFactors=FALSE)
  out.file <- rbind(out.file, file)
}

print (i)

setwd(extDataDir)
write.table(out.file, file = paste0("Merge_Step1_NoDup",organism_name,"_",i,"seq_", Sys.Date(),".fa"),sep=";", col.names = F,
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
    
    #DEFINE number of clusters nc as number distances >10 adaptation to human refs in all lines +1 (for Zero cluster)
   # nc <-1+ min(rowSums(Dist_matrix>8))
    
    suggestion <- 1+ min(rowSums(Dist_matrix>10))  
    
    print (paste0("tRNAs : ", aa,"  Suggested number of clusters : ", suggestion ,". Press ENTER to confirm or enter another number"))
    nc= readinteger(default = suggestion)
    
    
    #select the right number of clusters define nc value
    
   # print (paste0("tRNAs : ", aa,"  Suggested number of clusters : ", nc ,"  confirm or enter another number"))
   # nc= readinteger()
    
    graphics.off()
   
    ## Clustering ##
    d<-as.dist(Dist_matrix) #convert to distances
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
                                  cluster = "upgma",              # matrice de distance --> arbre enracin?
                                  substitutionMatrix = "iub")     # prend en compte la notation IUPAC. + 1.9 si match, 0 si non.
        
        consensusSeq <- consensusString(consensusMatrix(alignment2))
        consensusSeq
        
        if (countPattern(pattern = "?", subject = c2s(getSequence(consensusSeq))) >= 12) {   #12 adattation for human ref
          write.fasta(Fasta_RNA_col, names = getName(Fasta_RNA_col), file.out = paste0("Step2_Removed_sequences_", aa, "_cluster",c,".fa"))
        }
        else {
          consensusSeq <- gsub(pattern = "[?]", "N", consensusSeq)
          consensusSeq <- gsub(pattern = "[-]", "", consensusSeq)
          consensusSeq <- as.SeqFastadna(s2c(consensusSeq), Annot = paste0(aa, "_cluster",c," consensus sequence ", paste(gsub(organism_name, "", getName(Fasta_RNA_col)), collapse = "|")))
          
          write.fasta(consensusSeq, names = getAnnot(consensusSeq), file.out = paste0("Step2_Sequence_", aa, "_Consensus_cluster",c,".fa"))
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
write.table(out.file, file = paste0("Merge_Step2_Collapsed_",organism_name,"_",i,"seq_", Sys.Date(), ".fa"),sep=";", col.names = F,
              row.names = FALSE, quote= F,fileEncoding="windows-1252")


##########################################################################################################
#LAST GLOBAL HEATMAP for ALL sequences (to check cross aminoacid distances)
#HEATMAP for ALL sequences in Merge_Step2_Collapsed_ file
Fasta_file_RNA <- list.files(extDataDir, pattern="Merge_Step2_Collapsed_", recursive=FALSE, full.names = TRUE)
Fasta_RNA <-read.fasta(file = Fasta_file_RNA)
head(Fasta_RNA)

pattern <- character()
for (i in 1:length(Fasta_RNA)){pattern <- c(pattern, c2s(Fasta_RNA[[i]]))}

Dist_matrix <- matrix(nrow = length(Fasta_RNA), ncol = length(Fasta_RNA)) #EMPTY DISTANCE MATRIX
#FILL DISTANCE MATRIX by calculation of distances type = "overlap", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE), gapOpening = 1, gapExtension = 0
for (i in 1:length(Fasta_RNA)){
  alignment <- pairwiseAlignment(DNAStringSet(pattern),DNAStringSet(c2s(getSequence(Fasta_RNA[[i]]))), scoreOnly = TRUE, type = "overlap", substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE), gapOpening = 1, gapExtension = 0) 
  Dist_matrix[,i] <- alignment
  Dist_matrix[,i] <- getLength(Fasta_RNA[[i]]) - Dist_matrix[,i]
}
#Dist_matrix may be assymetric!
print(Dist_matrix)

#Order matrix by heatmap
heatmap.2(Dist_matrix, key = FALSE, cellnote = Dist_matrix, Rowv=F, Colv="Rowv", notecex=0.8, notecol="black", na.color=par("bg"), dendrogram="none", density.info="none",  labRow=getName(Fasta_RNA), labCol=getName(Fasta_RNA), margins = c(12, 12), trace=c("none"))

pdf.width=25
pdf.height=25
#SAVE file in pdf   
pdf(paste0("Step2_FinalFullHeatmap_",organism_name,".pdf"),width=pdf.width,height=pdf.height, paper='special') 
heatmap.2(Dist_matrix, key = FALSE, cellnote = Dist_matrix, Rowv=F, Colv="Rowv", notecex=0.8, notecol="black", na.color=par("bg"), dendrogram="none", density.info="none",  labRow=getName(Fasta_RNA), labCol=getName(Fasta_RNA), margins = c(10, 10), trace=c("none"))
dev.off()

###########################################################################################################
###########   REANNOTATION PART   #########################################################################

Fasta_file_RNA <- list.files(extDataDir, pattern="Merge_Step2_Collapsed_", recursive=FALSE, full.names = TRUE) #USE MERGED file Step2
Fasta_RNA <-read.fasta(file = Fasta_file_RNA)

Fasta_RNA_txt <-read.delim(file = Fasta_file_RNA, sep=" ", header= FALSE) #READ FASTA file in TEXT form, col V4 contains string to treat

Fasta_RNA_txt1 <-Fasta_RNA_txt[grep('>tRNA',Fasta_RNA_txt$V1),] #select lines format tRNA unique
Fasta_RNA_txt1$aa <-substr(Fasta_RNA_txt1$V1,7,9)
Fasta_RNA_txt1$ac <-substr(Fasta_RNA_txt1$V1,11,13)

Fasta_RNA_txt2 <-Fasta_RNA_txt[grep('>..._cluster',Fasta_RNA_txt$V1),] #select other
Fasta_RNA_txt3 <- Fasta_RNA_txt2[grep("[f|i]Met",Fasta_RNA_txt2$V4 ),] #select lines for tRNAfMet
Fasta_RNA_txt2 <- Fasta_RNA_txt2[-grep("[f|i]Met",Fasta_RNA_txt2$V4 ),] # take remaining lines
Fasta_RNA_txt2$aa<-substr(Fasta_RNA_txt2$V1,2,4) #extract 

#Merge anticodons in one string separator "_", check if this is goog for future
for(i in 1:nrow(Fasta_RNA_txt2)) {Fasta_RNA_txt2$ac[i]<-paste(unique(substr(unlist(strsplit(gsub("[^a-zA-Z]", "", as.character(Fasta_RNA_txt2$V4[i])), 'tRNA'))[-1], 4,6)),  collapse="_")
}
for(i in 1:nrow(Fasta_RNA_txt3)) {Fasta_RNA_txt3$ac[i]<-paste(unique(substr(unlist(strsplit(gsub("[^a-zA-Z]", "", as.character(Fasta_RNA_txt3$V4[i])), 'tRNA'))[-1], 5,7)),  collapse="_")
}
Fasta_RNA_txt3$aa<- "iMet"
#MERGE by rbind
Fasta_RNA_Annotated <-rbind(Fasta_RNA_txt1, Fasta_RNA_txt2, Fasta_RNA_txt3)
Fasta_RNA_Annotated$Name <- paste0(">zz_tRNA",Fasta_RNA_Annotated$aa,"_",Fasta_RNA_Annotated$ac) #CREATE tRNA name
Fasta_RNA_Annotated$iso <-" "

for (i in unique(Fasta_RNA_Annotated$Name)) {
  n<- nrow(Fasta_RNA_Annotated[grep(i,Fasta_RNA_Annotated$Name),])
  if (n>1) {Fasta_RNA_Annotated[grep(i,Fasta_RNA_Annotated$Name),]$iso<-c(1:n)}
  
}

Fasta_RNA_Annotated$Name<-paste0(Fasta_RNA_Annotated$Name, Fasta_RNA_Annotated$iso)

#REPLACE NAMES take column V1 from Fasta_RNA_Annotated, find equivalent in Fasta_RNA (object) and replace by  Fasta_RNA_Annotated$Name
output <- Fasta_RNA # To not erase Fasta_RNA by accident 

getName(output)
Fasta_RNA_Annotated$Name
Fasta_RNA_Annotated$V1

for (i in getName(output)){
  attr(output[[grep(i, attr(output, "name"))]], "name") <- gsub(">","",Fasta_RNA_Annotated$Name[grep(i, Fasta_RNA_Annotated$V1)])
}

write.fasta(sequences = getSequence(output), names = getName(output), file.out = paste0("Merge_Step2_Collapsed_and_Reannotated_",organism_name,"_",length(output),"seq_", Sys.Date(), ".fa"))
