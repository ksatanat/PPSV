#!/usr/bin/env R


# PPSV features
setwd('Path')
require(readxl)
require(readr)
require(igraph)
require(protr)
require(Biostrings)
require(testit)

data("BLOSUM62") # Biostrings
## Required data; drug-TP (DB), drug-disease (CTD), DP-disease (disGenet)
# Read drug (db) - TP (unp) ; DrugBank
db_tp <- read.table("./drug_target.txt", header = T, stringsAsFactors = F)
#head(db_tp)
# Read proteins' sequences
seq_df <- read.table('./sequence_unp.txt', header = T, stringsAsFactors = F)
# Read PPI data
string <- read.table('./ppi_data.txt', header = T, stringsAsFactors = F)
#head(string)
string800 <- string[which(string$combined_score >= 800), c(2,4,5)] # Choose only confident score >= 80%
#nrow(string800)
net<-graph.data.frame(string800,directed = F) # Generate network with 800 confident score
# Read metabolic pathways  
kegg_unp  <- read.table('./pathway_data.txt', header = T)
#head(kegg_unp)




# Nei score
nei_score <- function(net, p1, p2, margin){
  df <- as.data.frame(cbind(as.vector(net[p1,]), as.vector(net[p2,])))
  M_11 <- apply(df, margin, sum) == 2
  return(sum(M_11)/sqrt(sum(df[,1]) * sum(df[,2])))
}

# Closer score
closer_score<- function(net,p1,p2){
  shortstep <- suppressWarnings({length(as_ids(get.shortest.paths(net,from = as.character(p1), to = as.character(p2))$`vpath`[[1]]))-1})
  if (shortstep > 0) { 
    return(1/shortstep)
  }  else {
    return(0) # case of two disjoint proteins
  }
}

# Loc & Glo scores
alinment <- function(p1,p2){
  seq1 <- seq_df$seq[seq_df$tp_dp == p1]#getUniProt(p1)[[1]]
  seq2 <- seq_df$seq[seq_df$tp_dp == p2]#getUniProt(p2)[[1]]
  if (seq1 =='NA' | seq2 =='NA' | length(seq1)==0 | length(seq2)==0 |
      has_error(pairwiseAlignment(seq1, seq2
                                  , substitutionMatrix = "BLOSUM62",
                                  gapOpening = 10, gapExtension = 0.5, type = "local", scoreOnly =T),silent = TRUE)) {
    return(rep(0,2))
  } else {
    local_aln <- pairwiseAlignment(seq1, seq2
                                   , substitutionMatrix = "BLOSUM62",
                                   gapOpening = 10, gapExtension = 0.5, type = "local", scoreOnly =T)
    global_aln <- pairwiseAlignment(seq1, seq2
                                    , substitutionMatrix = "BLOSUM62",
                                    gapOpening = 10, gapExtension = 0.5, type = "global", scoreOnly =T)
    return(c(local_aln, global_aln))
  }
}

mainScores <- function(p1, p2) {
  feature_score <- vector()
  
  # Nei & Closer score
  if (p1 %in% as_ids(V(net)) & p2 %in% as_ids(V(net))) {
    if (p1 == p2) {
      feature_score <- c(feature_score, rep(1,2)) # score of self-proteins
    } else {
      feature_score <- c(feature_score, nei_score(net, p1, p2, 1)) 
      feature_score <- c(feature_score, closer_score(net, p1, p2))
    }
  } else {
    feature_score <- c(feature_score, rep(0,2)) 
  }
  # Conf score
  if (p1 == p2) {
    feature_score <- c(feature_score, 1000) # score of self-proteins
  } else {
    string_conf <- string[string$unp_p1 == p1 & string$unp_p2 == p2,]
    if (nrow(string_conf) !=0) {
      feature_score <- c(feature_score, string_conf$combined_score[1])
    } else {
      feature_score <- c(feature_score, 0)
    }
  }
  # Loc  & Glo scores
  feature_score <- c(feature_score, alinment(p1,p2))
  # Pathway score
  func_score <- length(intersect(kegg_unp$Entry[which(kegg_unp$Uniprot_id == p1)], kegg_unp$V1[which(kegg_unp$Uniprot_id == p2)]))
  feature_score <- c(feature_score, func_score)
  # ShareDr score
  sharedrug <- length(intersect(db_tp$DrugBankIDs[grepl(p1, db_tp$Target.uniprot)], db_tp$DrugBankIDs[grepl(p2, db_tp$Target.uniprot)]))
  feature_score <- c(feature_score, sharedrug)
  return(paste0(feature_score, collapse = ','))
}

p1 = 'P00734'
p2 = 'P11245'
print(paste('The PPSV scores:', mainScores(p1, p2)))

