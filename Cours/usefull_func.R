perm.geneList <- function(geneList) {
  ## perm.idx <- sample(seq_along(geneList), length(geneList), replace=FALSE)
  perm.idx <- sample.int(length(geneList))
  perm.geneList <- geneList
  names(perm.geneList) <- names(geneList)[perm.idx]
  return(perm.geneList)
}
perm.gseaEScore <- function(geneList, geneSets, exponent=1) {
  geneList <- perm.geneList(geneList)
  res <- sapply(1:length(geneSets), function(i)
    gseaScores(geneSet=geneSets[[i]],
               geneList=geneList,
               exponent=exponent)$ES
  )
  return(res)
}


gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
  ###################################################################
  ##    geneList                                                   ##
  ##                                                               ##
  ## 1. Rank order the N genes in D to form L = { g_1, ... , g_N}  ##
  ##    according to the correlation, r(g_j)=r_j,                  ##
  ##    of their expression profiles with C.                       ##
  ##                                                               ##
  ###################################################################
  
  ###################################################################
  ##    exponent                                                   ##
  ##                                                               ##
  ## An exponent p to control the weight of the step.              ##
  ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
  ##   the standard Kolmogorov-Smirnov statistic.                  ##
  ##   When p = 1, we are weighting the genes in S                 ##
  ##   by their correlation with C normalized                      ##
  ##   by the sum of the correlations over all of the genes in S.  ##
  ##                                                               ##
  ###################################################################
  
  ## genes defined in geneSet should appear in geneList.
  ## this is a must, see https://github.com/GuangchuangYu/DOSE/issues/23
  geneSet <- intersect(geneSet, names(geneList))
  
  N <- length(geneList)
  Nh <- length(geneSet)
  
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet ## logical
  
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  
  Pmiss[!hits] <-  1/(N-Nh)
  Pmiss <- cumsum(Pmiss)
  
  runningES <- Phit - Pmiss
  
  ## ES is the maximum deviation from zero of Phit-Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if( abs(max.ES) > abs(min.ES) ) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  
  df <- data.frame(x=seq_along(runningES),
                   runningScore=runningES,
                   position=as.integer(hits)
  )
  
  if(fortify==TRUE) {
    return(df)
  }
  
  df$gene = names(geneList)
  res <- list(ES=ES, runningES = df)
  return(res)
}

gsInfo <- function(object, geneSetID,shuffle=F) {
  
  geneList <- object@geneList
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  if(shuffle){
    geneList <- perm.geneList(geneList)
  }
  geneSet <- object@geneSets[[geneSetID]]
  
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  } else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

GetONTOLOGY <- function(OrgDb=org.Hs.eg.db::org.Hs.eg.db,ont="BP",keytype = "ENSEMBL"){
  kt <- keytypes(OrgDb)
  
  kk <- keys(OrgDb, keytype=keytype)
  
  goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
  if (ont != "ALL") {
    goterms <- goterms[goterms == ont]
  }
  go2gene <- suppressMessages(
    AnnotationDbi::mapIds(OrgDb, keys=names(goterms), column=keytype,
                          keytype="GOALL", multiVals='list')
  )
  goAnno <- stack(go2gene)
  colnames(goAnno) <- c(keytype, "GOALL")
  goAnno <- unique(goAnno[!is.na(goAnno[,1]), ])
  goAnno$ONTOLOGYALL <- goterms[goAnno$GOALL]
  return(goAnno)
}
