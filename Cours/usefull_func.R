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
    geneSet <- sample(names(geneList),size = length(object@geneSets[[geneSetID]]),replace = F)
  }else{
    geneSet <- object@geneSets[[geneSetID]]
  }
  
  
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


N = 10239
M = 447
K = 206
X = 33


mytable <- matrix(c(X, K-X, M-X, N-M-K+X), nr = 2) %>% as.data.frame()
colnames(mytable) <- c("Of_interest_yes", "Of_interest_yes")
rownames(mytable) <- c("in_pathway_yes","in_pathway_no")


desired_colnames <- colnames(mytable) |> 
  str_remove('Of_interest_') |> 
  str_to_title()
names(desired_colnames) <- actual_colnames

cols_label(.list = desired_colnames) |> 
  tab_spanner(
    label = md('**Adelie**'),
    columns = 3:4
  ) |> 
  tab_spanner(
    label = md('**Chinstrap**'),
    columns = c('Chinstrap_female', 'Chinstrap_male')
  )

fisher.test(mytable,
            alternative = "greater")$p.value


phyper(X,M, N-M, K,lower.tail = F)


d <- data.frame(gene.not.interest=c(2613, n-k), gene.in.interest=c(28, 29))
row.names(d) <- c("In_category", "not_in_category")
d



## GSEA


GOonto <- GetONTOLOGY()


geneList <- DE_result %>% dplyr::select(rowname,logFC) %>% arrange(desc(logFC)) %>% deframe() 

ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Hs.eg.db,
             ont          = "BP",
             keyType = "ENSEMBL",
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE)


my_go = "GO:0001944"
my_desc <- AnnotationDbi::Term(GO.db::GOTERM)[[my_go]]



gseaplot(ego, geneSetID = my_go, title = my_desc)

datatab <- gsInfo(ego,my_go)

p <- ggplot(datatab, aes_(x = ~x)) 
df2 <- data.frame(x = which(p$data$position == 1))
df2$y <- p$data$geneList[df2$x]
p <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                 color="black") + ggtitle(my_desc) + ylab("Ranked List Metric") + xlab("Position in the Ranked List of Genes")

p1 <- ggplot(datatab, aes(x = x,y=runningScore))   +
  geom_segment(data=df2, aes(x=x, xend=x, y=0.05, yend=0),
               color="black") + 
  geom_line() +
  ggtitle(my_desc) + ylab("Ranked List Metric") + xlab("Position in the Ranked List of Genes")


## Shuffle a bit

datatab.shuff <- 1:1000 %>% map(function(x){
  gsInfo(ego,my_go,shuffle=T)
}) %>% bind_rows(.id="Shuffle")

p2 <- ggplot(datatab.shuff)   +
  # geom_segment(data=df2, aes(x=x, xend=x, y=0.05, yend=0),
  #              color="black") + 
  geom_line(aes(x = x,y=runningScore,col=Shuffle),linetype="dashed") +
  geom_line(data=datatab, aes(x = x,y=runningScore)) +
  ggtitle(my_desc) + ylab("Ranked List Metric") + xlab("Position in the Ranked List of Genes") + theme(legend.position="none")


p5 <- datatab.shuff %>% group_by(Shuffle) %>% summarise(runningScore=max(runningScore)) %>% 
  ggplot(aes(y=runningScore)) + geom_histogram(bins=50) + geom_hline(yintercept = max(datatab$runningScore),col="red",linetype="dashed")
