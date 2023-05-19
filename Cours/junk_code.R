source("usefull_func.R")

N = 10239
M = 447
K = 206
X = 33


mytable <- matrix(c(X, K-X, M-X, N-M-K+X), nr = 2) %>% as.data.frame()
colnames(mytable) <- c("Of_interest_yes", "Of_interest_yes")
rownames(mytable) <- c("in_pathway_yes","in_pathway_no")

mytable <- d
mytable <- rbind(mytable,colSums(d))
mytable <- cbind(mytable,c(rowSums(d),sum(d)))
colnames(mytable) <- c("Of_interest_yes", "Of_interest_no","Total")
rownames(mytable) <- c("in_pathway_yes","in_pathway_no","Total")

desired_colnames <- colnames(mytable) |> 
  str_remove('Of_interest_') |> 
  str_to_title()
names(desired_colnames) <- colnames(mytable)

desired_rownames <- rownames(mytable) |> 
  str_remove('in_pathway_') |> 
  str_to_title()
names(desired_rownames) <- rownames(mytable)

mytable %>%  gt %>% cols_label(.list = desired_colnames) |> 
  tab_spanner(
    label = md('**Genes of interest**'),
    columns = 1:2
  ) %>% rows_label(.list = desired_rownames) |> 
  tab_spanner(
    label = md('**Pathway of interest**'),
    rows = 1:2
  )

fisher.test(mytable,
            alternative = "greater")$p.value


phyper(X,M, N-M, K,lower.tail = F)


d <- data.frame(gene.not.interest=c(2613, n-k), gene.in.interest=c(28, 29))
row.names(d) <- c("In_category", "not_in_category")
d




## GSEA


GOonto <- GetONTOLOGY()
saveRDS("GOonto","data/GOonto.rds")

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

datatab.shuff <- 1:100 %>% parallel::mclapply(function(x){
  gsInfo(ego,my_go,shuffle=T)
},mc.cores=10) %>% bind_rows(.id="Shuffle")

p2 <- datatab.shuff %>% ggplot()   +
  # geom_segment(data=df2, aes(x=x, xend=x, y=0.05, yend=0),
  #              color="black") + 
  geom_line(aes(x = x,y=runningScore,col=Shuffle),linetype="dashed") +
  geom_line(data=datatab, aes(x = x,y=runningScore)) +
  ggtitle(my_desc) + ylab("Ranked List Metric") + xlab("Position in the Ranked List of Genes") + theme(legend.position="none")


selected.gs <-  ego@geneSets[my_go]

observed_info <- lapply(selected.gs, function(gs)
  gseaScores(geneSet=gs,
             geneList=ego@geneList,
             exponent=1)
)
observedScore <- sapply(observed_info, function(x) x$ES)

permScores <- parallel::mclapply(1:10000, function(i) {
  perm.gseaEScore(ego@geneList, selected.gs, 1)
},mc.cores=10)

permScores <- do.call("cbind", permScores)

rownames(permScores) <- names(selected.gs)

pos.m <- apply(permScores, 1, function(x) mean(x[x >= 0]))
neg.m <- apply(permScores, 1, function(x) abs(mean(x[x < 0])))


normalized_ES <- function(ES, pos.m, neg.m) {
  s <- sign(ES)
  m <- numeric(length(ES))
  m[s==1] <- pos.m[s==1]
  m[s==-1] <- neg.m[s==-1]
  ES/m
}

NES <- normalized_ES(observedScore, pos.m, neg.m)

permScores <- apply(permScores, 2, normalized_ES, pos.m=pos.m, neg.m=neg.m)


p5 <- enframe(permScores[1,]) %>% 
  ggplot(aes(y=value)) + geom_histogram(bins=50) + geom_hline(yintercept = ES,col="red",linetype="dashed")



### Code for cours


DE_result <- read_tsv("../TP/Gene_DE_DIVA.tsv")

DE_result

de_genes <- DE_result %>% dplyr::filter(p.adj < 0.05,logFC>0.5)

de_genes %>% dim

geneList <- DE_result %>% dplyr::select(rowname,logFC) %>% arrange(desc(logFC)) %>% deframe() 


## ClusterProfiler Analysis
require(clusterProfiler)
require(org.Hs.eg.db)
org.db <- org.Hs.eg.db

## We can also use our own DB
hub <- AnnotationHub()
org.db <- hub[['AH111575']]
org.db



### 



kk <- enrichKEGG(gene         = res2$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

kk %>% as_tibble()

res2 <- bitr(DE_result$rowname, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb=org.db)
geneListENTREZ <- DE_result %>% left_join(res2,by=c("rowname"="ENSEMBL")) %>%
  dplyr::select(ENTREZID,logFC) %>%
  arrange(desc(logFC)) %>% 
  drop_na() %>% 
  deframe()

kk2 <- gseKEGG(geneList     = geneListENTREZ,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

kk2 %>% as_tibble()

mkk <- enrichMKEGG(gene = res2$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)





## Universal enrichment analysis
TERM2GENE <- select(org.db, keys = keys(org.db,"GOALL"), column="ENSEMBL",keytype="GOALL") %>% dplyr::select(GOALL,ENSEMBL)
enricher(de_genes$rowname, TERM2GENE = TERM2GENE,universe = DE_result$rowname) %>% as_tibble()


genesets <- split(TERM2GENE,TERM2GENE$GOALL) %>% map("ENSEMBL")


hyp_obj <- hypeR::hypeR(de_genes$rowname, genesets)
