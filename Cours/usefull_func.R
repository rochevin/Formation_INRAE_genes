
GetONTOLOGY <- function(OrgDb=org.Hs.eg.db::org.Hs.eg.db,ont="BP",keytpe = "ENSEMBL"){
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