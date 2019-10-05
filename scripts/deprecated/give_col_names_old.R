#Give column names old

givecolumnnames2 <- function(){
  sim_names <- c("N", "px", "iter", "nrepl", "y1", "y2", "y3")
  trueval_names <- c("Alfa1", "Alfa2", paste0("Beta", 1:px))
  em_names <- c("EMa1", "EMa2", paste0("EMb", 1:px))
  bi_names <- c("BIa1", paste0("BIb", 1:px))
  tri_names <- c("Tri_a1", "Tri_a2", paste0("Trib", 1:px))
  se_names <- c("bia1sd", "bib1sd", "ema1sd", "ema2sd", "emb1sd")
  columnnames <- c(sim_names, trueval_names, em_names, bi_names, tri_names, se_names)
  return(columnnames)
}

colnames(dfs2) <- givecolumnnames2()