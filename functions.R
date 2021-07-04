run_CellAssign <- function(ace, marker) {

  # Create sce object
  ace <- computeSumFactors(ace)
  sf <- sizeFactors(ace)
  
  mm = ACTIONet:::.preprocess_annotation_markers(lapply(markers, toupper), rownames(ace))
  mm = mm[fast_row_sums(mm)>0, ]
  
  # cellassign
  cellassign.fit <- cellassign(exprs_obj = ace[rownames(mm),], 
                                     marker_gene_info = as.matrix(mm), 
                                     s = sf, 
                                     learning_rate = 1e-2, 
                                     shrinkage = TRUE,
                                     verbose = TRUE)
  
  Labels <- celltypes(cellassign.fit)
  #scores = sapply(split(1:ncol(ace), Labels), function(idx) as.numeric(sparseVector(1, idx, ncol(ace))))
  scores = cellprobs(cellassign.fit)
  
  out = list(Labels = results$cell_labels, scores = scores)
  
  return(out)
}


run_Garnett <- function(ace, marker, species = "human") {
  library(garnett)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  
  # Construct CDS object
  suppressWarnings( {cds <- newCellDataSet(as(counts(ace), "dgCMatrix"), phenoData = new("AnnotatedDataFrame", data = as.data.frame(colData(ace))), featureData = new("AnnotatedDataFrame", data = as.data.frame(rowData(ace))))} )
  cds <- estimateSizeFactors(cds)
  
  # Export markers in a format Garnett likes
  lines = unlist(lapply(names(markers), function(n) sprintf(">%s\nexpressed: %s\n", n, paste(markers[[n]], collapse = ", "))))
  con <- file("markers.txt", "w")
  writeLines(lines, con)
  close(con)
  
  
  if(species == "human") {
    org.db = org.Hs.eg.db
  } else {
    org.db = org.Mm.eg.db
  }
  
  # marker_check <- check_markers(cds, "markers.txt",
  #                               db=org.db,
  #                               cds_gene_id_type = "SYMBOL",
  #                               marker_file_gene_id_type = "SYMBOL")
  
  
  classifier <- train_cell_classifier(cds = cds,
                                      marker_file = "markers.txt",
                                      db = org.db,
                                      cds_gene_id_type = "SYMBOL",
                                      num_unknown = 50,
                                      marker_file_gene_id_type = "SYMBOL")
  cds.out <- classify_cells(cds, classifier,
                            db = org.db,
                            cluster_extend = TRUE,
                            cds_gene_id_type = "SYMBOL")
  
  Labels = pData(cds.out)$cluster_ext_type
  scores = sapply(split(1:ncol(ace), Labels), function(idx) as.numeric(sparseVector(1, idx, ncol(ace))))
  
  out = list(Labels = Labels, scores = scores)
  
  return(out)  
}




run_SCINA <- function(ace, marker) {
  S = logcounts(ace)
  
  start_time <- Sys.time()
  results = SCINA(S, markers)
  end_time <- Sys.time()
  
  out = list(Labels = results$cell_labels, scores = Matrix::t(results$probabilities))
  
  return(out)
}

run_VISION <- function(ace, markers, thread_no = 6) {
  library(VISION)
  VISION.sigs = lapply(names(markers), function(nn) {
    gg = markers[[nn]]
    v = rep(1, length(gg))
    names(v) = gg
    ss = createGeneSignature(nn, v)
    return(ss)
  })
  
  vis <- Vision(ace, signatures = VISION.sigs)
  options(mc.cores = thread_no)
  vis <- analyze(vis)
  
  sigScores <- getSignatureScores(vis)
  Labels = names(markers)[apply(sigScores, 1, which.max)]
  res = list(Labels = Labels, scores = sigScores, obj = vis)
  
  return(res)
}

assess.annotation <- function(ace, scores, labels, true_labels, figures.path = "results/figures", exp_name = "eval") {
  if(!is.factor(labels)) {
    labels = as.factor(labels)
  }
  if(!is.factor(true_labels)) {
    true_labels = as.factor(true_labels)
  }
  
  gg0 = plot.ACTIONet(ace, true_labels)
  gps = lapply(colnames(scores), function(nn) {
    x = scores[, nn]
    gp = plot.ACTIONet.gradient(ace, x) + ggtitle(nn)
  })
  png(file.path(figures.path, sprintf("%s_scores.png", exp_name)), width = 16, height = 16, units = "in", res = 300)
  cowplot::plot_grid(plotlist = c(list(gg0), gps))
  dev.off()
  
  
  wilcox.out = presto::wilcoxauc(t(scores), true_labels)
  ii = match(wilcox.out$feature, colnames(scores))
  jj = match(wilcox.out$group, levels(true_labels))
  xx = wilcox.out$auc - 0.5
  enrichment = as.matrix(sparseMatrix(i = ii, j = jj, x = xx, dims = c(ncol(scores), length(levels(true_labels)))))
  rownames(enrichment) = colnames(scores)
  colnames(enrichment) = levels(true_labels)
  
  pdf(file.path(figures.path, sprintf("%s_alignment_heatmap.pdf", exp_name)), width = 12, height = 12)
  plot_pairwise_alignment(enrichment, rownames(enrichment), colnames(enrichment))
  dev.off()
  
  W = enrichment + 0.5
  meanAUC = sum(MWM_hungarian(W))/max(min(nrow(W), ncol(W)))
  
  NMI = ClusterR::external_validation(as.numeric(factor(labels)), as.numeric(factor(true_labels)), "nmi")
  ARI = ClusterR::external_validation(as.numeric(factor(labels)), as.numeric(factor(true_labels)), "adjusted_rand_index")
  
  out = c(ARI = ARI, NMI = NMI, meanAUC = meanAUC)
  
  return(out)
}








