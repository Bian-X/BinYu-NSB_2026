#!/usr/bin/env Rscript
# Optional exhaustive, exploratory comparison; does not overwrite chosen examples.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(Matrix))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),'common.R'))
opts <- parse_cli_options()
data_dir <- normalizePath(option_value(opts,'data_dir',detect_data_dir()))
dest <- option_value(opts,'output_dir',file.path(project_root(),'results/workflow_example_comparison'))
reference_dir <- option_value(opts,'reference_dir',file.path(project_root(),'results'))
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
for(folder in c('differential','GO','provenance'))dir.create(file.path(dest,folder),showWarnings=FALSE)
pairs <- combn(paste0('G',1:6),2,simplify=FALSE)
contrasts <- vapply(pairs,function(p)paste(p[2],p[1],sep='-'),character(1))
summaries <- list(); all_go <- list(); all_p <- list(); validation <- list()
write_gzip <- function(table,path) {
  con <- gzfile(path,'wt');on.exit(close(con));write.csv(table,con,row.names=FALSE,na='')
}
split_symbols <- function(x)unique(trimws(unlist(strsplit(x[!is.na(x)&nzchar(x)],';',fixed=TRUE))))
for(region in REGION_CODES) {
  cat(format(Sys.time()),'preprocessing',region,'\n');flush.console()
  prep <- preprocess_region_matrix(region,data_dir)
  data <- prep$data; n <- nrow(data)
  grp <- factor(prep$group_labels,levels=paste0('G',1:6))
  design <- model.matrix(~0+grp);colnames(design) <- levels(grp)
  fit <- eBayes(contrasts.fit(lmFit(prep$normalized,design),makeContrasts(contrasts=contrasts,levels=design)),trend=TRUE,robust=TRUE)
  symbols <- split_symbols(data$PG.Genes)
  mapping <- suppressMessages(select(org.Mm.eg.db,keys=symbols,keytype='SYMBOL',columns='ENTREZID'))
  mapping <- unique(mapping[!is.na(mapping$ENTREZID),c('SYMBOL','ENTREZID')])
  lookup <- split(mapping$ENTREZID,mapping$SYMBOL)
  ids <- lapply(data$PG.Genes,function(x)unique(unlist(lookup[split_symbols(x)],use.names=FALSE)))
  unique_gene <- lengths(ids)==1L
  row_gene <- rep(NA_character_,n);row_gene[unique_gene] <- vapply(ids[unique_gene],`[`,character(1),1)
  bg <- unique(row_gene[unique_gene])
  anno <- suppressMessages(select(org.Mm.eg.db,keys=bg,keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL')))
  anno <- unique(anno[!is.na(anno$GOALL),c('ENTREZID','GOALL','ONTOLOGYALL')])
  term_ids <- sort(unique(anno$GOALL));ont <- anno$ONTOLOGYALL[match(term_ids,anno$GOALL)]
  membership <- sparseMatrix(i=match(anno$GOALL,term_ids),j=match(anno$ENTREZID,bg),x=1,dims=c(length(term_ids),length(bg)))
  sizes <- as.numeric(rowSums(membership));eligible <- sizes>=5 & sizes<=500
  membership <- membership[eligible,,drop=FALSE];sizes <- sizes[eligible];term_ids <- term_ids[eligible];ont <- ont[eligible]
  go_test <- function(foreground) {
    hits <- as.numeric(membership %*% as.integer(bg %in% foreground))
    p <- phyper(hits-1,sizes,length(bg)-sizes,length(foreground),lower.tail=FALSE)
    list(count=hits,p=p,q=p.adjust(p,'BH'))
  }
  write.csv(data.frame(ProteinAccessions=data$PG.ProteinAccessions,Genes=data$PG.Genes,
    MappedGeneCount=lengths(ids),UniqueEntrezID=row_gene),file.path(dest,'provenance',paste0(region,'_unique_mapping.csv')),row.names=FALSE,na='')
  write.csv(prep$manifest,file.path(dest,'provenance',paste0(region,'_sample_manifest.csv')),row.names=FALSE)
  for(j in seq_along(pairs)) {
    a <- pairs[[j]][1]; b <- pairs[[j]][2];key <- paste(region,paste0(b,'_vs_',a),sep='_')
    ca <- which(prep$group_labels==a);cb <- which(prep$group_labels==b)
    oa <- rowSums(is.finite(prep$observed_counts[,ca,drop=FALSE]));ob <- rowSums(is.finite(prep$observed_counts[,cb,drop=FALSE]))
    lfc <- as.numeric(fit$coefficients[,j]);p <- as.numeric(fit$p.value[,j]);q <- p.adjust(p,'BH')
    selected <- q<.05 & abs(lfc)>log2(1.5)
    complete <- oa==length(ca) & ob==length(cb)
    fg <- unique(row_gene[selected & unique_gene]);fg_complete <- unique(row_gene[selected & unique_gene & complete])
    go <- go_test(fg);go_complete <- go_test(fg_complete)
    diff <- data.frame(Region=region,GroupA=a,GroupB=b,ProteinAccessions=data$PG.ProteinAccessions,Genes=data$PG.Genes,
      Log2FC=lfc,ModeratedT=as.numeric(fit$t[,j]),PValue=p,AdjPValue=q,Observed_A=oa,Observed_B=ob,Samples_A=length(ca),Samples_B=length(cb),
      Selected=selected,UniqueEntrezID=row_gene,FullyObservedContrast=complete,MeanLog2_A=rowMeans(prep$normalized[,ca,drop=FALSE]),MeanLog2_B=rowMeans(prep$normalized[,cb,drop=FALSE]))
    write_gzip(diff,file.path(dest,'differential',paste0(key,'.csv.gz')))
    tab <- data.frame(Region=region,GroupA=a,GroupB=b,GOID=term_ids,Ontology=ont,Count=go$count,TermSize=sizes,
      ForegroundSize=length(fg),BackgroundSize=length(bg),PValue=go$p,AdjPValue_BH=go$q,
      FullyObservedCount=go_complete$count,FullyObservedForeground=length(fg_complete),FullyObservedP=go_complete$p,FullyObservedBH=go_complete$q)
    all_go[[key]] <- tab;all_p[[key]] <- p
    summaries[[key]] <- data.frame(Region=region,GroupA=a,GroupB=b,Contrast=paste(b,a,sep='-'),Tested=n,Samples_A=length(ca),Samples_B=length(cb),
      Candidates=sum(selected),HighInB=sum(selected & lfc>0),HighInA=sum(selected & lfc<0),FullyObservedCandidates=sum(selected & complete),
      FullyObservedFraction=if(sum(selected)) mean(complete[selected]) else NA_real_,
      WeakCoverageCandidates=sum(selected & (oa<3 | ob<3)),UniqueForeground=length(fg),UniqueBackground=length(bg),
      GO_Tested=length(term_ids),GO_BH=sum(go$q<.05),GO_MinBH=min(go$q),
      GO_FullyObserved_BH=sum(go_complete$q<.05),GO_Both_BH=sum(go$q<.05 & go_complete$q<.05),
      MedianAbsLFC=if(sum(selected)) median(abs(lfc[selected])) else NA_real_)
    # Verify the screen against the primary PL contrast and GO results.
    if(region=='PL' & a=='G1' & b=='G5') {
      old <- read.csv(file.path(reference_dir,'PL_G5_vs_G1_differential_results_BH_FC1p5.csv'))
      idx <- match(diff$ProteinAccessions,old$ProteinAccessions)
      delta <- max(abs(diff$PValue-old$PValue[idx]));delta_fc <- max(abs(diff$Log2FC-old$Log2FC[idx]))
      oldgo <- read.csv(file.path(reference_dir,'PL_G5_vs_G1_GO_enrichment_results_BH_FC1p5.csv'))
      go_idx <- match(tab$GOID,oldgo$GOID)
      delta_go <- max(abs(tab$PValue-oldgo$PValue[go_idx]))
      stopifnot(!anyNA(idx),!anyNA(go_idx),delta<1e-12,delta_fc<1e-12,delta_go<1e-12,sum(selected)==251,length(fg)==242)
      validation <- data.frame(Reference='PL_G5_vs_G1',MaxPError=delta,MaxLog2FCError=delta_fc,MaxGO_PError=delta_go)
    }
    cat('  ',key,'DE',sum(selected),'GO',sum(go$q<.05),'complete-foreground GO',sum(go_complete$q<.05),'\n');flush.console()
  }
  # Save current progress even if a later region fails.
  write.csv(do.call(rbind,summaries),file.path(dest,'comparison_summary_partial.csv'),row.names=FALSE)
}
summary <- do.call(rbind,summaries)
# An extra exploratory correction over the whole search, in addition to each local family.
raw_go <- unlist(lapply(all_go,`[[`,'PValue'),use.names=FALSE)
raw_complete_go <- unlist(lapply(all_go,`[[`,'FullyObservedP'),use.names=FALSE)
global_go <- p.adjust(raw_go,'BH');global_complete_go <- p.adjust(raw_complete_go,'BH')
global_diff <- p.adjust(unlist(all_p,use.names=FALSE),'BH')
offset <- 0L;offset_diff <- 0L
for(key in names(all_go)) {
  tab <- all_go[[key]];count <- nrow(tab);idx <- offset+seq_len(count)
  tab$AcrossSearchBH <- global_go[idx];tab$FullyObservedAcrossSearchBH <- global_complete_go[idx]
  write_gzip(tab,file.path(dest,'GO',paste0(key,'.csv.gz')))
  summary[key,'GO_AcrossSearch_BH'] <- sum(tab$AcrossSearchBH<.05)
  summary[key,'GO_StableAcrossSearch_BH'] <- sum(tab$AcrossSearchBH<.05 & tab$FullyObservedAcrossSearchBH<.05)
  summary[key,'MinGO_AcrossSearchBH'] <- min(tab$AcrossSearchBH)
  nd <- length(all_p[[key]]);summary[key,'DifferentialAcrossSearch_BH'] <- sum(global_diff[offset_diff+seq_len(nd)]<.05)
  offset <- offset+count;offset_diff <- offset_diff+nd
}
summary$SearchComparisons <- length(pairs)*length(REGION_CODES)
write.csv(summary,file.path(dest,'comparison_summary.csv'),row.names=FALSE)
dictionary <- suppressMessages(select(GO.db::GO.db,keys=unique(unlist(lapply(all_go,`[[`,'GOID'))),keytype='GOID',columns=c('TERM','ONTOLOGY')))
write.csv(unique(dictionary),file.path(dest,'GO_term_dictionary.csv'),row.names=FALSE)
write.csv(validation,file.path(dest,'reference_validation.csv'),row.names=FALSE)
write.csv(data.frame(Comparisons=nrow(summary),DifferentialTests=length(global_diff),GOTests=length(global_go),
  Mapping='unique Entrez IDs; deduplicated',Family='all eligible BP/CC/MF 5-500 terms including zero overlap',
  Imputation='protein row median; normalized log2',Selection='within-analysis BH<0.05 and absolute log2FC>log2(1.5)',
  MissingnessSensitivity='same-background GO using selected genes from protein groups fully observed in both contrast groups; not a no-imputation model refit',
  SearchStatus='Exploratory exhaustive screen, not preregistered validation'),file.path(dest,'screen_metadata.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dest,'provenance/R_sessionInfo.txt'))
cat('Completed',nrow(summary),'comparisons;',length(global_go),'GO tests\n')
