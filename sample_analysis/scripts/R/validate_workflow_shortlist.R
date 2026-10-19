# Additional no-imputation refits, separate from the completed-matrix screen.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(Matrix))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),'common.R'))
opts <- parse_cli_options()
screen <- option_value(opts,'screen_dir',file.path(project_root(),'results/workflow_example_comparison'))
data_dir <- normalizePath(option_value(opts,'data_dir',detect_data_dir()))
dest <- file.path(screen,'shortlist_no_imputation')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
target_spec <- option_value(opts,'targets','SSC:G1:G5,BLA:G2:G4,SSC:G2:G4')
triplets <- strsplit(strsplit(target_spec,',',fixed=TRUE)[[1]],':',fixed=TRUE)
if(any(lengths(triplets)!=3L))stop('Use --targets REGION:GROUP_A:GROUP_B,REGION:GROUP_A:GROUP_B')
targets <- as.data.frame(do.call(rbind,triplets));colnames(targets) <- c('Region','GroupA','GroupB')
if(any(!targets$Region %in% REGION_CODES) || any(!targets$GroupA %in% paste0('G',1:6)) || any(!targets$GroupB %in% paste0('G',1:6)) || any(targets$GroupA==targets$GroupB))stop('Invalid validation target')
split_symbols <- function(x)unique(trimws(unlist(strsplit(x[!is.na(x)&nzchar(x)],';',fixed=TRUE))))
summaries <- list()
for(region in unique(targets$Region)) {
  prep <- preprocess_region_matrix(region,data_dir)
  # Retain only rows observed in every retained measurement, before any fill.
  keep <- rowSums(is.finite(prep$observed_counts))==ncol(prep$observed_counts)
  data <- prep$data[keep,];mat <- prep$normalized[keep,,drop=FALSE]
  stopifnot(all(is.finite(prep$observed_counts[keep,,drop=FALSE])))
  grp <- factor(prep$group_labels,levels=paste0('G',1:6));design <- model.matrix(~0+grp);colnames(design) <- levels(grp)
  fit0 <- lmFit(mat,design)
  symbols <- split_symbols(data$PG.Genes)
  mapping <- suppressMessages(select(org.Mm.eg.db,keys=symbols,keytype='SYMBOL',columns='ENTREZID'))
  mapping <- unique(mapping[!is.na(mapping$ENTREZID),c('SYMBOL','ENTREZID')]);lookup <- split(mapping$ENTREZID,mapping$SYMBOL)
  ids <- lapply(data$PG.Genes,function(x)unique(unlist(lookup[split_symbols(x)],use.names=FALSE)))
  unique_gene <- lengths(ids)==1L;row_gene <- rep(NA_character_,nrow(data));row_gene[unique_gene] <- vapply(ids[unique_gene],`[`,character(1),1)
  bg <- unique(row_gene[unique_gene])
  anno <- suppressMessages(select(org.Mm.eg.db,keys=bg,keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL')))
  anno <- unique(anno[!is.na(anno$GOALL),c('ENTREZID','GOALL','ONTOLOGYALL')])
  term <- sort(unique(anno$GOALL));ont <- anno$ONTOLOGYALL[match(term,anno$GOALL)]
  incidence <- sparseMatrix(i=match(anno$GOALL,term),j=match(anno$ENTREZID,bg),x=1,dims=c(length(term),length(bg)))
  sizes <- as.numeric(rowSums(incidence));ok <- sizes>=5 & sizes<=500
  incidence <- incidence[ok,,drop=FALSE];term <- term[ok];ont <- ont[ok];sizes <- sizes[ok]
  for(i in which(targets$Region==region)) {
    a <- targets$GroupA[i];b <- targets$GroupB[i];key <- paste(region,paste0(b,'_vs_',a),sep='_')
    fit <- eBayes(contrasts.fit(fit0,makeContrasts(contrasts=paste(b,a,sep='-'),levels=design)),trend=TRUE,robust=TRUE)
    p <- as.numeric(fit$p.value);q <- p.adjust(p,'BH');lfc <- as.numeric(fit$coefficients)
    selected <- q<.05 & abs(lfc)>log2(1.5);fg <- unique(row_gene[selected & unique_gene])
    counts <- as.numeric(incidence %*% as.integer(bg %in% fg));gp <- phyper(counts-1,sizes,length(bg)-sizes,length(fg),lower.tail=FALSE);gq <- p.adjust(gp,'BH')
    diff <- data.frame(ProteinAccessions=data$PG.ProteinAccessions,Genes=data$PG.Genes,Log2FC=lfc,PValue=p,AdjPValue=q,Selected=selected,UniqueEntrezID=row_gene)
    go <- data.frame(GOID=term,Ontology=ont,Count=counts,TermSize=sizes,ForegroundSize=length(fg),BackgroundSize=length(bg),PValue=gp,AdjPValue_BH=gq)
    write.csv(diff,file.path(dest,paste0(key,'_differential.csv')),row.names=FALSE,na='')
    write.csv(go,file.path(dest,paste0(key,'_GO.csv')),row.names=FALSE)
    main_con <- gzfile(file.path(screen,'differential',paste0(key,'.csv.gz')),'rt')
    main <- read.csv(main_con);close(main_con)
    overlap <- intersect(diff$ProteinAccessions[selected],main$ProteinAccessions[main$Selected])
    summary <- data.frame(Region=region,GroupA=a,GroupB=b,Tested=nrow(data),MainTested=prep$n_after_filter,
      Candidates=sum(selected),HighInB=sum(selected & lfc>0),HighInA=sum(selected & lfc<0),GO_Tested=length(term),GO_BH=sum(gq<.05),GO_MinBH=min(gq),
      Foreground=length(fg),Background=length(bg),CandidateOverlapWithMain=length(overlap),MainCandidates=sum(main$Selected),
      MainCandidateRetention=length(overlap)/sum(main$Selected),Strategy='No missing values in any retained measurement; same observed-value sample normalization; refit limma and rebuild GO background')
    summaries[[key]] <- summary
    cat(key,'complete-case candidates',sum(selected),'GO',sum(gq<.05),'\n');flush.console()
  }
}
updated <- do.call(rbind,summaries)
summary_file <- file.path(dest,'summary.csv')
if(file.exists(summary_file)) {
  old <- read.csv(summary_file)
  old_keys <- paste(old$Region,old$GroupA,old$GroupB)
  new_keys <- paste(updated$Region,updated$GroupA,updated$GroupB)
  updated <- rbind(old[!old_keys %in% new_keys,],updated)
}
write.csv(updated,summary_file,row.names=FALSE)
