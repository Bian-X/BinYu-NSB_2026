packages <- c('limma','AnnotationDbi','org.Mm.eg.db','GO.db','statmod','Matrix')
cat('R\t',as.character(getRversion()),'\n',sep='')
for(p in packages)cat(p,'\t',if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else'MISSING','\n',sep='')
