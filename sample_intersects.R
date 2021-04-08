cytokines.data = read.csv('Re__Obesity_data_set_2/overfeeding_Cytokine_Chemokine_luminexWithSSPG.csv', row.names = 1,
                         header = F, stringsAsFactors = FALSE, sep=',', dec='.')

meta.data = read.csv('Re__Obesity_data_set/overfeeding_metabolomics.csv', row.names = 1,
                     header = F, stringsAsFactors = FALSE, sep=',', dec='.')

cytof.data = read.csv('cytof/NewRAW.csv', row.names = 1,
                      header = T, stringsAsFactors = F, sep=',', dec='.')

cytof.names = sapply(rownames(cytof.data), function(x) {
  x.ls = strsplit(x, '_')[[1]]
  strsplit(x.ls[length(x.ls)], '\\.')[[1]][1]
}
  )


cytokines.data = data.frame(t(cytokines.data))
cytokines.names = sapply(as.vector(cytokines.data$nmID), function(x) {
  strsplit(x, '_')[[1]][2]
}
)
meta.data = data.frame(t(meta.data))
meta.names = sapply(as.vector(meta.data$SUBJECT.ID), function(x) {
  strsplit(x, '_')[[1]][2]
}
)

head(meta.data)
dim(meta.data)

length(cytokines.names)
length(unique(meta.names))
length(unique(cytof.names))

intersect(cytokines.names, meta.names)
intersect(unique(cytof.names), meta.names)
setdiff(unique(cytof.names), meta.names)
length(intersect(unique(cytof.names), meta.names))
intersect(unique(cytof.names), cytokines.names)

dim(cytof.data[,grep('Freq', colnames(cytof.data))])
length(grep('Median', colnames(cytof.data)))
cytof.data[,grep('Median', colnames(cytof.data))]
