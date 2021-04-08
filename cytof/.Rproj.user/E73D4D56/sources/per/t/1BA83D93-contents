rm(list=ls())

library(dplyr)
cytof.data = read.csv('NewRAW.csv', row.names = 1,
                      header = T, stringsAsFactors = F, sep=',', dec='.')

groups = read.csv('cytofs_w_phenotypes.csv', row.names = 1,
                  header = F, stringsAsFactors = F, sep=',', dec='.')
groups = data.frame(t(groups))
# groups[1:10,1:3]
groups = groups[,c(1,3)]
groups

controls = read.csv('controls.csv', row.names = 1,
                    header = T, stringsAsFactors = F, sep=',', dec='.')

## change control columns to numeric
controls[, sapply(controls, class) == 'character'] = as.numeric(as.matrix(controls[, sapply(controls, class) == 'character']))

str(controls)
dim(controls)
ctrl_correct = rbind(controls[2,]/controls[1,],
                     controls[3,]/controls[1,],
                     controls[4,]/controls[1,],
                     controls[5,]/controls[1,])

## change cytof columns to numeric
cytof.data[, sapply(cytof.data, class) == 'character'] = as.numeric(as.matrix(cytof.data[, sapply(cytof.data, class) == 'character']))

## make annotation matrix
annotations = data.frame(batches = as.factor(sapply(rownames(cytof.data), function(x) strsplit(x, '_')[[1]][2])),
                         timepoints = as.factor(cytof.data$Time))
cytof.data = cytof.data %>% select(-c('Time'))
dim(cytof.data)
annotations$ids = sapply(rownames(annotations), function(x) {
  ls.x = strsplit(x, '_')[[1]]
  # print(ls.x[length(ls.x)])
  ls.x = strsplit(ls.x[length(ls.x)], '\\.')[[1]][1]
  sub('^0', '', ls.x)
}
)
head(groups)
head(annotations)
# length(intersect(unique(annotations$ids), unique(groups$V1)))
annotations$groups = groups$insulin.status[match(annotations$ids, groups$V1)]

## delete columns with one values
non.unique.col = which(apply(cytof.data, 2, function(x) length(unique(x))) <= 1)
cytof.data = cytof.data[,-non.unique.col]
ctrl_correct = ctrl_correct[,-non.unique.col]

# delete na columns
col_na = which(colSums(is.na(cytof.data)) > 0)
cytof.data = cytof.data[,-col_na]
dim(cytof.data)
ctrl_correct = ctrl_correct[,-col_na]
ctrl_correct %>% dim()
cytof.data %>% dim()

### original pca
pca.cytof.res = prcomp(cytof.data,  scale = T)
library("factoextra")
fviz_eig(pca.cytof.res)
pca.plot = fviz_pca_ind(pca.cytof.res,
                        geom="point",
                        col.ind = annotations$batches, # Color by the quality of representation
                        repel = TRUE     # Avoid text overlapping
)

print(pca.plot)

library(sva)
cytof.combat = t(ComBat(t(cytof.data), batch = annotations$batch))

pca.cytof.combat.res = prcomp(cytof.combat,  scale = T)
library("factoextra")
fviz_eig(pca.cytof.combat.res)
pca.plot = fviz_pca_ind(pca.cytof.combat.res,
                        geom="point",
                        col.ind = annotations$batches, # Color by the quality of representation
                        repel = TRUE     # Avoid text overlapping
)
print(pca.plot)


library(stringr)
single.cyto = cytof.combat[,grep('Akt', colnames(cytof.combat))] %>% colnames() %>%
  # single.cyto = cytof.data[,grep('Akt', colnames(cytof.combat))] %>% colnames() %>%
  grep('Granu',., value = T) %>%
  # single.cyto = cytof.combat[,grep('pSTAT3', colnames(cytof.combat))] %>% colnames() %>% 
  #               grep('B.cells',., value = T) %>% 
  # single.cyto = cytof.combat[,grep('DC', colnames(cytof.combat))] %>% colnames() %>% 
  #               grep('DC',., value = T) %>% 
  cytof.combat[,.,drop = F] %>%
  # cytof.data[,.,drop = F] %>%
  data.frame()
# newname = (strsplit(colnames(single.cyto)[1], '\\.')[[1]] %>% last())
# single.cyto %>% dplyr::rename(newname = colnames(single.cyto)[1])
single.cyto$timepoint = strsplit(rownames(single.cyto), '_') %>% sapply(., function(x){x[length(x)-1]})
single.cyto$timepoint = single.cyto$timepoint %>% sub('v', '', .) %>% 
  as.factor()
# as.numeric()
single.cyto

single.cyto$ids = strsplit(rownames(single.cyto), '_') %>% sapply(., function(x){
  strsplit(x[length(x)],'\\.')[[1]][1]
}) %>% as.integer()


single.cyto[,1] = log(single.cyto[,1])
single.cyto

single.cyto %>%
  ggplot(aes_string(x = 'timepoint', y = colnames(single.cyto)[1])) + 
  # geom_boxplot()
  geom_boxplot() +
  geom_line(aes(group=ids), position = position_dodge(0.2)) +
  geom_point(aes(fill=timepoint,group=ids), position = position_dodge(0.2)) +
  theme(legend.position = "none") %>% 
  print()
# geom_jitter(shape=16, position=position_jitter(0.2))

# (single.cyto[grep('1', single.cyto$timepoint),1])
# (single.cyto[grep('2', single.cyto$timepoint),1])
# 
# mean(single.cyto[grep('1', single.cyto$timepoint),1])
# sd(single.cyto[grep('1', single.cyto$timepoint),1])
# mean(single.cyto[grep('2', single.cyto$timepoint),1])
# sd(single.cyto[grep('2', single.cyto$timepoint),1])
# mean(single.cyto[grep('3', single.cyto$timepoint),1])
# mean(single.cyto[grep('4', single.cyto$timepoint),1])
# 
# t.test(single.cyto[grep('1', single.cyto$timepoint),1], single.cyto[grep('4', single.cyto$timepoint),1]) 

##### plotting signifciant ones
tp4 = read.csv('logistic_volcano/4p.val__Median.csv', row.names = 1,
               header = T, stringsAsFactors = F, sep=',', dec='.')
tp3 = read.csv('logistic_volcano/3p.val__Median.csv', row.names = 1,
               header = T, stringsAsFactors = F, sep=',', dec='.')
tp2 = read.csv('logistic_volcano/2p.val__Median.csv', row.names = 1,
               header = T, stringsAsFactors = F, sep=',', dec='.')

tp4 = tp4[tp4$p.val < 0.05,]
tp3 = tp3[tp3$p.val < 0.05,]
tp2 = tp2[tp2$p.val < 0.05,]

tp43 = intersect(rownames(tp4), rownames(tp3))
tp32 = intersect(rownames(tp2), rownames(tp3))
tp42 = intersect(rownames(tp4), rownames(tp2))
tp43
tp32
tp42
intersect.sig = intersect(intersect(rownames(tp2), rownames(tp3)), rownames(tp4))
all.sig = c(rownames(tp2), rownames(tp3), rownames(tp4))
write(unique(all.sig), 'single_plot/sig_cytokines.csv')
unique(all.sig) %>% length()

### plot pipe
int.list = unique(c(tp43, tp32, tp42))
int.list = all.sig[grep('DC', all.sig)]
int.list
# pdf(paste0("single_plot/intersecting_cyto.pdf"), onefile = TRUE)
pdf(paste0("single_plot/DC.pdf"), onefile = TRUE)
for (cyto in int.list){
  single.cyto = cytof.combat[,cyto, drop = F] %>%
    data.frame()
  # newname = (strsplit(colnames(single.cyto)[1], '\\.')[[1]] %>% last())
  # single.cyto %>% dplyr::rename(newname = colnames(single.cyto)[1])
  single.cyto$timepoint = strsplit(rownames(single.cyto), '_') %>% sapply(., function(x){x[length(x)-1]})
  single.cyto$timepoint = single.cyto$timepoint %>% sub('v', '', .) %>% 
    as.factor()
  # as.numeric()
  head(single.cyto)
  
  single.cyto$ids = strsplit(rownames(single.cyto), '_') %>% sapply(., function(x){
    strsplit(x[length(x)],'\\.')[[1]][1]
  }) %>% as.integer()
  
  single.cyto[,1] = log(single.cyto[,1])
  head(single.cyto)
  
  # single.cyto %>%
  #   ggplot(aes_string(x = 'timepoint', y = colnames(single.cyto)[1])) + 
  #   # geom_boxplot()
  #   geom_boxplot() +
  #   geom_line(aes(group=ids), position = position_dodge(0.2)) +
  #   geom_point(aes(fill=timepoint,group=ids), position = position_dodge(0.2)) +
  #   theme(legend.position = "none") %>% 
  #   print()
  # 
  
  first.tp = single.cyto[single.cyto$timepoint == 1,]
  second.tp = single.cyto[single.cyto$timepoint == 2,]
  third.tp = single.cyto[single.cyto$timepoint == 3,]
  fourth.tp = single.cyto[single.cyto$timepoint == 4,]
  rest.tp = single.cyto[!(single.cyto$timepoint == 1),]
  rest.tp$facet_group = as.numeric(rest.tp$timepoint)
  faceted.tp = bind_rows(rest.tp, first.tp %>% mutate(facet_group = rep(2,nrow(first.tp))))
  faceted.tp = bind_rows(faceted.tp, first.tp %>% mutate(facet_group = rep(3,nrow(first.tp))))
  faceted.tp = bind_rows(faceted.tp, first.tp %>% mutate(facet_group = rep(4,nrow(first.tp))))
  faceted.tp = bind_rows(faceted.tp, second.tp %>% mutate(facet_group = rep(5,nrow(second.tp))))
  faceted.tp = bind_rows(faceted.tp, third.tp %>% mutate(facet_group = rep(5,nrow(third.tp))))
  faceted.tp = bind_rows(faceted.tp, second.tp %>% mutate(facet_group = rep(6,nrow(second.tp))))
  faceted.tp = bind_rows(faceted.tp, fourth.tp %>% mutate(facet_group = rep(6,nrow(fourth.tp))))
  faceted.tp = bind_rows(faceted.tp, third.tp %>% mutate(facet_group = rep(7,nrow(third.tp))))
  faceted.tp = bind_rows(faceted.tp, fourth.tp %>% mutate(facet_group = rep(7,nrow(fourth.tp))))
  faceted.tp
  
  print(
  ggplot(faceted.tp,aes_string(x= 'timepoint', y = colnames(faceted.tp)[1], group = 'timepoint', fill='timepoint')) +
    geom_boxplot() +
    geom_line(aes(group=ids), position = position_dodge(0.2)) +
    geom_point(aes(fill=timepoint,group=ids), shape = 21, position = position_dodge(0.2)) +
    theme(legend.position = "none") +
    facet_grid(facet_group ~ ., scales="free") 
  )
    
}
dev.off()

single.cyto %>%
  ggplot(aes_string(x = 'timepoint', y = colnames(single.cyto)[1])) +
  # geom_boxplot()
  geom_boxplot() +
  geom_line(aes(group=ids), position = position_dodge(0.2)) +
  geom_point(aes(fill=timepoint,group=ids), position = position_dodge(0.2)) +
  theme(legend.position = "none") %>%
  print()
