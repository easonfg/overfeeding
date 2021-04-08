library(dplyr)

groups = read.csv('cytof/cytofs_w_phenotypes.csv', row.names = 1,
                  header = F, stringsAsFactors = F, sep=',', dec='.')
groups = data.frame(t(groups))
groups$V1

# sspg = read.csv('Re__Obesity_data_set_2/obesityData_luminexWithSSPG.csv', row.names = 1)
sspg = read.csv('Re__Obesity_data_set_2/overfeeding_Cytokine_Chemokine_luminexWithSSPG.csv', skip = 1, row.names = 1)
sspg = t(sspg)
sspg_fch = sspg %>% data.frame()%>% select(c(1,2,4,5, grep('_FCH', colnames(sspg))))
str(sspg_fch)
sspg_base = sspg %>% data.frame()%>% select(c(grep('_B', colnames(sspg)))) 
str(sspg_base)
sspg_peak = sspg %>% data.frame()%>% select(c(grep('_S', colnames(sspg))))


res = sapply(1:ncol(sspg_base), function(x){
  test = t.test(sspg_base[,x], sspg_peak[,x], paired = TRUE, alternative = "two.sided")
  logFC = log(mean(sspg_peak[,x])/mean(sspg_base[,x]), 2)
  c(LogFC = logFC, p.val = test$p.value)
}
)
# res
# dim(res)
# str(data.frame(res))
res = data.frame(t(res))
rownames(res) = sub('_B', '', colnames(sspg_base))
colnames(sspg_base)
res[1:10,1:2]
# str(res)

res$padj = p.adjust(res$p.val, method = 'BH')
res[which(res$padj < 1),]
head(res)
rownames(res)

library(EnhancedVolcano)
jpeg(paste('luminex/', 'foldchange.jpeg', sep = ''),
     units="in", width=10, height=10, res=500)
volcano.plot = EnhancedVolcano(res,
                               lab = rownames(res),
                               # pCutoff = 0.05,
                               pCutoff = 0.1,
                               FCcutoff = 0,
                               drawConnectors = TRUE,
                               x = 'LogFC',
                               y = 'p.val',
                               legendPosition = 'bottom')
print(volcano.plot)
dev.off()
