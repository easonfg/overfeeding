kip.cytokines.CFS = read.csv('bmi/results/univariate_csv/1kip.cytokines.CFS_univariate.csv', row.names = 1)
kip.cytokines.SLV = read.csv('bmi/results/univariate_csv/1kip.cytokines.SLV_univariate.csv', row.names = 1)
# kip.cytokines = read.csv('bmi/results/univariate_csv/1kip.cytokines_univariate.csv', row.names = 1)
bmi.all = read.csv('bmi/results/univariate_csv/bmi_all_univariate.csv', row.names = 1)
# sspg = read.csv('bmi/results/univariate_csv/sspg_univariate.csv', row.names = 1)
sspg_bmi = read.csv('bmi/results/univariate_csv/sspg_bmi_univariate.csv', row.names = 1)
sspg_deltawt = read.csv('bmi/results/univariate_csv/sspg_deltaweight_univariate.csv', row.names = 1)
sspg_deltawt_controlled = read.csv('bmi/results/univariate_csv/sspg_deltaweight_controlled_univariate.csv', row.names = 1)
sspg_deltasspg = read.csv('bmi/results/univariate_csv/sspgsspg_deltasspg_univariate.csv', row.names = 1)
montoya = read.csv('bmi/results/univariate_csv/montoya_univariate.csv', row.names = 1)

sspg_bmi$genes = toupper(sub('_FCH', '', sspg_bmi$genes)) %>% sub('\\.', '', .)
sspg_deltawt$genes = toupper(sub('_FCH', '', sspg_deltawt$genes)) %>% sub('\\.', '', .)
sspg_deltawt_controlled$genes = toupper(sub('_FCH', '', sspg_deltawt_controlled$genes)) %>% sub('\\.', '', .)
sspg_deltasspg$genes = toupper(sub('_FCH', '', sspg_deltasspg$genes)) %>% sub('\\.', '', .)
bmi.all$genes = bmi.all$genes %>% toupper() %>% sub('\\.', '_', .)
kip.cytokines.CFS$genes = kip.cytokines.CFS$genes %>% toupper() %>% sub('\\.', '_', .)
kip.cytokines.SLV$genes = kip.cytokines.SLV$genes %>% toupper() %>% sub('\\.', '_', .)

kip.cytokines.CFS$dataset = rep('kip.cytokines.CFS', nrow(kip.cytokines.CFS))
kip.cytokines.SLV$dataset = rep('kip.cytokines.SLV', nrow(kip.cytokines.SLV))
bmi.all$dataset = rep('bmi.all', nrow(bmi.all))
sspg_bmi$dataset = rep('sspg_bmi', nrow(sspg_bmi))
sspg_deltawt$dataset = rep('sspg_deltawt', nrow(sspg_deltawt))
sspg_deltawt_controlled$dataset = rep('sspg_deltawt_controlled', nrow(sspg_deltawt_controlled))
sspg_deltasspg$dataset = rep('sspg_deltasspg', nrow(sspg_deltasspg))
montoya$dataset = rep('montoya', nrow(montoya))

all.cyto.data = bind_rows(kip.cytokines.CFS,kip.cytokines.SLV, bmi.all, sspg_bmi, sspg_deltawt, sspg_deltawt_controlled, sspg_deltasspg)
all.cyto.data
dim(all.cyto.data)
dim(kip.cytokines)
dim(bmi.all)
dim(sspg)

dim(montoya)
head(all.cyto.data)
  
gg.res = ggplot(all.cyto.data, aes(dataset, genes)) + 
  geom_point(aes(size = -log10(fdr), fill = coeff), 
             colour = 'black', shape = 21) + 
  geom_text(aes(label=round(coeff,2)),size = 3, hjust=-0.5, vjust=0.5) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(size = '-log10(FDR)', fill = 'coefficients') +
  ggtitle('Univariate Analysis')
print(gg.res)


### montoya
# gg.res = ggplot(montoya, aes(dataset, genes)) + 
#   geom_point(aes(size = -log10(fdr), fill = coeff), 
#              colour = 'black', shape = 21) + 
#   geom_text(aes(label=round(coeff,2)),size = 3, hjust=-0.5, vjust=0.5) +
#   scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
#   ylab('') + xlab('') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
#   labs(size = '-log10(FDR)', fill = 'coefficients') +
#   ggtitle('Univariate Analysis')
# print(gg.res)
# 
