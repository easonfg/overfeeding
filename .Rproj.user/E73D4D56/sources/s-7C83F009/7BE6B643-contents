kip.cytokines.CFS = read.csv('bmi/results/glmnet_csv/1kip.cytokines.CFS_glmnet.csv', row.names = 1)
kip.cytokines.SLV = read.csv('bmi/results/glmnet_csv/1kip.cytokines.SLV_glmnet.csv', row.names = 1)
# kip.cytokines = read.csv('bmi/results/glmnet_csv/1kip.cytokines_glmnet.csv', row.names = 1)
bmi.all = read.csv('bmi/results/glmnet_csv/bmi_all_glmnet.csv', row.names = 1)
# sspg = read.csv('bmi/results/glmnet_csv/sspg_glmnet.csv', row.names = 1)
sspg_bmi = read.csv('bmi/results/glmnet_csv/sspg_bmi_glmnet.csv', row.names = 1)
sspg_deltawt = read.csv('bmi/results/glmnet_csv/sspg_deltaweight_glmnet.csv', row.names = 1)
sspg_deltawt_controlled = read.csv('bmi/results/glmnet_csv/sspg_deltaweight_controlled_glmnet.csv', row.names = 1)
sspg_deltasspg = read.csv('bmi/results/glmnet_csv/sspgsspg_deltasspg_glmnet.csv', row.names = 1)
montoya = read.csv('bmi/results/glmnet_csv/montoya_glmnet.csv', row.names = 1)

sspg_bmi$name = toupper(sub('_FCH', '', sspg_bmi$name)) %>% sub('\\.', '', .)
sspg_deltawt$name = toupper(sub('_FCH', '', sspg_deltawt$name)) %>% sub('\\.', '', .)
sspg_deltawt_controlled$name = toupper(sub('_FCH', '', sspg_deltawt_controlled$name)) %>% sub('\\.', '', .)
sspg_deltasspg$name = toupper(sub('_FCH', '', sspg_deltasspg$name)) %>% sub('\\.', '', .)
bmi.all$name = bmi.all$name %>% toupper() %>% sub('\\.', '_', .)
kip.cytokines.CFS$name = kip.cytokines.CFS$name %>% toupper() %>% sub('\\.', '_', .)
kip.cytokines.SLV$name = kip.cytokines.SLV$name %>% toupper() %>% sub('\\.', '_', .)

kip.cytokines.CFS = kip.cytokines.CFS[-nrow(kip.cytokines.CFS),]
kip.cytokines.SLV = kip.cytokines.SLV[-nrow(kip.cytokines.SLV),]
bmi.all = bmi.all[-nrow(bmi.all),]
sspg_bmi = sspg_bmi[-nrow(sspg_bmi),]
sspg_deltawt = sspg_deltawt[-nrow(sspg_deltawt),]
sspg_deltawt_controlled = sspg_deltawt_controlled[-nrow(sspg_deltawt_controlled),]
sspg_deltasspg = sspg_deltasspg[-nrow(sspg_deltasspg),]

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

gg.res = ggplot(all.cyto.data, aes(dataset, name)) + 
  geom_point(aes(size = 4, fill = coefficient), 
             colour = 'black', shape = 21) + 
  geom_text(aes(label=round(coefficient,2)),size = 3, hjust=-0.5, vjust=0.5) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
  ylab('') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  labs(size = '-log10(FDR)', fill = 'coefficients') +
  ggtitle('LASSO Analysis')
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
