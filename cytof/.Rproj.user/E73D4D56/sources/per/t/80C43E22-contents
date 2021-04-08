rm(list=ls())

library(dplyr)
library(tibble)

all.timepoints = data.frame()
for (time in list(list(1,2), list(1,3), list(1,4), list(2,3), list(2,4), list(3,4))){
  time1 = time[[1]]
  time2 = time[[2]]
  assign(paste0('tp', time1, '_', time2), read.csv(paste0('logistic2_more_ctrl_batchCorrectinLM/', time1, '_', time2, 
                                                  '__Median_more_ctrl_BatchCinLM.csv'), row.names = 1,
                            header = T, stringsAsFactors = F, sep=',', dec='.'))
  
  tmp.df = get(paste0('tp', time1, '_', time2))
  tmp.df = tmp.df %>% rownames_to_column(., var = 'cytokines') %>% mutate(timepoint = paste0('tp', time1, '_', time2))
  all.timepoints = bind_rows(all.timepoints, tmp.df)
}

str(all.timepoints)

pdf('logistic2_more_ctrl_batchCorrectinLM/dotplots.res.pdf', onefile = T,
     width=20, height=8)
###plot all timepoints
print(
ggplot(all.timepoints, aes(x=timepoint, y = reorder(cytokines, Log2FC))) + 
  geom_point(shape = 21, aes(fill = Log2FC, size = -log(timepoint2.p.val))) +
  scale_fill_gradient2(low="blue",
                       # scale_color_gradient2(low="blue",
                       high="red", midpoint=mean(0)) + 
  # cowplot::theme_cowplot() + 
  # theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') 
  # theme(axis.ticks = element_blank())
)

## plot ones with multiple time point signficiant
count.cyto = table(all.timepoints$cytokines)
name.multi = count.cyto[count.cyto > 1] %>% names()
all.timepoints[all.timepoints$cytokines %in% name.multi,  ]

print(
ggplot(all.timepoints[all.timepoints$cytokines %in% name.multi,  ], aes(x=timepoint, y = reorder(cytokines, Log2FC))) + 
  geom_point(shape = 21, aes(fill = Log2FC, size = -log(timepoint2.p.val))) +
  scale_fill_gradient2(low="blue",
                       high="red", midpoint=mean(0)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') 
)

dev.off()

# write.csv(all.timepoints[all.timepoints$cytokines %in% name.multi,  ], 'logistic2_more_ctrl_batchCorrectinLM/multiple.signficant.tp.csv')
# write.csv(all.timepoints, 'logistic2_more_ctrl_batchCorrectinLM/all.signficant.tp.csv')
