sspg = c('sCD40L_FCH', 'RESISTIN_FCH', 'TRAIL_FCH', 'sICAM.1_FCH', 'IL.6_FCH', 'IFN.g_FCH', 'IL.12P40_FCH', 'IL.13_FCH', 'TGF.b_FCH', 'Age', 'LIF_FCH', 'MCP.3_FCH')
montoya = c('RPGRIP1','HLA.DMB','KIR2DL3','TNFRSF21','MSRB2','PADI4','SPRYD5','GENDER','POU2F2','C5orf32','C16orf35','FLNA','ALAS2','RBM38','MIR130A','PI3','XRN2','SLPI','KRT1','ADM')
bmi_all = c('Age','mod_085','pre.gmt','cd20.IL2.STAT3','IL.12P40','CD8.EM','mono.IL10.STAT5','cd20.IL7.STAT3','Fch.norm.Akt','mod_106','LEPTIN','sFAS.ligand','MIG','IL.17F','NK.cells','CD4.CM','CD8.CM','cd4.IL7.STAT1','TNF.b','IFN.g','TGF.b','mono.IL21.STAT1','mod_028','mod_072','cd8.IFNa.STAT3','CD4.CD28.','cd4.IL21.STAT1')
full_1kip = c('TRANSITIONALBCELLS','RANTES','VCAM1','NAiVE_CD4_TCELLS','IGDnegCD27negBCELLS','GMCSF','IL7','IL1B','PLASMABLASTS','BCELLS_S1_IL21','MIG','EOTAXIN','ENA78','IP10','CD8_TCELLS','CD4_S1_IL21','EFF_CD4_TCELLS','VEGF','MONO_S5_IL6','MONO_S3_IL6','AGE','CM_CD4_TCELLS','GCSF','GENDER','LEPTIN')

sspg = gsub('_FCH', '', sspg)
sspg = toupper(gsub('\\.', '', sspg))

montoya = toupper(gsub('\\.', '', montoya))

bmi_all = toupper(gsub('\\.', '', bmi_all))
full_1kip = toupper(gsub('\\.', '', full_1kip))

intersect(sspg, montoya)
intersect(sspg, bmi_all)
intersect(sspg, full_1kip)

intersect(montoya, bmi_all)
intersect(montoya, full_1kip)

intersect(bmi_all, full_1kip)
