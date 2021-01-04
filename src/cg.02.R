source("functions.R")
require(CAGEr)
require(BSgenome.Zmays.B73v4)
#require(CAGEfightR)
gen = "BSgenome.Zmays.B73v4"
dirw = file.path(dird, '02_cager')

#{{{ read RNA-seq object
yid = 'cg20a'
#res = rnaseq_cpm(yid)
#th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
#}}}


#{{{ meta plot for all samples
fi = file.path(dirr, '50_final/metaplot.tab')
ti = read_tsv(fi, comment='#')
cnames = colnames(ti)[-1]
ti2 = ti[-ncol(ti)]
colnames(ti2) = cnames
#
ti3 = ti2 %>% gather(sid, cnt) %>%
    filter(!is.na(cnt)) %>%
    separate(sid, c('sid','idx'), sep="_", fill='right') %>%
    separate(sid, c("sid",'srd'), sep='[\\.]') %>%
    mutate(idx = as.integer(idx)) %>%
    replace_na(list(idx=0)) %>%
    mutate(idx = idx+1)

tp = ti3 %>% group_by(sid, srd, idx) %>%
    summarise(cnt=sum(cnt)) %>% ungroup() %>%
    group_by(sid,srd) %>%
    mutate(pct = cnt / abs(sum(cnt))) %>%
    ungroup() %>%
    inner_join(tha, by=c('sid'='SampleID')) %>%
    mutate(lab = str_c(sid, lab, sep=': '))

tps = tibble(idx=c(51,150), lab=c('TSS','TES'))
p = ggplot(tp) +
    geom_line(aes(idx, pct, color=srd)) +
    #geom_vline(xintercept=tps$idx, linetype='dashed', size=.5) +
    scale_x_continuous(breaks=tps$idx, labels=tps$lab, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(expand=expansion(mult=c(.05,.05))) +
    scale_color_aaas(name='strand') +
    facet_wrap(~lab, ncol=5) +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T, strip.style='light',
           xgrid=T,xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T)
fo = file.path(dirw, '18.metaplot.pdf')
ggsave(p, file=fo, width=15, height=30)

fo = file.path(dirw, '18.metaplot.rds')
saveRDS(ti3, fo)
#}}}

#{{{ meta plot for picked samples
fi = file.path(dirw, '18.metaplot.rds')
ti3 = readRDS(fi)

tp = ti3 %>% group_by(sid, srd, idx) %>%
    summarise(cnt=sum(cnt)) %>% ungroup() %>%
    group_by(sid,srd) %>%
    mutate(pct = cnt / abs(sum(cnt))) %>%
    ungroup() %>%
    inner_join(thf, by=c('sid'='SampleID')) %>%
    filter(Treatment %in% c("normal",'cold_control','drought_control')) %>%
    mutate(xpan = str_c(Tissue, Genotype, sep='_')) %>%
    mutate(ypan = sprintf("%s: %s_%s", sid, Treatment, Replicate)) %>%
    mutate(lab = sprintf("%s (%s)", lab, sid))
labs = tp %>% distinct(Tissue,Genotype,Treatment,lab) %>% arrange(Tissue, Genotype, Treatment) %>% pull(lab)
tp = tp %>% mutate(lab = factor(lab, levels=labs))

tps = tibble(idx=c(51,150), lab=c('TSS','TES'))
p = ggplot(tp) +
    geom_line(aes(idx, pct, color=srd)) +
    #geom_vline(xintercept=tps$idx, linetype='dashed', size=.5) +
    scale_x_continuous(breaks=tps$idx, labels=tps$lab, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(expand=expansion(mult=c(.05,.05))) +
    scale_color_aaas(name='strand') +
    facet_wrap(~lab, ncol=4) +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T, strip.style='light',
           xgrid=T,xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T)
fo = file.path(dirw, '19.metaplot.picked.pdf')
ggsave(p, file=fo, width=10, height=20)
#}}}



#{{{ create CAGEexp object [long time]
#fis = sprintf("%s/34_ctss/%s.bw", dirr, thf$SampleID)
fis = sprintf("%s/20_bam/%s.bam", dirr, thf$SampleID)
ce = CAGEexp(genomeName=gen, inputFiles=fis, inputFilesType="bam",
             sampleLabels = thf$lab)
getCTSS(ce, sequencingQualityThreshold = 10, mappingQualityThreshold = 20,
    removeFirstG = F, correctSystematicG = F, useMulticore = F, nrCores = NULL)
#
fo = file.path(dirw, '01.rds')
saveRDS(ce, fo)
#}}}

fi = file.path(dirw, '01.rds')
ce = readRDS(fi)

#{{{ CTSS count stats
tc = CTSStagCountDF(ce) %>% as_tibble()
tc2 = tc %>% mutate(tss=1:n()) %>% gather(sid, cnt, -tss) %>%
    filter(cnt > 0) %>%
    dplyr::count(sid)

tp = tc2 %>% dplyr::rename(lab=sid) %>% inner_join(th, by='lab') %>%
    mutate(tis_gt = str_c(Tissue, Genotype, sep='_')) %>%
    mutate(cond = str_c(Treatment, Replicate, sep="_")) %>%
    select(tis_gt, cond, n) %>%
    spread(tis_gt, n)
fo = file.path(dirw, '03.ctss.cnt.tsv')
write_tsv(tp, fo, na='')
#}}}

#{{{ sample correlation
tis = 'root'
tis = 'shoot'
gt = 'B73'
gt = 'LH143'
#
labs = th %>% filter(Tissue==tis, Genotype==gt) %>% pull(lab)
fo = sprintf("%s/13.corr.%s.%s.pdf", dirw, tis, gt)
pdf(fo, height=8, width=8)
corr.m <- plotCorrelation2(ce, samples = labs,
    tagCountThreshold = 5, applyThresholdBoth = F, method = "pearson")
dev.off()
#}}}

fi = '/scratch.global/zhoux379/nf/work/rnaseq/cg18a/10/2941a9ca63561152b79267012c5425/s03.bam'
x = CAGEr:::import.bam.ctss(fi, filetype='bam', sequencingQualityThreshold=10,mappingQualityThreshold=20,removeFirstG=F,correctSystematicG=F,genome=coerceInBSgenome(gr))
x = CAGEr:::import.bam.ctss(fi, filetype='bam', sequencingQualityThreshold=10,mappingQualityThreshold=20,removeFirstG=F,correctSystematicG=F, genome=NULL)

fi = '~/Downloads/Follicle%20Associated%20Epithelium%2c%20pool2.CNhs13211.10262-104D1.mm10.nobarcode.ctss.bed.gz'
x = CAGEr:::import.bedCTSS(fi)

x2 = as.data.frame(x)
colnames(x2)[1]='chr'
x3 = as(x2, 'CAGEexp')
genomeName(x3) = 'BSgenome.Mmusculus.UCSC.mm10'
CAGEr:::exportCTSStoBedGraph(x3, values = "raw", format = "BigWig")
system(glue("mv score.CTSS.raw.plus.bw {fo1}"))
system(glue("mv score.CTSS.raw.minus.bw {fo2}"))
