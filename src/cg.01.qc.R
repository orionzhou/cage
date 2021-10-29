require(CAGEfightR)
require(GenomicFeatures)
source("functions.R")
dirw = glue("{dird}/01_qc")
#setwd(dirw)

#{{{ read RNA-seq object
yid = 'zm.cg20a.2'
r = rnaseq_cpm(yid)
names(r)
#}}}

#{{{ mapping stats
r2 = sum_mapping_stat(r)
tb1 = r2[[1]]; tb2 = r2[[2]]

tp = tb2 %>% inner_join(thf, by=c('sid')) %>%
    select(pnl=grp, tag1=rep, tag2=type, n=cnt) %>% mutate(n=n/1e6)

cols5 = pal_simpsons()(8)[c(3,1,2,4,5)]
p1 = cmp_proportion(tp, alph=.4, acc=.1,
    lab.size=2, nc=3, oneline=F, fills=cols5,
    expand.x=c(.02,.02), expand.y=c(.01,.04), legend.title='',
    legend.pos='top.center.out', legend.dir='h', legend.vjust=-1,
    margin= c(.5,.1,.1,.1), xtext=T, xtick=T) +
    theme(axis.text.x=element_text(size=8, color='red'))
fo = glue("{dirw}/02.map.rate.pdf")
ggsave(p1, filename = fo, width = 8, height = 8)
#}}}

#{{{ meta plot for all samples
ti = r$metaplot

tp = ti %>%
    separate(sid, c('sid','suf'), sep='-') %>%
    group_by(sid,srd) %>%
    mutate(pct = cnt / abs(sum(cnt))) %>%
    ungroup() %>%
    inner_join(tha, by='sid') %>%
    mutate(lab = glue("{sid} {pnl}"))

tps = tibble(idx=c(51,150), lab=c('TSS','TES'))
p = ggplot(tp) +
    geom_line(aes(idx, pct, color=srd)) +
    #geom_vline(xintercept=tps$idx, linetype='dashed', size=.5) +
    scale_x_continuous(breaks=tps$idx, labels=tps$lab, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(expand=expansion(mult=c(.05,.05))) +
    scale_color_aaas(name='strand') +
    facet_wrap(~lab, ncol=5) +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           strip.style='light', strip.compact=T,
           xgrid=T,xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T)
fo = file.path(dirw, '18.metaplot.raw.pdf')
ggsave(p, file=fo, width=10, height=10)
#}}}

#{{{ meta plot for picked samples
ti = r$metaplot

tp = ti %>%
    separate(sid, c('sid','suf'), sep='-') %>% filter(sid %in% thf$sid) %>%
    group_by(sid,srd) %>%
    mutate(pct = cnt / abs(sum(cnt))) %>%
    ungroup() %>%
    inner_join(tha, by='sid') %>%
    mutate(lab = glue("{sid} {pnl}"))

tps = tibble(idx=c(51,150), lab=c('TSS','TES'))
p = ggplot(tp) +
    geom_line(aes(idx, pct, color=srd)) +
    #geom_vline(xintercept=tps$idx, linetype='dashed', size=.5) +
    scale_x_continuous(breaks=tps$idx, labels=tps$lab, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(expand=expansion(mult=c(.05,.05))) +
    scale_color_aaas(name='strand') +
    facet_wrap(~lab, ncol=4) +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           strip.style='light', strip.compact=T,
           xgrid=T,xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T)
fo = file.path(dirw, '18.metaplot.filter.pdf')
ggsave(p, file=fo, width=9, height=9)
#}}}

#{{{ ase
th2 = th %>% rename(SampleID=sid, lab=pnl)
rg = r$ase_gene %>% separate(SampleID, c('SampleID','suf'), sep='-') %>%
    select(-suf) %>% filter(SampleID %in% th2$SampleID)
pa1 = plot_ase(rg, th2, val.col='gt', pal.col='aaas', min_rc=10)
fo = file.path(dirw, '31.afs_gene.pdf')
ggsave(fo, pa1, width=7, height=10)

rs = r$ase_snp %>% separate(SampleID, c('SampleID','suf'), sep='-') %>%
    select(-suf,-gt) %>% filter(SampleID %in% th2$SampleID)
pa2 = plot_ase(rs, th2, val.col='gt', pal.col='aaas', min_rc=10)
fo = file.path(dirw, '32.afs_site.pdf')
ggsave(fo, pa2, width=7, height=10)
#}}}

