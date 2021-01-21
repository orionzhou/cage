source("functions.R")
require(universalmotif)
dirw = glue("{dird}/12_tata")
setwd(dirw)

#{{{ make TATA motifs
fi = 'tata.txt'
mtfs = read_uniprobe(fi)
fo = 'tata.meme'
write_meme(mtfs, fo)
#}}}

#{{{ create promoter db
fi = glue("{dird}/03_qc/01.tss.gtss.rds")
r1 = readRDS(fi)
tss=r1$tss; SE_gtss=r1$SE_gtss; gtss=r1$gtss

fcfg = glue("{dirw}/config.xlsx")
tcfg = read_xlsx(fcfg)
#
chrom_size = read_tsv('00.sizes', col_names=c('chrom','size'))

get_itv <- function(offd, offu, tss, chrom_sizes) {
    #{{{
    tl = tss %>%
        mutate(b = ifelse(srd=='-', end-offd, start-offu)) %>%
        mutate(e = ifelse(srd=='-', end+offu, start+offd)) %>%
        inner_join(chrom_size, by='chrom') %>%
        mutate(b=pmax(0, b), e=pmin(e, size)) %>%
        mutate(score = '.') %>%
        select(chrom, b, e, tidx, score, srd) %>%
        arrange(chrom,b,e)
    tl
    #}}}
}
x = tcfg %>%
    mutate(tl = map2(offd,offu, get_itv, tss=tss, chrom_sizes=chrom_sizes)) %>%
    mutate(fo = glue("{dirw}/01_itv/{opt}.bed")) %>%
    mutate(x1 = map2(tl, fo, write_tsv, col_names=F))

db = "$ref/10.fasta"
system(glue("bedtools getfasta -fi {db} -bed 01.tss.up1k.bed -fo 01.tss.up1k.fas -nameOnly -s"))
#}}}

#{{{ using known TATA motif scan
#system("fimo.py locate tata.meme --motif tata1 01.tss.up1k.fas 04.tata.bed")
system("fimo.py locate tata.meme 01.tss.up1k.fas 04.tata.bed")

fi = glue("{dirw}/10.tata.bed")
ti = read_tsv(fi, col_names=c('hid','beg','end')) %>%
    separate(hid, c('mid','tidx'), sep="%")

itv = 20
tp = ti %>% mutate(pos = (beg+end)/2) %>%
    mutate(pos.bin=cut_width(pos, 20, labels=F, boundary=0)) %>%
    mutate(pos2 = (pos.bin-.5) * itv) %>%
    count(mid, pos2)

tpx = tibble(x=c(0,1000),lab=c('-1kb','TSS'))
p = ggplot(tp, aes(x=pos2,y=n, color=mid)) +
    geom_line(aes(), size=.5, na.rm = F) +
    geom_point(aes(), size=1, na.rm = F) +
    scale_x_continuous(name='position', expand=expansion(mult=c(.05,.05)), breaks=tpx$x, labels=tpx$lab) +
    scale_y_continuous(name="frequency", expand=expansion(mult=c(0,.03)), limits=c(0,NA)) +
    scale_color_npg() +
    facet_wrap(mid ~ ., ncol=3, scale='free') +
    otheme(legend.pos='none', strip.compact=T,
           xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T,
           legend.title=F)
#
fo = file.path(dirw, '10tata.pdf')
ggsave(fo, p, width=6, height=6)
#}}}

# run streme to find enriched motifs

#{{{ streme motifs
opt = 'up100'
opt = 'up1k'
opt = 'up500'
fi = glue("{dirw}/03_streme/{opt}.bed")
ti = read_tsv(fi, col_names=c('hid','beg','end')) %>%
    separate(hid, c('mid','tidx'), sep="%") %>%
    separate(mid, c('midx','conseq'), sep='-', remove=F) %>%
    mutate(midx = as.integer(midx))
tis = ti %>% distinct(midx, mid) %>% arrange(midx) %>%
    mutate(mid = as_factor(mid))# %>% slice(1:12)
ti = ti %>% filter(mid %in% tis$mid) %>%
    mutate(mid = factor(mid, levels=tis$mid))

itv = 10; size = 100
itv = 20; size = 1000
itv = 20; size = 500
tpx = tibble(x=c(0,size),lab=c(glue('-{size}'),'TSS'))
tp = ti %>%
    mutate(pos = (beg+end)/2) %>%
    mutate(pos.bin=cut_width(pos, itv, labels=F, boundary=0)) %>%
    mutate(pos2 = (pos.bin-.5) * itv) %>%
    count(mid, midx, pos2)

nc = 8; wd = 12; ht = 9
nc = 3; wd = 6; ht = 4
nc = 7; wd = 10; ht = 8
p = ggplot(tp, aes(x=pos2,y=n, color=as.factor(midx %% 3))) +
    geom_line(aes(), size=.5, na.rm = F) +
    geom_point(aes(), size=1, na.rm = F) +
    scale_x_continuous(name='position', expand=expansion(mult=c(.1,.1)), breaks=tpx$x, labels=tpx$lab) +
    scale_y_continuous(name="frequency", expand=expansion(mult=c(0,.05)), limits=c(0,NA)) +
    scale_color_npg() +
    facet_wrap(mid ~ ., ncol=nc, scale='free') +
    otheme(legend.pos='none', strip.compact=T,
           xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T,
           legend.title=F)
#
fo = glue('{dirw}/21.mtfs.{opt}.pdf')
ggsave(fo, p, width=wd, height=ht)
#}}}

tss = tss %>% rename(iqr=IQR) %>%
    mutate(shape = map_chr(iqr, iqr2shape)) %>%
    mutate(shape.s = map_chr(iqr, iqr2shape, opt='s')) %>%
    mutate(shape = factor(shape, levels=shapes)) %>%
    mutate(shape.s = factor(shape.s, levels=shapess))
gtss = gtss %>% rename(iqr=IQR) %>%
    mutate(shape = map_chr(iqr, iqr2shape)) %>%
    mutate(shape.s = map_chr(iqr, iqr2shape, opt='s')) %>%
    mutate(shape = factor(shape, levels=shapes)) %>%
    mutate(shape.s = factor(shape.s, levels=shapess))
g2t = gtss %>% select(gidx,tidx) %>% unnest(tidx)

#{{{ obtain TATA
itv = 20
x0 = ti %>% mutate(pos = (beg+end)/2) %>%
    mutate(pos.bin=cut_width(pos, 20, labels=F, boundary=0)) %>%
    mutate(pos2 = (pos.bin-.5) * itv) %>%
    filter(mid == 'tata1', pos2==970) %>%
    select(tidx) %>%
    separate(tidx, c('tidx','srd'), sep='\\(') %>%
    pull(tidx)

itv = 10
x = ti %>%
    mutate(pos = (beg+end)/2) %>%
    mutate(pos.bin=cut_width(pos, itv, labels=F, boundary=0))  %>%
    filter(pos.bin>=7, pos.bin <=9, midx==1) %>%
    select(tidx) %>%
    separate(tidx, c('tidx','srd'), sep='\\(') %>%
    pull(tidx)


tss1 = tss %>% filter(tidx %in% x0)
tss1 %>% count(peakType) %>% mutate(prop = n/sum(n))
tss1 %>% count(shape) %>% mutate(prop = n/sum(n))
#}}}

tg1 = stb %>% select(gid) %>% mutate(tag='stable')
cmps = c('Ar Br', 'As Bs')
tg2 = tg %>% filter(cmp %in% cmps) %>%
    filter(tag=='highly variable') %>%
    count(gid) %>% filter(n==length(cmps)) %>%
    select(gid) %>% mutate(tag='variable')
tg0 = tgl %>% select(gid) %>% mutate(tag = 'all')
tgc = rbind(tg1,tg2,tg0) %>%
    mutate(tag = as_factor(tag))
tgcs = tgc %>% count(tag) %>% mutate(pnl = glue("{tag} ({n})")) %>%
    mutate(pnl = as_factor(pnl))
tgc = tgc %>% inner_join(tgcs, by='tag')

#{{{ characterize TATA-box genes
tata = tss1 %>% inner_join(g2t, by='tidx') %>%
    distinct(gid) %>%
    mutate(tag = 'tata') %>%
    inner_join(tgl, by='gid') %>% print(n=30)
tge = stb %>% select(gid,tag) %>%
    bind_rows(tata %>% select(gid,tag))

#{{{ GO
tgo = read_go(src='Interproscan5')
tgo = read_go(src='arabidopsis')
tgo = read_go(src='uniprot.plants')
tgrp = tgo %>% select(gotype, grp=goid, gid, note=goname) %>%
    group_by(gotype) %>% nest() %>% rename(tgrp = data)

x0 = tge %>% rename(group = tag) %>%
    group_by(group) %>% summarise(gids = list(gid)) %>% ungroup() %>%
    crossing(tgrp) %>%
    mutate(x = map2(gids, tgrp, hyper_enrich))
x1 = x0 %>% select(group,gotype,x) %>% unnest(x) %>%
    arrange(group, pval.raw) %>%
    select(-grp) %>%
    group_by(group) %>% slice(1:15) %>% ungroup() %>% select(-group,-gotype) %>%
    #filter(group=='stable') %>%
    print(n=50)
x1$note[[1]]
#}}}

#{{{ shape
fi = glue("{dird}/03_qc/02.tss.gtss.cond.rds")
r2 = readRDS(fi)
iqr2shape <- function(iqr, opt='l') {
    #{{{
    if (opt == 's') {
        vals = shapess
    } else {
        vals = shapes
    }
    ifelse(iqr==0, vals[1], ifelse(iqr <= 10, vals[2], vals[3]))
    #}}}
}
gtss = r2$gtss %>% rename(iqr=IQR) %>%
    mutate(shape = map_chr(iqr, iqr2shape)) %>%
    mutate(shape.s = map_chr(iqr, iqr2shape, opt='s')) %>%
    mutate(shape = factor(shape, levels=shapes)) %>%
    mutate(shape.s = factor(shape.s, levels=shapess))

gshape = gtss %>%
    inner_join(gcfg$gene %>% select(gid,ttype), by='gid') %>%
    filter(ttype=='mRNA') %>% select(gid,tag2=shape.s)

tp0 = gshape
tpc = tp0 %>% mutate(tag1='all')
tp = tge %>% select(gid,tag1=tag) %>%
    inner_join(tp0, by='gid') %>%
    bind_rows(tpc) %>% mutate(tag1 = as_factor(tag1)) %>%
    count(tag1, tag2)
p1 = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = pal_npg()(5)) +
    o_margin(.1,.5,.1,.5)
#}}}
fo = glue("{dirw}/12.shape.pdf")
ggsave(p1, file=fo, width=4, height=4)

#{{{ tis-specificity
tis_spec = read_tis_spec() %>%
    inner_join(gcfg$gene %>% select(gid,ttype), by='gid') %>%
    filter(ttype=='mRNA') %>% select(gid,tag2=tag)

tp0 = tis_spec
tpc = tp0 %>% mutate(tag1='all')
tp = tge %>% select(gid,tag1=tag) %>%
    inner_join(tp0, by='gid') %>%
    bind_rows(tpc) %>% mutate(tag1 = as_factor(tag1)) %>%
    count(tag1, tag2)
p = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = brewer.pal(5,'Set2')) +
    o_margin(.3,.3,.3,.3)
#}}}
fo = glue("{dirw}/12.tis.pdf")
ggsave(p, file=fo, width=4, height=4)

#{{{ syn
syn = read_syn(gcfg) %>% arrange(gid, -ftype) %>% group_by(gid) %>%
    slice(1) %>% ungroup() %>%
    inner_join(gcfg$gene %>% select(gid,ttype), by='gid') %>%
    filter(ttype=='mRNA') %>% select(gid,tag2=ftype)

tp0 = syn
tpc = tp0 %>% mutate(tag1='all')
tp = tge %>% select(gid,tag1=tag) %>%
    inner_join(tp0, by='gid') %>%
    bind_rows(tpc) %>% mutate(tag1 = as_factor(tag1)) %>%
    count(tag1, tag2)
p = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = brewer.pal(5,'Set3')) +
    o_margin(.1,.1,.1,.1)
fo = glue("{dirw}/12.syn.pdf")
ggsave(p, file=fo, width=4, height=4)

x = read_syn(gcfg, opt=2) %>%
    left_join(tata %>% select(maize1=gid,tag1=tag), by='maize1') %>%
    left_join(tata %>% select(maize2=gid,tag2=tag), by='maize2') %>%
    replace_na(list(tag1='non-tata', tag2='non-tata'))
x %>% count(ftype, tag1,tag2)
#}}}

#{{{ utr5 length
tg0 = gcfg$gene.loc %>% group_by(gid, ttype) %>%
    mutate(size = end - start + 1) %>%
    summarise(size.utr5 = sum(size[etype=='five_prime_UTR'])) %>%
    ungroup() %>% filter(ttype=='mRNA') %>% select(-ttype) %>%
    mutate(tag2 = ifelse(size.utr5 == 0, "UTR5 = 0", "UTR5 > 0"))
#tg0 = tg0 %>% filter(size.utr5 > 0)

tp = tg0 %>%
    left_join(tge %>% select(gid,tag1=tag), by='gid') %>%
    replace_na(list(tag1='non-tata'))
tp %>% count(tag1,tag2)
tp %>% filter(size.utr5 > 0) %>%
    group_by(tag1) %>% skim(size.utr5)

p1 = cmp_proportion1(tp,xangle=0, oneline=T,legend.title='', barwidth=.8)
fo = file.path(dirw, "33.utr5.1.pdf")
ggsave(p1, file=fo, width=3, height=4)

tpc = tg0 %>% mutate(tag1='all genes') %>% select(gid,tag1,size.utr5)
tp = tg0 %>% mutate(tag1='stable') %>%
    inner_join(tgl, by='gid') %>%
    bind_rows(tpc) %>%
    mutate(tag1 = as_factor(tag1))

#{{{ boxplot
tps = tp %>% filter(size.utr5>0) %>% rename(ctag=tag1, score=size.utr5) %>% group_by(ctag) %>%
  summarise(ng = n(), y = median(score),
            y25 = quantile(score,.25), y75=quantile(score, .75),
            y5 = quantile(score,.05), y95=quantile(score,.95)) %>% ungroup()
pv = ggplot(tps) +
    #geom_violin(aes(x=ctag, y=score), trim=T, alpha=.8) +
    #geom_boxplot(aes(x=ctag, y=score), outlier.shape=NA, width=.3) +
    geom_boxplot(aes(x=ctag, ymin=y5,lower=y25, middle=y, upper=y75, ymax=y95, fill=ctag),
                 stat='identity', width=.3, position=position_dodge(.5)) +
    scale_x_discrete(expand=expansion(mult=c(.2,.2))) +
    scale_y_continuous(name='UTR5 length (bp)', expand=expansion(mult=c(.05,.05))) +
    scale_fill_manual(values=pal_jco()(8)) +
    scale_color_aaas(name='') +
    otheme(legend.pos='none', strip.style='white',
           xtext=T, xtick=T, ytext=T, ytick=T, ytitle=T, ygrid=T)
fo = file.path(dirw, '33.utr5.2.pdf')
ggsave(pv, filename=fo, width=3, height=4)
#}}}
#}}}

#{{{ nearest TE distance
ft = file.path(dirw, '36.TE.bed')
tt = read_tsv(ft, col_names=F)
tt2 = tt %>% filter(str_detect(X6, "^(RL|DT)"))
fo = glue("{dirw}/36.TE.1.bed")
write_tsv(tt2, fo, col_names=F)
fo = glue("{dirw}/36.TE.bed")
system(glue("closestBed -D a -t first -a 36.gene.bed -b {fo} > 36.bed"))

fd = file.path(dirw, '36.bed')
td = read_tsv(fd, col_names=F)
tg0 = td %>%
    select(tid=X4, dst = X16) %>%
    separate(tid, c('gid','iso'), sep='_') %>% select(-iso) %>%
    mutate(score = pmin(abs(dst), 100000)) %>%
    select(gid, score) %>%
    mutate(tag2 = ifelse(score==0, 'ovlp TE', 'not ovlp TE'))

tpc = tg0 %>% mutate(tag1='ctrl-bg') %>% select(gid,tag1,tag2)
tp = tg0 %>% mutate(tag1='stable') %>%
    inner_join(tgl, by='gid') %>%
    replace_na(list(tag1 = 'ctrl-bg')) %>%
    bind_rows(tpc) %>%
    count(tag1, tag2) %>%
    mutate(tag1 = as_factor(tag1)) %>%
    mutate(tag2 = as_factor(tag2))
p1 = cmp_proportion1(tp,xangle=0, oneline=T,legend.title='', barwidth=.8)
fo = file.path(dirw, "35.TE.1.pdf")
ggsave(p1, file=fo, width=3, height=4)

tpc = tg0 %>% mutate(tag1='ctrl-bg') %>% select(gid,tag1,score)
tp = tg0 %>% mutate(tag1='stable') %>%
    inner_join(tgl, by='gid') %>%
    bind_rows(tpc) %>%
    mutate(tag1 = as_factor(tag1))
#{{{ boxplot
tps = tp %>% filter(score>0) %>% rename(ctag=tag1) %>% group_by(ctag) %>%
  summarise(ng = n(), y = median(score),
            y25 = quantile(score,.25), y75=quantile(score, .75),
            y5 = quantile(score,.05), y95=quantile(score,.95)) %>% ungroup()
pv = ggplot(tps) +
    #geom_violin(aes(x=ctag, y=score), trim=T, alpha=.8) +
    #geom_boxplot(aes(x=ctag, y=score), outlier.shape=NA, width=.3) +
    geom_boxplot(aes(x=ctag, ymin=y5,lower=y25, middle=y, upper=y75, ymax=y95, fill=ctag),
                 stat='identity', width=.3, position=position_dodge(.5)) +
    scale_x_discrete(expand=expansion(mult=c(.2,.2))) +
    scale_y_continuous(name='distance to nearest upstream TE (bp)', expand=expansion(mult=c(.05,.05))) +
    scale_fill_manual(values=pal_jco()(8)) +
    scale_color_aaas(name='') +
    otheme(legend.pos='none', strip.style='white',
           xtext=T, xtick=T, ytext=T, ytick=T, ytitle=T, ygrid=T)
fo = file.path(dirw, '35.TE.2.pdf')
ggsave(pv, filename=fo, width=3, height=4)
#}}}


#}}}
#}}}

#{{{ check prop. genes w. TATA boxes
tp = tgc %>% select(gid,tag1 = pnl) %>%
    left_join(tata, by='gid') %>%
    rename(tag2 = tag) %>%
    replace_na(list(tag2 = 'non-tata')) %>%
    mutate(tag2 = factor(tag2, levels=c('tata','non-tata'))) %>%
    count(tag1, tag2)

p1 = cmp_proportion1(tp, ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = pal_npg()(5)) +
    o_margin(.1,.5,.1,.5)
fo = glue("{dirw}/13.stvar.tata.pdf")
ggsave(p1, file=fo, width=4, height=4)
#}}}

#{{{ TATA shape
tp = gtss %>% select(gid,tag1 = shape.s) %>%
    left_join(tata, by='gid') %>%
    rename(tag2 = tag) %>%
    replace_na(list(tag2 = 'non-tata')) %>%
    mutate(tag2 = factor(tag2, levels=c('tata','non-tata'))) %>%
    count(tag1, tag2)

p1 = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = pal_npg()(5)) +
    o_margin(.1,.5,.1,.5)
fo = glue("{dirw}/13.shape.tata.pdf")
ggsave(p1, file=fo, width=4, height=4)
#}}}

#{{{ st var chharacterz

#{{{ shape
gshape = gtss %>%
    inner_join(gcfg$gene %>% select(gid,ttype), by='gid') %>%
    filter(ttype=='mRNA') %>% select(gid,tag2=shape.s)

tp0 = gshape
tp = tgc %>% select(gid,tag1=pnl) %>%
    inner_join(tp0, by='gid') %>%
    count(tag1, tag2)
p1 = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = pal_npg()(5)) +
    o_margin(.1,.5,.1,.5)
#}}}
fo = glue("{dirw}/14.st.var.shape.pdf")
ggsave(p1, file=fo, width=4, height=4)

#{{{ tis-specificity
tis_spec = read_tis_spec() %>%
    inner_join(gcfg$gene %>% select(gid,ttype), by='gid') %>%
    filter(ttype=='mRNA') %>% select(gid,tag2=tag)

tp0 = tis_spec
tp = tgc %>% select(gid,tag1=pnl) %>%
    inner_join(tp0, by='gid') %>%
    count(tag1, tag2)
p = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = brewer.pal(5,'Set2')) +
    o_margin(.3,.3,.3,.3)
#}}}
fo = glue("{dirw}/14.st.var.tis.pdf")
ggsave(p, file=fo, width=4, height=4)
#}}}

