source("functions.R")
dirw = glue("{dird}/10_shift")
#setwd(dirw)

#{{{ add stype, IBD, prop_tpm columns to rtss tibble
fi = glue('{dirw}/03.rep.merged.rds')
r3 = readRDS(fi)
rtss=r3$rtss; tss=r3$tss; gtss=r3$gtss

# stype & IBD
stypes=c('identical','1_snp','2_snp','other')
fi = glue('{dird}/06_snp/11.tss.vnt.rds')
tv = readRDS(fi)$tv %>% select(tidx=tid,stype,ibd) %>%
    mutate(stype=factor(stype,levels=stypes)) %>%
    mutate(ibd = factor(ibd))

# proportion TPM in multi-TSS gene 
dtss = rtss %>% filter(!peakType %in% c("intergenic", "antisense")) %>%
    select(tidx, tpm, gid) %>% group_by(gid) %>%
    mutate(prop_tpm = tpm/sum(tpm)) %>% ungroup() %>%
    select(tidx, prop_tpm)

rtss2 = rtss %>% inner_join(tv, by='tidx') %>%
    left_join(dtss, by='tidx') %>%
    mutate(dominant=ifelse(is.na(prop_tpm), NA, prop_tpm>=.15)) %>%
    select(tidx, tpm, IQR, entropy, support, tid, peakType, gid,
        npos, chrom, start, end, srd, peakPos=domi, 
        vtype=stype, ibd, prop_tpm, dominant,
        tpm.grp, cidxs, cnts)
fo = glue("{dirw}/04.rtss.rds")
saveRDS(rtss2, fo)
#}}}

#{{{ read shifting scores - tg
fi = glue("{dirw}/04.rtss.rds")
rtss = readRDS(fi)
fi = glue('{dird}/03_qc/05.shift.rds')
r5 = readRDS(fi)
ss = r5$rtss; 
gene = gcfg$gene %>% mutate(loc=glue("{chrom}:{start}-{end}")) %>% select(gid,loc)

itvs = c(seq(0,1,by=.1))
tags = c("conserved",'variable','highly variable')
tg = tss %>% select(i,gid, cmp) %>% unnest(cmp) %>%
    mutate(cmp = glue("{grp1} {grp2}")) %>%
    mutate(cmp = as.character(cmp)) %>%
    arrange(cmp, desc(shift)) %>%
    group_by(cmp) %>%
    mutate(rank.shift=1:n()) %>% ungroup() %>%
    #arrange(comp, pval.ks.raw) %>%
    #group_by(comp) %>%
    #mutate(rank.ks = 1:n()) %>%
    #mutate(padj = p.adjust(pval.ks.raw)) %>% ungroup() %>%
    #mutate(logp = -log10(pval.ks.raw)) %>%
    #mutate(logp = ifelse(is.infinite(logp), 18, logp)) %>%
    #mutate(shift = pmax(-1, shift)) %>%
    mutate(bin.shift = cut(shift,breaks=itvs, right=T, include.lowest=T, ordered_result=T)) %>%
    #mutate(sig.ks = padj < 0.05) %>%
    mutate(tag = ifelse(shift>.5, tags[3], ifelse(shift<=.2, tags[1], tags[2]))) %>%
    mutate(tag = factor(tag, levels=tags))
tg %>% count(grp1,grp2,tag) %>%
    group_by(grp1,grp2) %>% mutate(p=n/sum(n)) %>% ungroup() %>%
    mutate(txt = glue("{n} ({percent(p,accuracy=.1)})")) %>%
    select(grp1, grp2, tag, txt) %>% spread(tag, txt)
x = tss %>% select(i,tpm.grp) %>% unnest(tpm.grp) %>% select(i,grp,tpm)
tg2 = tg %>% inner_join(x, by=c('i'='i','grp1'='grp')) %>% rename(tpm1=tpm) %>%
    inner_join(x, by=c('i'='i','grp2'='grp')) %>% rename(tpm2=tpm)

cmps = c('root.B root.A', 'shoot.B shoot.A', 'root.B shoot.B', 'root.A shoot.A')
min_tpm = 3
tp = tg2 %>% filter(cmp %in% cmps) %>%
    filter(i %in% tids_domi) %>%
    filter(tpm1>min_tpm, tpm2>=min_tpm) %>%
    inner_join(tv, by='i') %>%
    count(cmp, stype, tag) %>% select(pnl=cmp, tag1=stype, tag2=tag, n)
    #filter(stype=='identical') %>% count(cmp, ibd, tag) %>% select(pnl=cmp, tag1=ibd, tag2=tag, n)
tp %>% group_by(pnl, tag1) %>% mutate(p=n/sum(n)) %>% ungroup() %>%
    mutate(txt = glue("{n} ({percent(p,accuracy=.1)})")) %>%
    select(pnl,tag1,tag2,txt) %>% spread(tag2, txt)
tp %>% group_by(pnl) %>% summarise(n=sum(n))

cols5 = pal_simpsons()(8)#[c(3,1,2,4,5)]
p1 = cmp_proportion(tp, alph=.8, acc=1,
    lab.size=2, nc=4, oneline=T, fills=cols5, strip.compact=F, panel.spacing=.2,
    expand.x=c(.02,.02), expand.y=c(.01,.04), legend.title='',
    legend.pos='top.center.out', legend.dir='h', legend.vjust=-1,
    margin= c(.5,.1,.1,.1), xtext=T, xtick=T) +
    theme(axis.text.x=element_text(size=8, color='black'))
fo = glue("{dirw}/05.3.ibd.pdf")
fo = glue("{dirw}/05.1.tpm3.pdf")
fo = glue("{dirw}/05.0.all.pdf")
fo = glue("{dirw}/05.2.domi.pdf")
ggsave(p1, filename=fo, width=9, height=4)
#}}}

#{{{ shifting score distribution
#{{{ scatter plot - f2a & b
require(hexbin)

x='Br Bs'; y='Ar As'; ctag='genotype-wise'
x='As Bs'; y='Ar Br'; ctag='tissue-wise'
plot_point_dens <- function(tg, x, y, ctag, cmps) {
    #{{{
    tp = tg %>% filter(cmp %in% c(x,y)) %>%
        select(i,gid,cmp,shift) %>%
        spread(cmp, shift) %>%
        select(i, x=eval(x), y=eval(y)) %>% mutate(ctag=!!ctag)
    tpl = tp %>% group_by(ctag) %>% nest() %>% ungroup() %>%
        mutate(fit = map(data, ~ lm(y~x, data=.x))) %>%
        mutate(tidied = map(fit, glance))%>%
        unnest(tidied) %>%
        mutate(lab = glue("adjusted R<sup>2</sup> = {number(adj.r.squared,accuracy=.01)}"))
    xlab = cmps %>% filter(cmp==x) %>% pull(cmp.l) %>% as.character()
    ylab = cmps %>% filter(cmp==y) %>% pull(cmp.l) %>% as.character()
    #
    ggplot(tp, aes(x=x,y=y)) +
        #geom_point(aes(x=x, y=y), alpha=.8, size=1) +
        #geom_text(data=tps, aes(x=bin.shift, y=n+100, label=n), vjust=0, size=2.5) +
        #geom_vline(xintercept=, linetype='dashed', size=.5) +
        geom_hex(bins=80) +
        geom_richtext(data=tpl, aes(x=.5,y=1,label=lab,hjust=.5,vjust=1), size=2.5) +
        scale_x_continuous(name=xlab,breaks=c(.1,.5),expand=expansion(mult=c(.02,.02))) +
        scale_y_continuous(name=ylab,breaks=c(.1,.5),expand=expansion(mult=c(.02,.02))) +
        scale_fill_viridis(name='density', direction=-1) +
        #scale_fill_manual(name='', values=pal_npg()(3)[c(3,2,1)]) +
        facet_wrap(ctag~., nrow=1) +
        otheme(legend.pos='bottom.right', legend.title=T, legend.dir='v',
               legend.box='h', legend.vjust=-.5, panel.spacing=.3,
               xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
               xgrid=T, ygrid=T) + o_margin(.3,.3,.3,.3)
    #}}}
}
p1 = plot_point_dens(tg, x='Br Bs', y='Ar As', ctag='genotype-wise', cmps)
p2 = plot_point_dens(tg, x='As Bs', y='Ar Br', ctag='tissue-wise', cmps)

fo = glue("{dirw}/10.shift.comp.pdf")
#ggsave(p, file=fo, width=8, height=4.5)
fo = glue("{dirf}/f2a.rds")
saveRDS(p1, fo)
fo = glue("{dirf}/f2b.rds")
saveRDS(p2, fo)
#}}}

#{{{ bar plot - f2c
tp = tg %>% filter(cmp %in% cmps5 | cmp %in% cmps4) %>%
    count(cmp, tag) %>%
    inner_join(cmps %>% select(cmp,cmp.l), by='cmp') %>%
    rename(tag1=cmp.l,tag2=tag)
p = cmp_proportion1(tp, xangle=30, ytext=T, legend.title='TC type:',
                    lab.size=2, fills = cols_shift) +
    o_margin(.1,.3,.1,.3) +
    theme(legend.position='none')
fo = glue("{dirf}/f2c.rds")
saveRDS(p, fo)
#}}}

#{{{ bar plot - sf0x
tp = tg %>% filter(cmp %in% cmps5 | cmp %in% cmps4) %>%
    count(cmp, bin.shift, tag) %>%
    inner_join(cmps %>% select(cmp,cmp.l), by='cmp')
tps = tp %>% group_by(cmp.l,bin.shift) %>% summarise(n=sum(n)) %>% ungroup()
p = ggplot(tp) +
    geom_col(aes(x=bin.shift, y=n, fill=tag), width=.8, alpha=1) +
    geom_text(data=tps, aes(x=bin.shift, y=n+100, label=n), vjust=0, size=2.5) +
    #geom_vline(xintercept=, linetype='dashed', size=.5) +
    scale_x_discrete(name='shifting score', expand=expansion(add=c(1,1))) +
    scale_y_continuous(name='number genes', expand=expansion(mult=c(.02,.05))) +
    scale_fill_manual(name='TC type', values=cols_shift) +
    facet_wrap(cmp.l~., nrow=2) +
    otheme(legend.pos='bottom.right', legend.title=T, legend.dir='v',
           legend.box='h', legend.vjust=-.5, panel.spacing=.3,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=F, ygrid=F) + o_margin(.3,.3,.3,.3) +
    theme(axis.text.x=element_text(angle=30, hjust=1,vjust=1))
fo = glue("{dirw}/10.shift.gene.pdf")
ggsave(p, file=fo, width=8, height=4)
fo = glue("{dirf}/sf0x.pdf")
#ggsave(p, file=fo, width=9, height=6)
#}}}

tp0 %>% group_by(cmp) %>% summarise(n_pass=sum(shift>=.5), n_total=n()) %>%
    ungroup() %>% mutate(pct_pass = percent(n_pass/n_total, accuracy=.1))
tp0 %>% filter(cmp %in% c("Br_Ar",'Bs_As'), shift>=.5) %>% count(i) %>% count(n)
tp0 %>% filter(cmp %in% c("As_Ar",'Bs_Br'), shift>=.5) %>% count(i) %>% count(n)
tpv = tp0 %>% filter(cmp %in% c("Br_Ar",'Bs_As'), shift>=.5) %>%
    mutate(cond = str_sub(cond1, 2, 2)) %>%
    group_by(i,gid) %>% summarise(conds = str_c(cond, collapse=',')) %>%
    ungroup()

tp0 %>% filter(shift<.05, padj>=.5) %>% count(gid) %>% count(n)
#}}}

#{{{ call stable & variable genes
tg1 = tg %>% filter(cmp %in% cmps5) %>%
    filter(tag == 'conserved') %>%
    count(i, gid) %>% filter(n==5) %>%
    select(gid) %>% mutate(tag='stable')
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
tgc %>% count(pnl)

ss = tg %>% select(gid,cmp,shift) %>%
    mutate(cmp = str_replace(cmp, 'C','M')) %>%
    mutate(cmp = str_replace(cmp, ' ','_')) %>%
    mutate(cmp = glue("SS_{cmp}")) %>%
    spread(cmp, shift)

to = tgc %>% filter(tag != 'all') %>%
    select(pnl, gid) %>%
    inner_join(ss, by='gid')
to %>% count(pnl)
fo = glue("{dird}/91_share/15.stable.variable.tsv")
write_tsv(to, fo)
#}}}

#{{{ characterize  stable TSSs
stb = tg %>% filter(cmp %in% cmps5) %>%
    filter(tag == 'conserved') %>%
    count(i, gid) %>% filter(n==5) %>%
    mutate(tag = 'stable') %>%
    inner_join(tgl, by='gid') %>% print(n=30)
stb %>% count(ttype)

#{{{ GO
tgo = read_go(src='Interproscan5')
tgo = read_go(src='arabidopsis')
tgo = read_go(src='uniprot.plants')
tgrp = tgo %>% select(gotype, grp=goid, gid, note=goname) %>%
    group_by(gotype) %>% nest() %>% rename(tgrp = data)

x0 = stb %>% rename(group = tag) %>%
    group_by(group) %>% summarise(gids = list(gid)) %>% ungroup() %>%
    crossing(tgrp) %>%
    mutate(x = map2(gids, tgrp, hyper_enrich))
x1 = x0 %>% select(group,gotype,x) %>% unnest(x) %>%
    arrange(group, pval.raw) %>%
    select(-grp) %>%
    group_by(group) %>% slice(1:15) %>% ungroup() %>% select(-group,-gotype) %>%
    #filter(group=='stable') %>%
    print(n=50)
x$note[[1]]
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

tp = gtss %>% select(gid,tag2=shape.s) %>%
    left_join(stb %>% select(gid,tag1=tag), by='gid') %>%
    replace_na(list(tag1='non-stable')) %>%
    count(tag1,tag2)
p1 = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = pal_npg()(5)) +
    o_margin(.1,.5,.1,.5)
fo = file.path(dirw, "32.shape.pdf")
ggsave(p1, file=fo, width=4, height=4)
#}}}

#{{{ tis-specificity
ti = read_tsv('~/projects/rnaseq/data/11_qc/Zmays_B73/rnc01/30.tis.expression.tsv.gz')
tiss = unique(ti$etag)[c(4,3,1,2)]
tp0 = ti %>% select(tag2=etag, gid) %>% mutate(tag2=factor(tag2, levels=tiss))

tpc = tp0 %>% count(tag2) %>% mutate(tag1='all genes')
tp = stb %>% rename(tag1=tag) %>%
    inner_join(tp0, by='gid') %>%
    count(tag1, tag2) %>%
    bind_rows(tpc) %>%
    mutate(tag1 = as_factor(tag1))
p1 = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = brewer.pal(5,'Set2')) +
    o_margin(.1,.1,.1,.1)
fo = file.path(dirw, "32.tis.pdf")
ggsave(p1, file=fo, width=4, height=4)
#}}}

#{{{ syn
t_syn = read_syn(gcfg) %>% arrange(gid, -ftype) %>% group_by(gid) %>% slice(1) %>% ungroup()
tp0 = stb %>% rename(ctag=tag)

tpc = t_syn %>% count(ftype) %>% mutate(tag1='all genes') %>% rename(tag2=ftype)
tp = tp0 %>%
    inner_join(t_syn, by='gid') %>%
    rename(tag1=ctag, tag2=ftype) %>% count(tag1, tag2) %>%
    bind_rows(tpc) %>%
    mutate(tag1 = as_factor(tag1))
p1 = cmp_proportion1(tp,ytext=T, oneline=T, ypos='right',legend.pos='none',
    fills = brewer.pal(5,'Set3')) +
    o_margin(.1,.1,.1,.1)
fo = file.path(dirw, "31.syn.pdf")
ggsave(p1, file=fo, width=4, height=4)

x = read_syn(gcfg, opt=2) %>%
    left_join(stb %>% select(maize1=gid,tag1=tag), by='maize1') %>%
    left_join(stb %>% select(maize2=gid,tag2=tag), by='maize2') %>%
    replace_na(list(tag1='non-stable', tag2='non-stable'))
x %>% count(ftype, tag1,tag2)
#}}}

#{{{ utr5 length
tg0 = gcfg$gene.loc %>% group_by(gid, ttype) %>%
    mutate(size = end - start + 1) %>%
    summarise(size.utr5 = sum(size[etype=='five_prime_UTR'])) %>%
    ungroup() %>% filter(ttype=='mRNA') %>% select(-ttype) %>%
    mutate(tag2 = ifelse(size.utr5 == 0, "UTR5 = 0", "UTR5 > 0"))
#tg0 = tg0 %>% filter(size.utr5 > 0)

tpc = tg0 %>% mutate(tag1='all genes') %>% select(gid,tag1,tag2)
tp = tg0 %>% mutate(tag1='stable') %>%
    inner_join(stb, by='gid') %>%
    replace_na(list(tag1 = 'all genes')) %>%
    bind_rows(tpc) %>%
    count(tag1, tag2) %>%
    mutate(tag1 = as_factor(tag1)) %>%
    mutate(tag2 = as_factor(tag2))
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


#{{{ characterize conserved/variable gene patterns in F1
cmps = c('Ar Br', 'As Bs')
tg1 = tg %>% filter(cmp %in% cmps) %>%
    filter(tag=='conserved') %>%
    count(gid) %>% filter(n==length(cmps)) %>%
    select(gid) %>% mutate(tag='conserved in root + shoot')
cmps = c('Ar As', 'Br Bs')
tg2 = tg %>% filter(cmp %in% cmps) %>%
    filter(tag=='conserved') %>%
    count(gid) %>% filter(n==length(cmps)) %>%
    select(gid) %>% mutate(tag='conserved in B73 + A632')
cmps = c('Ar Br', 'As Bs', 'Ar As', 'Br Bs')
tg3 = tg %>% filter(cmp %in% cmps) %>%
    filter(tag=='conserved') %>%
    count(gid) %>% filter(n==length(cmps)) %>%
    select(gid) %>% mutate(tag='conserved in 4 cmps')
tg4 = stb %>% select(gid) %>% mutate(tag='conserved in 5 cmps')
tgx = rbind(tg1,tg2,tg3,tg4) %>% mutate(tag = as_factor(tag))
tgxs = tgx %>% count(tag) %>%
    mutate(pnl = glue("{tag} ({n})")) %>%
    mutate(pnl = as_factor(pnl))
tgx = tgx %>% inner_join(tgxs, by='tag')

tp0 = tg %>% filter(cmp %in% c("Hs As", "Hs Bs", 'Hs Cs')) %>%
    select(gid, cmp, shift) %>%
    spread(cmp, shift) %>%
    rename(x=2, y=3, d=4) %>%
    filter(!is.na(x), !is.na(y)) %>%
    inner_join(tgx, by='gid')

#{{{
tp = tp0
tpl = tp %>% group_by(pnl) %>% nest() %>% ungroup() %>%
    mutate(fit = map(data, ~ lm(y~x, data=.x))) %>%
    mutate(tidied = map(fit, glance))%>%
    unnest(tidied) %>%
    mutate(lab = glue("adjusted R<sup>2</sup> = {number(adj.r.squared,accuracy=.01)}"))
xlab = 'shifting score: F1 shoot vs. A632 shoot'
ylab = 'shifting score: F1 shoot vs. B73 shoot'
#
p = ggplot(tp, aes(x=x,y=y)) +
    #geom_point(aes(x=x, y=y), alpha=.8, size=1) +
    #geom_text(data=tps, aes(x=bin.shift, y=n+100, label=n), vjust=0, size=2.5) +
    #geom_vline(xintercept=, linetype='dashed', size=.5) +
    geom_hex(bins=80) +
    geom_richtext(data=tpl, aes(x=.5,y=1,label=lab,hjust=.5,vjust=1), size=2.5) +
    scale_x_continuous(name=xlab,breaks=c(.1,.5),expand=expansion(mult=c(.02,.02)), limits=c(0,1)) +
    scale_y_continuous(name=ylab,breaks=c(.1,.5),expand=expansion(mult=c(.02,.02)), limits=c(0,1)) +
    scale_fill_viridis(name='density', direction=-1) +
    #scale_fill_manual(name='', values=pal_npg()(3)[c(3,2,1)]) +
    facet_wrap(pnl~., nrow=1) +
    otheme(legend.pos='bottom.right', legend.title=T, legend.dir='v',
           legend.box='h', legend.vjust=-.5, panel.spacing=.3,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) + o_margin(.3,.3,.3,.3)
#}}}
fo = glue("{dirw}/41.shift.pdf")
ggsave(p, file=fo, width=10, height=3)

tp1 = tp0 %>% filter(tag=='conserved in B73 + A632')
tgf = tg %>% filter(gid %in% tp1$gid, cond1=='As', cond2=='Bs') %>%
    select(gid, tag)
tp = tp1 %>% select(gid,x,y,d) %>%
    inner_join(tgf, by='gid')
tps = tp %>% count(tag) %>% arrange(tag) %>%
    mutate(pnl = glue("{tag} ({n})")) %>%
    mutate(pnl = as_factor(pnl))
tp = tp %>% inner_join(tps, by='tag')
#{{{
tpl = tp %>% group_by(pnl) %>% nest() %>% ungroup() %>%
    mutate(fit = map(data, ~ lm(y~x, data=.x))) %>%
    mutate(tidied = map(fit, glance))%>%
    unnest(tidied) %>%
    mutate(lab = glue("adjusted R<sup>2</sup> = {number(adj.r.squared,accuracy=.01)}"))
xlab = 'shifting score: F1 shoot vs. A632 shoot'
ylab = 'shifting score: F1 shoot vs. B73 shoot'
#
p = ggplot(tp, aes(x=x,y=y)) +
    geom_point(aes(x=x, y=y, color=d), size = .5) +
    #geom_point(aes(x=x, y=y), alpha=.8, size=1) +
    #geom_text(data=tps, aes(x=bin.shift, y=n+100, label=n), vjust=0, size=2.5) +
    #geom_vline(xintercept=, linetype='dashed', size=.5) +
    #geom_hex(bins=80) +
    geom_richtext(data=tpl, aes(x=.5,y=1,label=lab,hjust=.5,vjust=1), size=2.5) +
    scale_x_continuous(name=xlab,breaks=c(.1,.5),expand=expansion(mult=c(.02,.02)), limits=c(0,1)) +
    scale_y_continuous(name=ylab,breaks=c(.1,.5),expand=expansion(mult=c(.02,.02)), limits=c(0,1)) +
    scale_fill_viridis(name='density', direction=-1) +
    #scale_fill_manual(name='', values=pal_npg()(3)[c(3,2,1)]) +
    scale_color_viridis(name="shift_F1_merged", direction=-1) +
    facet_wrap(pnl~., nrow=1) +
    otheme(legend.pos='right', legend.title=T, legend.dir='v',
           legend.box='h', legend.vjust=-.5, panel.spacing=.3,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) + o_margin(.3,.3,.3,.3)
#}}}
fo = glue("{dirw}/41.shift.sub.pdf")
ggsave(p, file=fo, width=8, height=3)

#{{{ save for Nathan & Erich
fi = glue("{dird}/03_qc/02.tss.gtss.cond.rds")
r2 = readRDS(fi)
to0 = r2$gtss %>% select(gid,tpm.cond) %>%
    unnest(tpm.cond) %>% select(gid, tpm, cond) %>%
    filter(cond %in% c("Bs",'Br','As','Ar','Hs','Hr')) %>%
    mutate(cond = glue("tpm_{cond}")) %>%
    spread(cond,tpm)
to1 = tp %>% select(pnl,gid,SS_Hs_As=x,SS_Hs_Bs=y,SS_Hs_Ms=d,tag)
to2 = tg %>% filter(cmp == 'Ar As', gid %in% tp$gid) %>%
    select(gid, SS_Ar_As=shift)
to3 = tg %>% filter(cmp == 'Br Bs', gid %in% tp$gid) %>%
    select(gid, SS_Br_Bs=shift)
to4 = tg %>% filter(cmp == 'As Bs', gid %in% tp$gid) %>%
    select(gid, SS_As_Bs=shift)
lk = 'https://s3.msi.umn.edu/zhoup-igv/index.html?sessionURL=https://s3.msi.umn.edu/zhoup-igv-data/sessions/cage/{gid}.json'
to = to1 %>% inner_join(to0, by='gid') %>%
    inner_join(to2, by='gid') %>%
    inner_join(to3, by='gid') %>%
    inner_join(to4, by='gid') %>%
    select(panel=pnl,gid,starts_with("tpm"), SS_Ar_As,SS_Br_Bs,SS_As_Bs,
           SS_Hs_As,SS_Hs_Bs,SS_Hs_Ms,tag) %>%
    arrange(panel, gid)
links = to %>% glue_data(lk)
to = to %>% mutate(link = links)

fo = glue("{dird}/91_share/05.shift.3584.tsv")
write_tsv(to, fo)
#}}}

x = tp %>% filter(tag == 'conserved') %>%
    inner_join(tgl, by='gid') %>%
    filter(x>.4,y>.4) %>%
    print(n=20, width=Inf)
#}}}


#{{{ prepare meta table to share
fi = glue("{dirw}/04.rtss.rds")
rtss = readRDS(fi)
fi = glue('{dird}/03_qc/05.shift.rds')
r5 = readRDS(fi)
ss = r5$rtss; 

to1 = rtss %>% select(tidx, tpm.grp) %>% unnest(tpm.grp) %>%
    select(tidx, grp, tpm) %>%
    mutate(grp = fct_relabel(grp, ~ paste0("tpm.", .x))) %>%
    spread(grp, tpm)

cmps = ss$cmp[[1]] %>% mutate(cmp = glue("ss-{grp1}-{grp2}")) %>% pull(cmp)
to2 = ss %>% select(tidx=i, cmp) %>% unnest(cmp) %>%
    mutate(cmp = glue("ss-{grp1}-{grp2}")) %>%
    mutate(cmp = factor(cmp, levels=cmps)) %>%
    select(tidx, cmp, shift) %>% spread(cmp, shift)

to = rtss %>% select(-tpm.grp, -cidxs, -cnts) %>%
    inner_join(to1, by='tidx') %>%
    left_join(to2, by='tidx')

fo = glue("{dird}/91_share/01.rtss.tsv")
write_tsv(to, fo, na='')
#}}}

#{{{ [old] varaible TSS btw. B, A in F1
tis = 'root'
tv1 = tg %>% filter(cmp == 'Ar Br') %>% select(gid, tag, shift)
tv2 = tg %>% filter(cmp %in% c('Hr Br','Hr Ar','Hr Cr')) %>%
    mutate(shift = pmax(0, shift)) %>%
    select(gid, cmp, shift)
gids = tv2 %>% count(gid) %>% filter(n==3) %>% pull(gid)
tv2 = tv2 %>% filter(gid %in% gids) %>% spread(cmp, shift)
tv = tv1 %>% inner_join(tv2, by='gid') %>% rename(x=`Hr Br`, y=`Hr Ar`, d=`Hr Cr`)

tis = 'shoot'
tv1 = tg %>% filter(comp == 'Bs_As') %>% select(gid, tag, shift)
tv2 = tg %>% filter(comp %in% c('Hs_Bs','Hs_As','Hs_Cs')) %>%
    mutate(shift = pmax(0, shift)) %>%
    select(gid, comp, shift)
gids = tv2 %>% count(gid) %>% filter(n==3) %>% pull(gid)
tv2 = tv2 %>% filter(gid %in% gids) %>% spread(comp, shift)
tv = tv1 %>% inner_join(tv2, by='gid') %>% rename(x=Hs_Bs, y=Hs_As, d=Hs_Cs)

#{{{ plot
require(ggExtra)
plot_panel <- function(tag, ng, tp) {
    #{{{
    xtitle = 'Shifting score (F1 vs. B73)'
    ytitle = 'Shifting score (F1 vs. A632)'
    p1 = ggplot(tp) +
        geom_point(aes(x=x, y=y, color=d), size = .5) +
        #geom_abline(intercept = 0, slope = 1) +
        scale_x_continuous(name=xtitle,breaks=c(.1,.5), expand = c(.01, 0), limits = c(0,1)) +
        scale_y_continuous(name=ytitle,breaks=c(.1,.5), expand = c(.01, 0), limits = c(0,1)) +
        scale_color_viridis(direction=-1) +
        ggtitle(glue("B73 vs A632: {tag} ({ng})")) +
        otheme(legend.pos='none', legend.dir='v', legend.box='v',legend.title=F,
             margin = c(0,.1,.2,.1), xgrid=T, ygrid=T,
             xtick=T, ytick=T, xtitle=T, ytitle=T, xtext=T, ytext=T) +
        theme(plot.title=element_text(hjust=.5, size=10))
    if (tag != 'conserved') p1 = p1 + no_y_axis()
    ggMarginal(p1, type = 'histogram', margins = 'both', size = 3,
               #groupColour = T, groupFill = T,
               xparams = list(size=0), yparams = list(size=0))
    #}}}
}

to = tv %>% group_by(tag) %>% nest() %>% ungroup() %>% rename(tp=data) %>%
    mutate(ng = map_int(tp, nrow)) %>% arrange(tag) %>%
    mutate(p = pmap(list(tag,ng,tp), plot_panel))
#
fo = glue("{dirw}/14.{tis}.pdf")
ggarrange(to$p[[1]], to$p[[2]], to$p[[3]], nrow=1, ncol=3,
          labels = '', widths=c(1.1,1,1)) %>%
ggexport(filename=fo, width=10, height=4)
#}}}
#}}}

