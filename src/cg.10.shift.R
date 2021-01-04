source("functions.R")
dirw = glue("{dird}/10_shift")
setwd(dirw)

#{{{ read shifting scores - tg
fi = glue('{dird}/03_qc/05.shift.rds')
r5 = readRDS(fi)
tss = r5$tss; gtss = r5$gtss
conds4 = c('Br','Bs','Ar','As')
gene = gcfg$gene %>% mutate(loc=glue("{chrom}:{start}-{end}")) %>% select(gid,loc)
#
#tss %>% mutate(x=map_int(cmp,length)) %>% count(x, ncond4)
gtss %>% mutate(x=map_int(cmp,length)) %>% count(x, ncond4)
#
itvs = c(seq(0,1,by=.1))
tags = c("conserved",'somewhat similar','variable')
tg = gtss %>% select(i,gid, cmp) %>% unnest(cmp) %>%
    mutate(comp = str_c(cond1,cond2,sep="_")) %>%
    mutate(comp = as_factor(comp)) %>%
    arrange(comp, desc(shift)) %>%
    group_by(comp) %>%
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
    mutate(tag = ifelse(shift>.5, tags[3], ifelse(shift<=.1, tags[1], tags[2]))) %>%
    mutate(tag = factor(tag, levels=tags))
#}}}

#{{{ shifting score distribution
comps4 = c("Bs_Br","As_Ar", "Bs_As", "Br_Ar")
#{{{ scatter plot
require(hexbin)
tp = tg %>% filter(comp %in% comps4) %>% select(i,gid,comp,shift) %>%
    spread(comp, shift)
tp1 = tp %>% select(i, x=Bs_As,y=Br_Ar) %>% mutate(ctag='x=Bs_As y=Br_Ar')
tp2 = tp %>% select(i, x=Bs_Br,y=As_Ar) %>% mutate(ctag='x=Bs_Br y=As_Ar')
tp = rbind(tp1,tp2) %>% filter(!is.na(x), !is.na(y))
tpl = tp %>% group_by(ctag) %>% nest() %>% ungroup() %>%
    mutate(fit = map(data, ~ lm(y~x, data=.x))) %>%
    mutate(tidied = map(fit, glance))%>%
    unnest(tidied) %>%
    mutate(lab = glue("adjusted R<sup>2</sup> = {number(adj.r.squared,accuracy=.01)}"))

p = ggplot(tp, aes(x=x,y=y)) +
    #geom_point(aes(x=x, y=y), alpha=.8, size=1) +
    #geom_text(data=tps, aes(x=bin.shift, y=n+100, label=n), vjust=0, size=2.5) +
    #geom_vline(xintercept=, linetype='dashed', size=.5) +
    geom_hex(bins=80) +
    geom_richtext(data=tpl, aes(x=.5,y=1,label=lab,hjust=.5,vjust=1), size=2.5) +
    scale_x_continuous(name='',breaks=c(.1,.5),expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='',breaks=c(.1,.5),expand=expansion(mult=c(.02,.02))) +
    scale_fill_viridis(direction=-1) +
    #scale_fill_manual(name='', values=pal_npg()(3)[c(3,2,1)]) +
    facet_wrap(ctag~., nrow=1) +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           legend.box='h', legend.vjust=-.5, panel.spacing=.3,
           xtext=T, xtick=T, xtitle=F, ytitle=F, ytext=T, ytick=T,
           xgrid=T, ygrid=T) + o_margin(.3,.3,.3,.3)
fo = glue("{dirw}/10.shift.comp.pdf")
ggsave(p, file=fo, width=8, height=4.5)
#}}}

tp = tg %>%
    filter(comp %in% comps4) %>%
    count(comp, bin.shift, tag)
tps = tp %>% group_by(comp,bin.shift) %>% summarise(n=sum(n)) %>% ungroup()
#{{{ bar plot
p = ggplot(tp) +
    geom_col(aes(x=bin.shift, y=n, fill=tag), alpha=.8) +
    geom_text(data=tps, aes(x=bin.shift, y=n+100, label=n), vjust=0, size=2.5) +
    #geom_vline(xintercept=, linetype='dashed', size=.5) +
    scale_x_discrete(name='shifting score', expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name='num. genes', expand=expansion(mult=c(0,.05))) +
    scale_fill_manual(name='', values=pal_npg()(3)[c(3,2,1)]) +
    facet_wrap(comp~., nrow=2) +
    otheme(legend.pos='top.center.out', legend.title=T, legend.dir='h',
           legend.box='h', legend.vjust=-.5, panel.spacing=.3,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=F, ygrid=F) + o_margin(1,.3,.3,.3) +
    theme(axis.text.x=element_text(angle=30, hjust=1,vjust=1))
fo = glue("{dirw}/10.shift.gene.pdf")
ggsave(p, file=fo, width=6, height=6)
#}}}

tp0 %>% group_by(comp) %>% summarise(n_pass=sum(shift>=.5), n_total=n()) %>%
    ungroup() %>% mutate(pct_pass = percent(n_pass/n_total, accuracy=.1))
tp0 %>% filter(comp %in% c("Br_Ar",'Bs_As'), shift>=.5) %>% count(i) %>% count(n)
tp0 %>% filter(comp %in% c("As_Ar",'Bs_Br'), shift>=.5) %>% count(i) %>% count(n)
tpv = tp0 %>% filter(comp %in% c("Br_Ar",'Bs_As'), shift>=.5) %>%
    mutate(cond = str_sub(cond1, 2, 2)) %>%
    group_by(i,gid) %>% summarise(conds = str_c(cond, collapse=',')) %>%
    ungroup()

tp0 %>% filter(shift<.05, padj>=.5) %>% count(gid) %>% count(n)
#}}}

#{{{ characterize  stable TSSs
comps4 = c("Bs_Br","As_Ar", "Bs_As", "Br_Ar")
x = tg %>%
    filter(comp %in% comps4) %>%
    filter(tag == 'conserved') %>%
    count(i, gid) %>% filter(n==4) %>%
    inner_join(tgl, by='gid') %>% print(n=30)

#{{{ GO
tgo = read_go(src='Interproscan5')
tgo = read_go(src='arabidopsis')
tgo = read_go(src='uniprot.plants')
tgrp = tgo %>% select(gotype, grp=goid, gid, note=goname) %>%
    group_by(gotype) %>% nest() %>% rename(tgrp = data)

tg = tgl %>% rename(group=tag) %>%
    group_by(group) %>% summarise(gids = list(gid)) %>% ungroup() %>%
    crossing(tgrp) %>%
    mutate(x = map2(gids, tgrp, hyper_enrich))
x = tg %>% select(group,gotype,x) %>% unnest(x) %>%
    arrange(group, pval.raw) %>%
    select(-grp) %>%
    group_by(group) %>% slice(1:15) %>% ungroup() %>% select(-group,-gotype) %>%
    #filter(group=='stable') %>%
    print(n=50)
x$note[[1]]
#}}}

#{{{ tis-specificity
ti = read_tsv('~/projects/rnaseq/data/11_qc/Zmays_B73/rnc01/30.tis.expression.tsv.gz')
tiss = unique(ti$etag)[c(4,3,1,2)]
tp0 = ti %>% select(tag2=etag, gid) %>% mutate(tag2=factor(tag2, levels=tiss))

tpc = tp0 %>% count(tag2) %>% mutate(tag1='ctrl-bg')
tp = tgl %>% rename(tag1=tag) %>%
    inner_join(tp0, by='gid') %>%
    count(tag1, tag2) %>%
    bind_rows(tpc) %>%
    mutate(tag1 = as_factor(tag1))
p1 = cmp_proportion1(tp,xangle=0, oneline=T,legend.title='', barwidth=.8)
fo = file.path(dirw, "32.tis.pdf")
ggsave(p1, file=fo, width=4, height=4)
#}}}

#{{{ syn
t_syn = read_syn(gcfg) %>% arrange(gid, -ftype) %>% group_by(gid) %>% slice(1) %>% ungroup()
tp0 = tgl %>% rename(ctag=tag)

tpc = t_syn %>% count(ftype) %>% mutate(tag1='ctrl-bg') %>% rename(tag2=ftype)
tp = tp0 %>%
    inner_join(t_syn, by='gid') %>%
    rename(tag1=ctag, tag2=ftype) %>% count(tag1, tag2) %>%
    bind_rows(tpc) %>%
    mutate(tag1 = as_factor(tag1))
p1 = cmp_proportion1(tp,xangle=0, oneline=T,legend.title='synteny:') +
    otheme(legend.pos='top.center.out', legend.dir='v', xtext=T) +
    o_margin(4.5,.2,0,.2)
fo = file.path(dirw, "31.syn.pdf")
ggsave(p1, file=fo, width=4, height=5)
#}}}

#{{{ utr5 length
tg0 = gcfg$gene.loc %>% group_by(gid, ttype) %>%
    mutate(size = end - start + 1) %>%
    summarise(size.utr5 = sum(size[etype=='five_prime_UTR'])) %>%
    ungroup() %>% filter(ttype=='mRNA') %>% select(-ttype) %>%
    mutate(tag2 = ifelse(size.utr5 == 0, "UTR5 = 0", "UTR5 > 0"))
#tg0 = tg0 %>% filter(size.utr5 > 0)

tpc = tg0 %>% mutate(tag1='ctrl-bg') %>% select(gid,tag1,tag2)
tp = tg0 %>% mutate(tag1='stable') %>%
    inner_join(tgl, by='gid') %>%
    replace_na(list(tag1 = 'ctrl-bg')) %>%
    bind_rows(tpc) %>%
    count(tag1, tag2) %>%
    mutate(tag1 = as_factor(tag1)) %>%
    mutate(tag2 = as_factor(tag2))
p1 = cmp_proportion1(tp,xangle=0, oneline=T,legend.title='', barwidth=.8)
fo = file.path(dirw, "33.utr5.1.pdf")
ggsave(p1, file=fo, width=3, height=4)

tpc = tg0 %>% mutate(tag1='ctrl-bg') %>% select(gid,tag1,size.utr5)
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



#{{{ varaible TSS btw. B, A in F1
tis = 'root'
tv1 = tg %>% filter(comp == 'Br_Ar') %>% select(gid, tag, shift)
tv2 = tg %>% filter(comp %in% c('Hr_Br','Hr_Ar','Hr_Cr')) %>%
    mutate(shift = pmax(0, shift)) %>%
    select(gid, comp, shift)
gids = tv2 %>% count(gid) %>% filter(n==3) %>% pull(gid)
tv2 = tv2 %>% filter(gid %in% gids) %>% spread(comp, shift)
tv = tv1 %>% inner_join(tv2, by='gid') %>% rename(x=Hr_Br, y=Hr_Ar, d=Hr_Cr)

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

x = tv %>% inner_join(tgl, by='gid') %>% select(gid,shift,x,y,d,loc)

x %>% filter(shift < .1, x < .2, y < .2) %>%  print(n=30)

x %>% filter(shift>.1,shift < .3, x >=.5, y >= .5) %>% print(n=30)

x %>% filter(shift > .6, d < .1) %>% print(n=30)



