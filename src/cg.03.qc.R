require(CAGEfightR)
require(GenomicFeatures)
source("functions.R")
dirw = glue("{dird}/03_qc")
setwd(dirw)

#{{{ mapping rates
yid = 'cg20a'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

tp = res$bamstat %>% separate(sid, c('sid','suf'), sep='[\\.]') %>%
    inner_join(thf, by=c('sid'='SampleID')) %>%
    mutate(map.rate = unpair_map / unpair) %>%
    select(gt=Genotype, tis=Tissue, score=map.rate) %>%
    filter(gt %in% c("B73",'A632'), tis %in% c("root",'shoot'))
tps = tp %>%
    group_by(gt,tis) %>%
    summarise(avg=mean(score), std=sd(score)) %>%
    ungroup() %>% print(n=40)

p = ggplot(tp, aes(x=gt, y=score)) +
    geom_boxplot(aes(group=gt), size=.5, width=.7, outlier.shape=NA) +
    geom_jitter(aes(color=gt), width=.3, size=2, alpha=.9) +
    scale_x_discrete(expand=expansion(mult=c(.5,.5))) +
    scale_y_continuous(name='Mapping Rate', expand=expansion(mult=c(.02,.02)),limits=c(0,1)) +
    scale_color_aaas() +
    facet_grid(.~tis, scale='free') +
    #ggtitle(tp$tit[1]) +
    otheme(ytitle=T,ytext=T,ytick=T, xtick=T, xtext=T,
           legend.pos='none', panel.spacing=.2, strip.compact=F) +
    theme(plot.title=element_text(size=8))
fo = glue("{dirw}/07.mapping.rate.pdf")
ggsave(p, file=fo, width=4, height=4)
fo = glue("{dirf}/sf03.pdf")
ggsave(p, file=fo, width=4, height=4)
#}}}

#{{{ call TSS & gTSS
txdb = load_txdb("Zmays_B73")
dsg = thf %>% rename(Name=SampleID)
diri = glue("{dird}/tracks/raw_bw")

SE_ctss = read_cage_bws(diri, dsg, minSupport=2)
#SE_tss = make_tss(SE_ctss, unexpressed=.5, minSamples=2)
SE_tss = make_tss(SE_ctss, txdb, unexpressed=1, minSamples=3)

# characterize condition-wise TPM
x = thf %>% rename(cond=cond.s) %>%
    group_by(cond) %>% summarise(sids = list(SampleID)) %>% ungroup() %>%
    mutate(x = map(sids, subset_samples, ex=SE_tss)) %>%
    select(-sids) %>% unnest(x) %>%
    group_by(i) %>% nest() %>% rename(tpm.cond=data)

tss = rowRanges(SE_tss) %>% as_tibble() %>%
    mutate(i = 1:n()) %>% mutate(tpm = score/ncol(assay(SE_tss))) %>%
    inner_join(x, by='i') %>%
    mutate(tidx = glue("t{i}")) %>%
    select(tidx, chrom=seqnames,start,end,srd=strand,tpm,domi=thick.start,
           IQR, support, tid=txID, peakType=peakTxType, gid=geneID, tpm.cond)
tss %>% count(peakType)
tss %>% select(tpm.cond) %>% unnest(tpm.cond) %>% group_by(cond) %>% summarise(n_expressed=sum(tpm>=1))

#{{{ create gene-level TSS clusters
#{{{ filter TSS and create GRanges object
x = tss %>% filter(peakType %in% c("promoter","proximal",'fiveUTR'))
x %>% mutate(size=end-start+1) %>% count(size) %>% arrange(desc(size))
x1 = x %>% group_by(gid) %>%
    summarise(chrom=chrom[1], start=min(start), end=max(end), srd=srd[1],
              n_tss = n(), tidx = list(tidx)) %>%
    ungroup() %>%
    arrange(chrom,start,end)
x1 %>% mutate(size=end-start+1) %>% count(size) %>% arrange(desc(size))
gr = with(x1, GRanges(seqnames=chrom,ranges=IRanges(start,end=end), strand=srd))
isDisjoint(gr)

#{{{ fix overlapping genes
idxs = which(disjointBins(gr) != 1)
x1[(idxs-1):idxs,]
gr = gr[-idxs]
isDisjoint(gr)
#}}}
seqinfo(gr) = seqinfo(SE_ctss)
#}}}
SE_gtss = SE_ctss %>% quantifyClusters(clusters=gr, inputAssay='counts') %>%
    calcTPM(totalTags='totalTags') %>%
    calcPooled() %>%# assignGeneID(geneModels=txdb, outputColumn='geneID')
    calcShape(pooled=SE_ctss, outputColumn='IQR', shapeFunction=shapeIQR,
        lower=.1, upper=.9)

x = thf %>% rename(cond=cond.s) %>%
    group_by(cond) %>% summarise(sids = list(SampleID)) %>% ungroup() %>%
    mutate(x = map(sids, subset_samples, ex=SE_gtss)) %>%
    select(-sids) %>% unnest(x) %>%
    group_by(i) %>% nest() %>% rename(tpm.cond=data)

gtss = rowRanges(SE_gtss) %>% as_tibble() %>%
    mutate(i = 1:n()) %>% mutate(tpm = score/ncol(assay(SE_gtss))) %>%
    select(i, chrom=seqnames,start,end,srd=strand,tpm,IQR) %>%
    inner_join(x1, by=c('chrom','start','end','srd')) %>%
    inner_join(x, by='i') %>%
    rename(gidx=i) %>% mutate(gidx = glue("g{gidx}"))
gtss %>% count(n_tss)
gtss %>% count(n_tss > 1)
gtss %>% select(tpm.cond) %>% unnest(tpm.cond) %>% group_by(cond) %>% summarise(n_expressed=sum(tpm>=1))
#}}}

r1 = list(SE_ctss=SE_ctss, SE_tss=SE_tss, tss=tss,
          SE_gtss=SE_gtss, gtss=gtss)
fo = glue("{dirw}/01.tss.gtss.rds")
saveRDS(r1, fo)
#}}}

#{{{ subset to B73 TSS & gTSS
fi = glue("{dirw}/01.tss.gtss.rds")
r1 = readRDS(fi)
SE_ctss=r1$SE_ctss; SE_tss=r1$SE_tss; tss=r1$tss; SE_gtss=r1$SE_gtss; gtss=r1$gtss
#
tcond = thf %>% rename(cond=cond.s,gt=Genotype,tis=Tissue) %>% distinct(cond,gt,tis)

ti1 = tss %>% select(i=tidx,tpm.cond) %>% unnest(tpm.cond) %>%
    mutate(s = ifelse(support/nrep >= .2, "E", "x")) %>%
    inner_join(tcond, by='cond') %>%
    select(i,cond,gt, tis, s) %>% mutate(tag = 'tss')
ti2 = gtss %>% select(i=gidx,tpm.cond) %>% unnest(tpm.cond) %>%
    mutate(s = ifelse(support/nrep >= .2, "E", "x")) %>%
    inner_join(tcond, by='cond') %>%
    select(i,cond,gt, tis, s) %>% mutate(tag = 'gene')
t_exp = ti1 %>% rbind(ti2) %>% mutate(tag = factor(tag, levels=c('tss','gene')))

tb = t_exp %>% filter(gt == 'B73', s=='E') %>% distinct(i, tag)
tb %>% count(tag)
tss.b = tss %>% filter(tidx %in% tb$i)
gtss.b = gtss %>% filter(gidx %in% tb$i)

r2 = list(tss=tss,gtss=gtss,tss.b=tss.b,gtss.b=gtss.b,tb=tb)
fo = glue("{dirw}/02.tss.gtss.cond.rds")
saveRDS(r2, fo)

# B73 stats
nrow(gtss.b)
gtss.b %>% count(n_tss == 1)
tss.b %>% count(peakType== 'intergenic')
tss.b %>% count(peakType)
gtss.b %>% inner_join(gcfg$gene %>% select(gid,ttype), by='gid') %>%count(ttype)

#{{{ B73 only stats
ti1 %>% filter(s=='E') %>% count(cond)
ti2 = ti1 %>% spread(cond, s) %>%
    rename(Broot=Br, Bshoot=Bs) %>%
    count(Broot,Bshoot,Bstem,Bhusk) %>% print(n=20)

ti2a = ti2 %>% slice(1:12)
ti2b = ti2 %>% filter(Broot=='x',Bshoot=='x',Bstem=='E'|Bhusk=='E')
ti2c = ti2 %>% slice(1:15)
cat(glue("{sum(ti2$n)} total TSSs/gTSSs\n"))
cat(glue("{sum(ti2a$n)} root/shoot TSSs/gTSSs\n"))
cat(glue("{sum(ti2b$n)} stem/husk-specific TSSs/gTSSs\n"))
cat(glue("{sum(ti2c$n)} B73 expressed TSSs/gTSSs\n"))
cat(glue("{ti2 %>% slice(16) %>% pluck('n',1)} total B73 TSSs/gTSSs\n"))
#}}}
#}}}

#{{{ shape characterization
fi = glue("{dirw}/02.tss.gtss.cond.rds")
r2 = readRDS(fi)
tss=r2$tss.b; gtss=r2$gtss.b
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

#{{{ different feature types
tp = tss %>% rename(size=iqr) %>% count(peakType, size, shape)
p1 = tp %>% ggplot(aes(x=size, y=n, fill=shape)) +
    geom_col(alpha=1) +
    #scale_x_continuous(breaks=tpx$iqr2, labels=tpx$iqr, expand=expansion(mult=c(.01,.01))) +
    scale_x_continuous(name='10-90% IQR', expand=expansion(mult=c(.02,.02)), limits=c(NA,100)) +
    scale_y_continuous(name='Frequency', expand=expansion(mult=c(0,.03))) +
    scale_fill_npg(name='peakType') +
    facet_wrap(peakType~., nrow=3, scale='free_y') +
    otheme(legend.pos='top.right', legend.spacing.x=.5, legend.spacing.y=.5,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) +
    guides(color=F)
#}}}
fo = file.path(dirw, '11.ftype.panel.pdf')
ggsave(p1, file=fo, width=8, height=8)

#{{{ f1c & sf02b
tp = tss %>% filter(peakType != 'exon') %>%
    rename(tag2=shape) %>% rename(tag1=peakType) %>% count(tag1, tag2)
p = cmp_cnt1(tp, ytext=F, ypos='right', legend.title='shape:') +
    o_margin(.1,.3,.1,.3) +
    theme(legend.position='none')
fo = glue("{dirw}/11.ftype.iqr.cnt.pdf")
ggsave(p, file=fo, width=5, height=5)
fo = glue("{dirf}/f1c.rds")
saveRDS(p, fo)

p = cmp_proportion1(tp, ytext=T, ypos='right', legend.title='shape:') +
    o_margin(.1,.3,.1,.3) +
    theme(legend.position='none')
fo = glue("{dirw}/11.ftype.iqr.prop.pdf")
ggsave(p, file=fo, width=5, height=5)
fo = glue("{dirf}/sf02b.rds")
saveRDS(p, fo)
#}}}

#{{{ piechart of gene TSS shapes - sf02c
tp = gtss %>% count(shape.s) %>% rename(tag=shape.s)
p = pie1(tp, lab.size=2.5)
fo = glue("{dirw}/13.shape.pie.pdf")
ggsave(p, file=fo, width=4, height=4)
fo = glue("{dirf}/sf02c.rds")
saveRDS(p, fo)
#}}}

#{{{ expression level
tp = gtss %>% select(shape.s, tpm.cond) %>% rename(tag=shape.s) %>%
    unnest(tpm.cond) %>% inner_join(thfs, by=c('cond'='cond.s')) %>%
    filter(Genotype=='B73')
p = ggplot(tp) +
    geom_density(aes(x=tpm, color=tag,fill=tag), alpha=.5) +
    scale_x_continuous(name='TPM',limits=c(0,30), expand=expansion(mult=c(.03,.03))) +
    scale_y_continuous(expand=expansion(mult=c(.00,.03))) +
    scale_fill_npg(name='shape') +
    scale_color_npg(name='shape') +
    facet_wrap(~Tissue, nrow=2, scale='free') +
    otheme(xtitle=T,xtext=T,xtick=T, legend.pos='top.right')
fo = glue("{dirw}/13.shape.tpm.pdf")
ggsave(p, file=fo, width=6, height=6)
#fo = glue("{dirf}/sf02c.rds")
#saveRDS(p, fo)
#}}}

#{{{ GO
tgo = read_go(src='Interproscan5')
tgo = read_go(src='arabidopsis')
tgo = read_go(src='uniprot.plants')
tgrp = tgo %>% select(gotype, grp=goid, gid, note=goname) %>%
    group_by(gotype) %>% nest() %>% rename(tgrp = data)

tg = tsh %>% rename(group=tag) %>%
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
ti$etag = factor(ti$etag, levels=tiss)

tpc = ti %>% count(etag) %>% mutate(tag='all genes')
tp = gtss %>% rename(tag=shape.s) %>%
    inner_join(ti, by='gid') %>%
    count(tag, etag) %>%
    #bind_rows(tpc) %>%
    mutate(tag = as_factor(tag)) %>%
    rename(tag1=etag, tag2=tag)
p1 = cmp_proportion1(tp, oneline=T, ytext=T,ypos='right',legend.pos='none') +
    theme(axis.text.x = element_text(angle=10, hjust=.5,vjust=1, size=8))
fo = glue("{dirw}/13.shape.tis.pdf")
ggsave(p1, file=fo, width=4, height=4)
#}}}

#{{{ syn
t_syn = read_syn(gcfg) %>% arrange(gid, -ftype) %>% group_by(gid) %>% slice(1) %>% ungroup()
tp0 = gtss %>% rename(tag=shape.s)

#tpc = t_syn %>% count(ftype) %>% mutate(tag1='ctrl-bg') %>% rename(tag2=ftype)
tp = tp0 %>%
    inner_join(t_syn, by='gid') %>%
    count(ftype,tag) %>%
    rename(tag1=ftype, tag2=tag)
p1 = cmp_proportion1(tp, oneline=T, ytext=T,ypos='right',legend.pos='none') +
    theme(axis.text.x = element_text(angle=10, hjust=.5,vjust=1, size=8))
fo = glue("{dirw}/13.shape.syn.pdf")
ggsave(p1, file=fo, width=5, height=4)
#}}}

#{{{ utr5 length
tg0 = gcfg$gene.loc %>% group_by(gid, ttype) %>%
    mutate(size = end - start + 1) %>%
    summarise(size.utr5 = sum(size[etype=='five_prime_UTR'])) %>%
    ungroup() %>% filter(ttype=='mRNA') %>% select(-ttype) %>%
    mutate(tag2 = ifelse(size.utr5 == 0, "UTR5 = 0", "UTR5 > 0"))
#tg0 = tg0 %>% filter(size.utr5 > 0)

tpc = tg0 %>% mutate(tag1='ctrl-bg') %>% select(gid,tag1,size.utr5,tag2)
tp0 = tsh %>% select(gid,tag1=tag) %>%
    inner_join(tg0, by='gid') %>%
    bind_rows(tpc) %>%
    mutate(tag1 = factor(tag1, levels=c(shapes,'ctrl-bg'))) %>%
    mutate(tag2 = as_factor(tag2))

tp = tp0 %>% count(tag1, tag2)
p1 = cmp_proportion1(tp,xangle=0, oneline=T,legend.title='', barwidth=.8)
fo = file.path(dirw, "13.utr5.1.pdf")
ggsave(p1, file=fo, width=3, height=4)

#{{{ boxplot
tps = tp0 %>% filter(size.utr5>0) %>% rename(ctag=tag1, score=size.utr5) %>% group_by(ctag) %>%
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
fo = file.path(dirw, '13.utr5.2.pdf')
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

#{{{ cross-species TSS shape characterization
#{{{ merge 4 species
fi = glue("{dirw}/02.tss.gtss.cond.rds")
tss1 = readRDS(fi)$tss.b %>% select(-tpm.cond) %>%
    rename(i=tidx,strand=srd,pTPM=tpm) %>% mutate(org='Maize')
fi = glue("{dird}/90_refs/Kurihara2018/10.tss.rds")
tss2 = readRDS(fi)$tss %>% mutate(org='Arabidopsis')
fi = glue("{dird}/90_refs/Phantom5_human/10.tss.rds")
tss3 = readRDS(fi)$tss %>% mutate(org='Human')
fi = glue("{dird}/90_refs/Phantom5_mouse/10.tss.rds")
tss4 = readRDS(fi)$tss %>% mutate(org='Mouse')
#
orgs= c('Human','Mouse','Arabidopsis','Maize')
tt = rbind(tss1,tss2,tss3,tss4) %>% rename(iqr=IQR) %>%
    mutate(shape = map_chr(iqr, iqr2shape)) %>%
    mutate(shape.s = map_chr(iqr, iqr2shape, opt='s')) %>%
    mutate(shape = factor(shape, levels=shapes)) %>%
    mutate(shape.s = factor(shape.s, levels=shapess)) %>%
    group_by(org) %>%
    mutate(TPM = pTPM / sum(pTPM) * 1e6) %>% ungroup() %>%
    mutate(org=factor(org, levels=orgs))
tts = tt %>% count(org) %>% arrange(org) %>%
    mutate(pnl = glue("{org} ({number(n,accuracy=1)})")) %>%
    mutate(pnl = as_factor(pnl))
tt = tt %>% inner_join(tts, by='org')
#}}}

#{{{ iqr dist in 4 species - sf02a
tp = tt %>% rename(size=iqr) %>% count(pnl, size, shape)
p = tp %>% ggplot(aes(x=size, y=n, fill=shape)) +
    geom_col(alpha=1) +
    #scale_x_continuous(breaks=tpx$iqr2, labels=tpx$iqr, expand=expansion(mult=c(.01,.01))) +
    scale_x_continuous(name='10-90% IQR', expand=expansion(mult=c(.02,.02)), limits=c(NA,100)) +
    scale_y_continuous(name='Frequency', expand=expansion(mult=c(0,.03))) +
    scale_fill_npg(name='shape') +
    facet_wrap(pnl~., nrow=2, scale='free_y') +
    otheme(legend.pos='top.right', legend.spacing.x=.5, legend.spacing.y=.5,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) +
    guides(color=F)
#}}}
fo = file.path(dirw, '12.cross.panel.pdf')
ggsave(p, file=fo, width=8, height=6)
fo = glue("{dirf}/sf02a.rds")
saveRDS(p, fo)

#{{{ iqr dist in maize - f1b
tp = tt %>% rename(size=iqr) %>% filter(org=='Maize') %>%
    count(pnl, size, shape)
p = tp %>% ggplot(aes(x=size, y=n, fill=shape)) +
    geom_col(alpha=1) +
    #scale_x_continuous(breaks=tpx$iqr2, labels=tpx$iqr, expand=expansion(mult=c(.01,.01))) +
    scale_x_continuous(name='10-90% IQR', expand=expansion(mult=c(.02,.02)), limits=c(NA,100)) +
    scale_y_continuous(name='Frequency', expand=expansion(mult=c(0,.03))) +
    scale_fill_npg(name='shape') +
    facet_wrap(pnl~., nrow=2, scale='free_y') +
    otheme(legend.pos='top.right', legend.spacing.x=.5, legend.spacing.y=.5,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) +
    guides(color=F)
#}}}
fo = glue("{dirf}/f1b.rds")
saveRDS(p, fo)

#{{{ f1d
tp = tt %>%
    rename(tag2=shape, tag1=org) %>% count(tag1, tag2)
p = cmp_proportion1(tp, ytext=F, oneline=T, legend.title='shape:') +
    o_margin(.3,.3,.01,.3) +
    theme(legend.position='none')
#}}}
fo = glue("{dirw}/12.cross.prop.pdf")
ggsave(p, file=fo, width=5, height=5)
fo = glue("{dirf}/f1d.rds")
saveRDS(p, fo)

#{{{ shape of different annotation features - sf02b
tp = tt %>% count(org, peakType, shape.s) %>%
    rename(pnl=org, tag1=peakType, tag2=shape.s)
p = cmp_proportion(tp, xangle=40, oneline=F, legend.title='shape:', lab.size=1.5,
    strip.compact=T, expand.x=c(.1,.1), expand.y=c(.01,.08), nc=2, alph=.7) +
    o_margin(.5,.5,.2,.8)
#tp = tt %>%
    #rename(pnl=org, tag1=shape, tag2=peakType) %>% count(pnl, tag1, tag2)
#p1 = cmp_proportion(tp, xangle=30, oneline=T, legend.title='', pal='simpsons',
    #margin.l=.2, nc=4) + o_margin(1,.5,.2,.5) +
    #guides(fill=guide_legend(nrow=1))
#}}}
fo = glue("{dirw}/12.cross.ptype.pdf")
ggsave(p, file=fo, width=6, height=6)
fo = glue("{dirf}/sf02b.rds")
saveRDS(p, fo)

#{{{ shape TPM
tp = tt %>% rename(pnl=org,ctag = shape) %>% mutate(score = log10(TPM))
tps = tp %>% group_by(pnl, ctag) %>%
  summarise(y = median(score),
            y25 = quantile(score,.25), y75=quantile(score, .75),
            y5 = quantile(score,.05), y95=quantile(score,.95)) %>% ungroup()
pv = ggplot(tps) +
    #geom_violin(aes(x=ctag, y=score), trim=T, alpha=.8) +
    #geom_boxplot(aes(x=ctag, y=score), outlier.shape=NA, width=.3) +
    geom_boxplot(aes(x=ctag, ymin=y5,lower=y25, middle=y, upper=y75, ymax=y95, fill=pnl),
                 stat='identity', width=.3, position=position_dodge(.5)) +
    scale_x_discrete(expand=expansion(mult=c(.2,.2))) +
    scale_y_continuous(name='log10(TPM)', expand=expansion(mult=c(.05,.05))) +
    scale_fill_manual(values=pal_jco()(8)) +
    scale_color_aaas(name='') +
    otheme(legend.pos='top.left', strip.style='white',
           xtext=T, xtick=T, ytext=T, ytick=T, ytitle=T, ygrid=T) +
    theme(axis.text.x = element_text(angle = 20, hjust=.8, vjust=1.1))
fo = file.path(dirw, '11.shape.tpm.pdf')
ggsave(pv, filename=fo, width=4, height=4)
#}}}

x3 = x2 %>% select(shape,freq, meta) %>% unnest(meta)
x3 %>% count(peakType, shape) %>% spread(shape, n)
x3 %>% count(peakType, shape) %>% group_by(shape) %>% mutate(p=n/sum(n)) %>% ungroup() %>% select(-n) %>% spread(shape, p)
x3 %>% count(peakType, shape) %>% group_by(peakType) %>% mutate(p=n/sum(n)) %>% ungroup() %>% select(-n) %>% spread(shape, p)
x3 %>% count(peakType, freq) %>% spread(freq, n)
x3 %>% count(peakType, freq) %>% group_by(freq) %>% mutate(p=n/sum(n)) %>% ungroup() %>% select(-n) %>% spread(freq, p)
#}}}

#{{{ sample clustering t-SNE - sf03b
fi = glue("{dirw}/02.tss.gtss.cond.rds")
r2 = readRDS(fi)
tss=r2$tss; gtss=r2$gtss
tcond = thf %>% rename(cond=cond.s,gt=Genotype,tis=Tissue) %>% distinct(cond,gt,tis)
#
th = thf %>% mutate(grp = glue("{Genotype} {Tissue}")) %>%
    group_by(grp) %>% mutate(idx = 1:n()) %>% ungroup() %>%
    mutate(clab=ifelse(idx==1, grp, ''))

res = rnaseq_cpm_raw(yid)
#tm = assay(SE_tss, 'TPM') %>% as_tibble() %>% mutate(gid = glue("g{1:n()}")) %>%
    #gather(SampleID, value, -gid) %>%
    #mutate(value=asinh(value))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

p = plot_tsne(tm,th,pct.exp=.6,perp=3,iter=1500, seed=2,
    var.shape='Genotype', var.col='Tissue', var.ellipse='grp', var.lab='clab',
    shapes=c(15,16,4), legend.pos='top.right', legend.dir='v', pal.col='aaas') +
    theme(legend.box='horizontal', legend.margin=margin(.5,.5,0,0,'lines')) +
    guides(shape=guide_legend(order=1))
fp = glue('{dirw}/10.tsne.pdf')
ggsave(p, file=fp, width=5, height=5)
fp = glue('{dirf}/sf03b.pdf')
ggsave(p, file=fp, width=5, height=5)
#}}}

#{{{ TSS support in tissue/genotypes
fi = glue("{dirw}/02.tss.gtss.cond.rds")
r2 = readRDS(fi)
tss=r2$tss; gtss=r2$gtss
tcond = thf %>% rename(cond=cond.s,gt=Genotype,tis=Tissue) %>% distinct(cond,gt,tis)

ti1 = tss %>% select(i=tidx,tpm.cond) %>% unnest(tpm.cond) %>%
    mutate(s = ifelse(support/nrep >= .2, "E", "x")) %>%
    inner_join(tcond, by='cond') %>%
    select(i,cond,gt, tis, s) %>% mutate(tag = 'tss')
ti2 = gtss %>% select(i=gidx,tpm.cond) %>% unnest(tpm.cond) %>%
    mutate(s = ifelse(support/nrep >= .2, "E", "x")) %>%
    inner_join(tcond, by='cond') %>%
    select(i,cond,gt, tis, s) %>% mutate(tag = 'gene')
ti = ti1 %>% rbind(ti2) %>% mutate(tag = factor(tag, levels=c('tss','gene')))

#{{{ st01 stats
tis = ti %>% filter(cond=='Bs') %>% select(-cond,-gt,-tis) %>% rename(s0=s)
to1 = ti %>% inner_join(tis, by=c("i",'tag')) %>%
    filter(s == 'E') %>%
    group_by(tag, cond,gt, tis) %>%
    summarise(n=n(), n1 = sum(s0=='x')) %>%
    ungroup() %>% pivot_wider(names_from=tag,values_from=c(n,n1)) %>%
    rename(n_tss_non_shoot=n1_tss, n1_gene_non_shoot=n1_gene)
#
to2 = ti %>% inner_join(tis, by=c("i",'tag')) %>%
    filter(s == 'E', s0 == 'x', tag=='tss') %>% select(tidx=i, cond,gt,tis)
g2t = gtss %>% select(gidx,tidx) %>% unnest(tidx)
to2b = ti %>% filter(tag=='gene', cond=='Bs', s=='E') %>%
    select(gidx=i) %>% inner_join(g2t, by='gidx')
to3 = to2 %>% inner_join(to2b, by='tidx') %>% count(cond,gt,tis) %>%
    rename(n_alt_tss = n)
#
to2 = to1 %>% left_join(to3, by=c('cond','gt','tis')) %>% select(-cond) %>%
    replace_na(list(n_alt_tss = 0))

x = to2 %>%
    rename(Genotype=1,Tissue=2,`# TCs`=3,`# Genes w. TCs`=4,
        `# novel TCs`=5, `# Genes w. novel TCs`=6,
        `# Genes w. Alternate TCs`=7) %>%
    kbl(format='latex', escape=T, longtable=F, booktabs=T, linesep="",
        format.args = list(big.mark = ",")) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
        full_width=F, font_size = 9, position='left')
fo = file.path(dirf, 'st1.rds')
saveRDS(x, file=fo)
fo = file.path(dirf, 'st1.pdf')
save_kable(x,file=fo)
#}}}

#{{{ ms stats
ti %>% filter(gt == 'B73', s=='E') %>% distinct(i,tag) %>% count(tag)
tib = ti %>% filter(gt == 'B73', s=='E', tag=='gene')
tib %>% count(cond) %>% pull(n) %>% mean()
tib %>% distinct(i)

tib = ti %>% filter(gt == 'B73', s=='E')
tia = ti %>% filter(gt == 'A632', s=='E') %>% distinct(i, tag) %>%
    filter(! i %in% tib$i) %>% count(tag) %>% print(n=4)
#}}}
#}}}


#{{{ merge replicates & make counts table for TSSs & gTSSs
fi = glue("{dirw}/01.tss.gtss.rds")
r1 = readRDS(fi)
SE_ctss=r1$SE_ctss; SE_tss=r1$SE_tss; tss=r1$tss; SE_gtss=r1$SE_gtss; gtss=r1$gtss
tcond = thf %>% rename(cond=cond.s,gt=Genotype,tis=Tissue) %>% distinct(cond,gt,tis)

#{{{ make replicate-merged ctss object
rowsum1 <- function(sids, ctss) rowSums(assay(ctss[,sids],'counts'))
x = thf %>% rename(cond=cond.s) %>% group_by(cond) %>%
    summarise(sids=list(SampleID)) %>% ungroup() %>%
    mutate(cnts = map(sids, rowsum1, ctss=SE_ctss))
ns = nrow(thfs)
xm = matrix(x %>% select(cnts) %>% unnest(cnts) %>% pull(cnts), ncol=ns, byrow=F)
xh = thfs %>% rename(Name=cond.s)
colnames(xm) = xh$Name
#
SE_ctss_m = SummarizedExperiment(assays=SimpleList(counts=xm),
                         rowRanges=rowRanges(SE_ctss)[,-c(1,2)],
                         colData=xh) %>%
    calcTotalTags(inputAssay="counts", outputColumn="totalTags")
ctss_m = rowRanges(SE_ctss)[,-c(1,2)] %>% as.data.frame() %>% as_tibble() %>%
    mutate(cidx=1:n()) %>% select(cidx,chrom=1,pos,strand)
#}}}

#{{{ TSS - merge replicates
ti = findOverlaps(SE_ctss_m, rowRanges(SE_tss)) %>% as_tibble() %>%
    select(cidx=1, i=2)
tic = assay(SE_ctss_m, 'counts') %>% as_tibble() %>%
    mutate(cidx = 1:n()) %>%
    inner_join(ti, by='cidx') %>%
    gather(cond, cnt, -i, -cidx) %>%
    arrange(i, cond, cidx) %>%
    group_by(i, cond) %>% summarise(cnts=list(cnt), cidxs=list(cidx)) %>%
    ungroup() %>%
    mutate(cond = factor(cond, levels=thfs$cond.s)) %>%
    arrange(i, cond) %>%
    mutate(npos = map_int(cnts, length)) %>%
    group_by(i, npos, cidxs) %>% nest() %>% rename(cnts=data) %>% ungroup() %>%
    select(i,npos,cidxs,cnts) %>%
    rename(tidx = i) %>% mutate(tidx=glue("t{tidx}"))
#}}}
tss2 = tss %>% inner_join(tic, by=c('tidx')) %>%
    select(tidx,tpm,IQR,support,tid,peakType,gid,npos,everything())

#{{{ gTSS - merge replicates
ti = findOverlaps(SE_ctss_m, rowRanges(SE_gtss)) %>% as_tibble() %>%
    select(cidx=1, i=2)
tic = assay(SE_ctss_m, 'counts') %>% as_tibble() %>%
    mutate(cidx = 1:n()) %>%
    inner_join(ti, by='cidx') %>%
    gather(cond, cnt, -i, -cidx) %>%
    arrange(i, cond, cidx) %>%
    group_by(i, cond) %>% summarise(cnts=list(cnt), cidxs=list(cidx)) %>%
    ungroup() %>%
    mutate(cond = factor(cond, levels=thfs$cond.s)) %>%
    arrange(i, cond) %>%
    mutate(npos = map_int(cnts, length)) %>%
    group_by(i, npos, cidxs) %>% nest() %>% rename(cnts=data) %>% ungroup() %>%
    select(i,npos,cidxs,cnts) %>%
    rename(gidx = i) %>% mutate(gidx=glue("g{gidx}"))
#}}}
gtss2 = gtss %>% inner_join(tic, by=c('gidx')) %>%
    select(gidx,tpm,gid,n_tss,npos,everything())

r3 = list(tss = tss2, gtss=gtss2, SE_ctss_m=SE_ctss_m, ctss_m=ctss_m)
fo = glue('{dirw}/03.rep.merged.rds')
saveRDS(r3, fo)
#}}}
# run cg.job.03.R


#{{{ export to IGV
to = tss %>%
    mutate(start=start-1, tstart=domi-1, tend=domi, width=end-start,
           score=support,
           id = sprintf("t%05d", i),
           rgb = ifelse(is.na(gid), '0,128,128', '0,0,0'),
           blockCnt=1, blockSizes=sprintf("%d,",width), blockStarts="0,") %>%
    mutate(note1 = sprintf("width=%d", width)) %>%
    mutate(note2 = sprintf("IQR=%d", IQR)) %>%
    mutate(note3 = sprintf("peakType=%s", peakType)) %>%
    select(chrom, start, end, id, score, srd, tstart, tend, rgb,
           blockCnt, blockSizes, blockStarts,
           note1, note2, note3, gid) %>%
    arrange(chrom, start, end)
fo = glue('{dirw}/91.ftss.bed')
write_tsv(to,fo, col_names=F)
system("bgzip -c 91.ftss.bed > 91.ftss.bed.gz")
system("tabix -p bed 91.ftss.bed.gz")

t_enh = rowRanges(enh) %>% as_tibble()
to = t_enh %>%
    mutate(start=start-1, tstart=thick.start-1, tend=thick.end,
           score=ifelse(score>1000, 1000, round(score)),
           id=sprintf("enh%04d", 1:n()),
           rgb = 0,
           blockCnt=1, blockSizes=sprintf("%d,",width), blockStarts="0,") %>%
    mutate(strand=ifelse(strand=='*', '.', strand)) %>%
    mutate(note1 = sprintf("width=%d", width)) %>%
    mutate(note3 = sprintf("peakType=%s", peakTxType)) %>%
    select(chrom=seqnames, start, end, id, score, strand, tstart, tend, rgb,
           blockCnt, blockSizes, blockStarts,
           note1, note3) %>%
    arrange(chrom, start, end)
fo = file.path(dirw, '92.enh.bed')
write_tsv(to,fo, col_names=F)
# bgzip -c 92.enh.bed > 92.enh.bed.gz
# tabix -p bed 92.enh.bed.gz
#}}}

##### OBSOLETE #####
#{{{ [old] tss meta plots (sample-wise)
x0 = rowRanges(SE_tss) %>% as_tibble() %>%
    mutate(tstart=thick.start-1, tend=thick.start, start=start-1) %>%
    select(chrom=1,tstart,tend,start,end,strand,width,tpm=score,support)
x1 = assay(ftss, 'TPM') %>% as_tibble()
x = x0 %>% bind_cols(x1)
fo = file.path(dirw, 'tss.bed')
write_tsv(x, fo, col_names=F)
system("intersectBed -wao -s -a tss.bed -b $genome/data2/Zmays_B73/35.metaplot.bin.bed > x.bed")

fi = file.path(dirw, 'x.bed')
ti = read_tsv(fi,col_names=F)
colnames(ti) = c('chrom','tstart','pos','start','end','srd','size','tpm','support',
           colnames(ftss), 'chrom2','start2','end2','gid','i','srd2','bp')

ti2 = ti %>% select(gid,i,starts_with('s')) %>%
    filter(srd != '.') %>%
    select(-start,-size,-srd,-support,-srd2) %>%
    mutate(i = as.integer(i)) %>% filter(i != -1) %>%
    gather(sid, tpm, -i, -gid) %>%
    group_by(sid, i) %>% summarise(tpm = sum(tpm)) %>% ungroup()

#{{{ meta plot for picked samples
tp = ti2 %>%
    inner_join(th, by=c('sid'='SampleID')) %>%
    filter(Treatment %in% c("normal",'cold_control','drought_control')) %>%
    mutate(lab = sprintf("%s (%s)", lab, sid))
labs = tp %>% distinct(Tissue,Genotype,Treatment,lab) %>% arrange(Tissue, Genotype, Treatment) %>% pull(lab)
tp = tp %>% mutate(lab = factor(lab, levels=labs))

tps = tibble(idx=c(51,150), lab=c('TSS','TES'))
p = ggplot(tp) +
    geom_line(aes(i, tpm)) +
    #geom_vline(xintercept=tps$idx, linetype='dashed', size=.5) +
    scale_x_continuous(breaks=tps$idx, labels=tps$lab, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(expand=expansion(mult=c(.05,.05))) +
    facet_wrap(~lab, ncol=4) +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T, strip.style='light',
           xgrid=T,xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T)
fo = file.path(dirw, '13.metaplot.pdf')
ggsave(p, file=fo, width=12, height=10)
#}}}
#}}}

#{{{ DE TSS calling
require(DESeq2)
fi = file.path(dirw, '03.cnts.rds')
x = readRDS(fi)
th1 = thf %>% rename(cond=cond.s)
th2 = thf %>% rename(cond=Genotype)
th3 = thf %>% rename(cond=Tissue)
tg = x %>% select(i, meta) %>% unnest(meta) %>% rename(gid0 = gid) %>%
    mutate(gid = sprintf("c%05d", i))
tm = assay(ftss, 'counts') %>% as_tibble() %>% mutate(gid=tg$gid) %>%
    select(gid, everything()) %>%
    gather(SampleID, ReadCount, -gid)

get_ds <- function(condR, cond, dds, gids) {
    #{{{
    res1 = results(dds, contrast = c("cond",condR,cond), pAdjustMethod='fdr')
    stopifnot(rownames(res1) == gids)
    tibble(gid = gids, padj = res1$padj, log2fc = res1$log2FoldChange) %>%
        replace_na(list(padj = 1))
    #}}}
}
run_deseq2 <- function(th, tm, ct) {
    #{{{
    #cat(sprintf('--> working on %s - %s\n', gene_alias, Tissue))
    th1 = th
    tm1 = tm %>% filter(SampleID %in% th1$SampleID)
    #{{{ prepare data
    vh = th1 %>% arrange(SampleID)
    vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
    gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 5)) %>%
        filter(n.sam > .1 * nrow(vh)) %>% pull(gid)
    vm = tm1 %>% filter(gid %in% gids) %>%
        select(SampleID, gid, ReadCount)
    #x = readcount_norm(vm)
    #mean.lib.size = mean(x$tl$libSize)
    #vm = x$tm
    vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
    vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
    stopifnot(identical(rownames(vh.d), colnames(vm.d)))
    #}}}
    # DESeq2
    dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design=~cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds, fitType = 'parametric')
    disp = dispersions(dds)
    #dds = nbinomLRT(dds, reduced = ~ 1)
    dds = nbinomWaldTest(dds)
    resultsNames(dds)
    res = ct %>% mutate(ds = map2(condR, cond, get_ds, dds = dds, gids = gids))
    res
    #}}}
}

ct1 = tibble(cond  = c("Ar", 'As', 'Ar', 'Br'),
             condR = c("Br", 'Bs', 'As', 'Bs'))
ds1 = run_deseq2(th1, tm, ct1)
ct2 = tibble(cond = 'A632', condR = 'B73')
ds2 = run_deseq2(th2, tm, ct2)
ct3 = tibble(cond = 'root', condR = 'shoot')
ds3 = run_deseq2(th3, tm, ct3)
ds = rbind(ds1,ds2,ds3)

fo = file.path(dirw, '20.de.tss.rds')
saveRDS(ds, fo)
#}}}

#{{{ DE TSS analysis
fi = file.path(dirw, '20.de.tss.rds')
ds = readRDS(fi)

#{{{ volcano plot
tp = ds %>% unnest(ds) %>% mutate(logp = -log10(padj)) %>%
    mutate(logp = ifelse(logp > 30, 30, logp)) %>%
    mutate(log2fc = ifelse(log2fc > 10, 10, log2fc)) %>%
    mutate(log2fc = ifelse(log2fc < -10, -10, log2fc)) %>%
    mutate(sig = ifelse(padj < .05 & abs(log2fc) >= 1, 'y', 'n')) %>%
    mutate(pan = str_c(cond, condR, sep=' vs '))
tps = tp %>% count(cond, condR, pan, sig) %>% filter(sig=='y') %>%
    mutate(txt = sprintf("%d dTSSs", n))
p = ggplot(tp) +
    geom_point(aes(x=log2fc, y=logp, col=sig), size=1, alpha=.5) +
    geom_vline(xintercept=c(-1,1), linetype='dashed', color='gray') +
    geom_hline(yintercept=-log10(.05), linetype='dashed', color='gray') +
    geom_text(tps, mapping=aes(label=txt),x=0,y=30, color='black', size=3) +
    scale_x_continuous(name='log2FoldChange',expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(name='-log10(pval.adj)', expand=expansion(mult=c(.05,.05))) +
    scale_color_manual(values=pal_aaas()(10)[c(9,2)]) +
    facet_wrap(~pan, ncol=2) +
    otheme(legend.pos='none', legend.dir='v', legend.title=T, strip.style='light',
           xtick=T,ytick=T, xtitle=T, ytitle=T, xtext=T, ytext=T)
fo = file.path(dirw, '21.volcano.pdf')
ggsave(p, file=fo, width=8, height=8)
#}}}

x1 = ds %>% unnest(ds) %>%
    mutate(sig = ifelse(padj < .05 & abs(log2fc) >= 1, 'y', 'n')) %>%
    mutate(pan = str_c(cond, condR, sep=' vs ')) %>%
    rename(pid = gid) %>%
    filter(sig == 'y') %>% select(-sig) %>%
    #inner_join(tg, by='pid') %>%
    mutate(drc = ifelse(log2fc > 0, "+", "-"))
x2 = x1 %>%
    group_by(cond,condR,pan,gid) %>%
    summarise(ndts=n(), drc=str_c(sort(unique(drc)),collapse='')) %>%
    ungroup() %>%
    mutate(ndts = ifelse(ndts > 1, "2+", ndts))
x2 %>%
    count(pan, ndts, drc) %>%
    spread(drc, n) %>% print(n=50)

x = dtu %>% unnest(ds) %>% mutate(dtu = padj<.05 & abs(log2fc)>=1)
x %>% filter(dtu) %>% count(cond, condR, dtu)
x %>% filter(padj < .05, abs(log2fc) >= 1) %>%
    separate(cond, c('tis','gt'), sep='_', remove=F) %>%
    separate(condR, c('tisR','gtR'), sep='_', remove=F) %>%
    filter(tis==tisR) %>% count(gid)  %>% count(n)
x %>% filter(padj < .05, abs(log2fc) >= 1) %>%
    separate(cond, c('tis','gt'), sep='_', remove=F) %>%
    separate(condR, c('tisR','gtR'), sep='_', remove=F) %>%
    filter(gt==gtR) %>% count(gid) %>% count(n)
#}}}

#{{{ # enhancer candidates
bc = clusterBidirectionally(fctss, balanceThreshold=0.95) %>%
    calcBidirectionality(samples=fctss)
table(bc$bidirectionality)
fbc = bc %>% subset(bidirectionality>2)
#subsetBySupport(fbc, inputAssay="counts", unexpressed=0, minSamples=2)

enh = quantifyClusters(fctss, clusters=fbc, inputAssay='counts') %>%
    assignTxType(txModels=txdb, outputColumn='txType') %>%
    assignTxType(txModels=txdb, outputColumn='peakTxType', swap='thick')
fenh = enh %>% subset(peakTxType %in% c("intron","intergenic"))
#}}}

#{{{ [old] check TSS ovlp w. TEs
fo = glue("{dirw}/t1.bed")
to = tss %>% select(chrom,start,end,i) %>% mutate(start=start-1)
write_tsv(to, fo, col_names=F)

system(glue('intersectBed -u -f 0.5 -a t1.bed -b ~/projects/genome/data2/TE/11.B73.bed > t2.bed'))
ti2 = read_tsv('t2.bed', col_names=c('chrom','start','end','i')) %>%
    select(i) %>% mutate(te = 'TE')

ttypes = gcfg$gene %>% distinct(ttype) %>% arrange(ttype) %>% pull(ttype)
ttypes = c(ttypes[c(3,1,2,4:12)], 'intergenic','antisense')
ptypes = c(levels(tss$peakType), 'TE', 'non-TE')
tss2 = tss %>% left_join(gcfg$gene %>% select(gid,ttype), by='gid') %>%
    left_join(ti2, by=c('i')) %>%
    replace_na(list(te='non-TE')) %>%
    mutate(ttype = ifelse(peakType %in% c('intergenic','antisense'), as.character(peakType), ttype)) %>%
    mutate(ptype = ifelse(peakType %in% c('intergenic','antisense'), te, as.character(peakType))) %>%
    mutate(ttype = factor(ttype, levels=ttypes)) %>%
    mutate(ptype = factor(ptype, levels=ptypes))
tss2 %>% count(ttype, ptype) %>% spread(ttype, n) %>% print(n=20, width=Inf)

tp = tss2 %>%
    rename(tag1=ttype,tag2=ptype) %>% count(tag1, tag2)
p = cmp_proportion1(tp, xangle=15, oneline=F, legend.title='', expand.x=c(.05,.05),
    margin.l=2, fills=pal_simpsons()(11)) +
    o_margin(1.5,.3,.01,1.5)
#
fo = glue("{dirw}/08.te.prop.pdf")
ggsave(p, file=fo, width=5, height=5)
p = cmp_cnt1(tp, xangle=15, oneline=F, legend.title='', expand.x=c(.05,.05),
    margin.l=2, fills=pal_simpsons()(11)) +
    o_margin(1.5,.3,.01,1.5)
#
fo = glue("{dirw}/08.te.cnt.pdf")
ggsave(p, file=fo, width=5, height=5)
#fo = glue("{dirf}/sf02b.rds")
#saveRDS(p, fo)
#}}}

