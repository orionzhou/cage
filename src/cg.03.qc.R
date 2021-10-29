require(CAGEfightR)
require(GenomicFeatures)
source("functions.R")
dirw = glue("{dird}/03_qc")
#setwd(dirw)

org='Zmays_B73v5'
#{{{ call TSS & gTSS
txdb = load_txdb(org)
dsg = thf %>% mutate(Name=sid)
diri = glue("~/projects/s3/zhoup-nfo/zm.cg20a.2/34_ctss")

SE_ctss = read_cage_bws(diri, dsg, minSupport=2, genome=org)
seqlevels(SE_ctss) = seqlevels(txdb)
rowRanges(SE_ctss) = sort(rowRanges(SE_ctss))
#SE_rtss = make_tss(SE_ctss, unexpressed=.5, minSamples=2)
SE_rtss = make_tss(SE_ctss, txdb, unexpressed=1, minSamples=3)

#{{{ all TSS clusters
#characterize condition-wise TPM
x = thf %>%
    group_by(grp) %>% summarise(sids = list(sid)) %>% ungroup() %>%
    mutate(x = map(sids, subset_samples, ex=SE_rtss)) %>%
    select(-sids) %>% unnest(x) %>%
    group_by(i) %>% nest() %>% rename(tpm.grp=data)
rtss = rowRanges(SE_rtss) %>% as_tibble() %>%
    mutate(i = 1:n()) %>% mutate(tpm = score/ncol(assay(SE_rtss))) %>%
    inner_join(x, by='i') %>%
    mutate(tidx = glue("t{i}")) %>%
    select(tidx, chrom=seqnames,start,end,srd=strand,tpm,domi=thick.start,
           IQR,entropy, support, tid=txID, peakType=peakTxType, gid=geneID, tpm.grp)
rtss %>% count(peakType)
rtss %>% select(tpm.grp) %>% unnest(tpm.grp) %>% group_by(grp) %>% summarise(n_expressed=sum(tpm>=1))
rtss %>% select(tpm.grp) %>% unnest(tpm.grp) %>% group_by(grp) %>% summarise(tpm=sum(tpm))
#}}}

SE_tss = SE_rtss %>% subset(peakTxType %in% peakTypes) #%>% calcTotalTags() %>% calcTPM() %>% calcPooled()
#{{{ filter to promoter TSS clusters 
# characterize condition-wise TPM
#x = thf %>%
    #group_by(grp) %>% summarise(sids = list(sid)) %>% ungroup() %>%
    #mutate(x = map(sids, subset_samples, ex=SE_tss)) %>%
    #select(-sids) %>% unnest(x) %>%
    #group_by(i) %>% nest() %>% rename(tpm.grp=data)
#tss = rowRanges(SE_tss) %>% as_tibble() %>%
    #mutate(i = 1:n()) %>% mutate(tpm = score/ncol(assay(SE_tss))) %>%
    #inner_join(x, by='i') %>%
    #mutate(tidx = glue("t{i}")) %>%
    #select(tidx, chrom=seqnames,start,end,srd=strand,tpm,domi=thick.start,
           #IQR,entropy, support, tid=txID, peakType=peakTxType, gid=geneID, tpm.grp)
tss = rtss %>% filter(peakType %in% types)
tss %>% count(peakType)
tss %>% select(tpm.grp) %>% unnest(tpm.grp) %>% group_by(grp) %>% summarise(n_expressed=sum(tpm>=1))
tss %>% select(tpm.grp) %>% unnest(tpm.grp) %>% group_by(grp) %>% summarise(tpm=sum(tpm))
#}}}

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
idx = idxs[1]
x1[(idx-1):idx,]
gr = gr[-idxs]
isDisjoint(gr)
#}}}
seqinfo(gr) = seqinfo(SE_ctss)
#}}}

SE_gtss = SE_ctss %>% quantifyClusters(clusters=gr, inputAssay='counts')
#colData(SE_gtss)[,'totalTags'] = colData(SE_tss)$totalTags
SE_gtss = SE_gtss %>%
    calcTotalTags() %>% calcTPM() %>% calcPooled() %>%
    # assignGeneID(geneModels=txdb, outputColumn='geneID')
    calcShape(pooled=SE_ctss, outputColumn='IQR', shapeFunction=shapeIQR,
        lower=.1, upper=.9) %>%
    calcShape(pooled=SE_ctss, shapeFunction=shapeEntropy, outputColumn='entropy')

x = thf %>%
    group_by(grp) %>% summarise(sids = list(sid)) %>% ungroup() %>%
    mutate(x = map(sids, subset_samples, ex=SE_gtss)) %>%
    select(-sids) %>% unnest(x) %>%
    group_by(i) %>% nest() %>% rename(tpm.grp=data)
gtss = rowRanges(SE_gtss) %>% as_tibble() %>%
    mutate(i = 1:n()) %>% mutate(tpm = score/ncol(assay(SE_gtss))) %>%
    select(i, chrom=seqnames,start,end,srd=strand,tpm,IQR, entropy) %>%
    inner_join(x1, by=c('chrom','start','end','srd')) %>%
    inner_join(x, by='i') %>%
    rename(gidx=i) %>% mutate(gidx = glue("g{gidx}"))
gtss %>% count(n_tss)
gtss %>% count(n_tss > 1)
gtss %>% select(tpm.grp) %>% unnest(tpm.grp) %>% group_by(grp) %>% summarise(n_expressed=sum(tpm>=1))
gtss %>% select(tpm.grp) %>% unnest(tpm.grp) %>% group_by(grp) %>% summarise(tpm=sum(tpm))
#}}}

r1 = list(SE_ctss=SE_ctss, SE_rtss=SE_rtss, SE_tss=SE_tss, SE_gtss=SE_gtss,
          rtss=rtss, tss=tss, gtss=gtss)
fo = glue("{dirw}/01.tss.gtss.rds")
saveRDS(r1, fo)
#}}}

#{{{ subset to B73 TSS & gTSS
fi = glue("{dirw}/01.tss.gtss.rds")
r1 = readRDS(fi)
SE_ctss=r1$SE_ctss; SE_rtss=r1$SE_rtss; SE_tss=r1$SE_tss; SE_gtss=r1$SE_gtss;
rtss=r1$rtss; tss=r1$tss; gtss=r1$gtss
#
tcond = thf %>% distinct(grp,gt,tis)

ti1 = rtss %>% select(i=tidx,tpm.grp) %>% unnest(tpm.grp) %>%
    mutate(s = ifelse(support/nrep >= .2, "E", "x")) %>%
    inner_join(tcond, by='grp') %>%
    select(i,grp,gt, tis, s) %>% mutate(tag = 'tss')
ti2 = gtss %>% select(i=gidx,tpm.grp) %>% unnest(tpm.grp) %>%
    mutate(s = ifelse(support/nrep >= .2, "E", "x")) %>%
    inner_join(tcond, by='grp') %>%
    select(i,grp,gt, tis, s) %>% mutate(tag = 'gene')
t_exp = ti1 %>% rbind(ti2) %>% mutate(tag = factor(tag, levels=c('tss','gene')))

tb = t_exp %>% filter(gt == 'B73', s=='E') %>% distinct(i, tag)
tb %>% count(tag)
tidxs = tb %>% filter(tag == 'tss') %>% pull(i)
gidxs = tb %>% filter(tag == 'gene') %>% pull(i)
tss %>% filter(tidx %in% tidxs)
tss %>% filter(tidx %in% tidxs) %>% distinct(gid)

r2 = list(SE_ctss=SE_ctss, SE_rtss=SE_rtss, SE_tss=SE_tss, SE_gtss=SE_gtss,
          rtss=rtss, tss=tss, gtss=gtss,
          tidxs=tidxs, gidxs=gidxs)
fo = glue("{dirw}/02.tss.gtss.B73.rds")
saveRDS(r2, fo)

# B73 stats
tss.b = tss %>% filter(tidx %in% tidxs)
gtss.b = gtss %>% filter(gidx %in% gidxs)
nrow(gtss.b)
gtss.b %>% count(n_tss == 1)
tss.b %>% count(peakType== 'intergenic')
tss.b %>% count(peakType)
gtss.b %>% inner_join(gcfg$gene %>% select(gid,ttype), by='gid') %>% count(ttype)

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

#{{{ merge replicates & make counts table for TSSs & gTSSs
fi = glue("{dirw}/02.tss.gtss.B73.rds")
r2 = readRDS(fi)
SE_ctss=r2$SE_ctss; SE_rtss=r2$SE_rtss; SE_tss=r2$SE_tss; SE_gtss=r2$SE_gtss;
rtss=r2$rtss; tss=r2$tss; gtss=r2$gtss
tcond = thf %>% distinct(grp,gt,tis)

#{{{ make replicate-merged ctss object
rowsum1 <- function(sids, ctss) rowSums(assay(ctss[,sids],'counts'))
x = thf %>% group_by(grp) %>%
    summarise(sids=list(sid)) %>% ungroup() %>%
    mutate(cnts = map(sids, rowsum1, ctss=SE_ctss))
ns = nrow(tcond)
xm = matrix(x %>% select(cnts) %>% unnest(cnts) %>% pull(cnts), ncol=ns, byrow=F)
xh = tcond %>% rename(Name=grp)
colnames(xm) = xh$Name
#
SE_ctss_m = SummarizedExperiment(assays=SimpleList(counts=xm),
                         rowRanges=rowRanges(SE_ctss)[,-c(1,2)],
                         colData=xh) %>%
    calcTotalTags(inputAssay="counts", outputColumn="totalTags")
ctss_m = rowRanges(SE_ctss)[,-c(1,2)] %>% as.data.frame() %>% as_tibble() %>%
    mutate(cidx=1:n()) %>% select(cidx,chrom=1,pos,strand)
#}}}

#{{{ rTSS - merge replicates
ti = findOverlaps(SE_ctss_m, rowRanges(SE_rtss)) %>% as_tibble() %>%
    select(cidx=1, i=2)
tic = assay(SE_ctss_m, 'counts') %>% as_tibble() %>%
    mutate(cidx = 1:n()) %>%
    inner_join(ti, by='cidx') %>%
    gather(grp, cnt, -i, -cidx) %>%
    arrange(i, grp, cidx) %>%
    group_by(i, grp) %>% summarise(cnts=list(cnt), cidxs=list(cidx)) %>%
    ungroup() %>%
    mutate(grp = factor(grp, levels=tcond$grp)) %>%
    arrange(i, grp) %>%
    mutate(npos = map_int(cnts, length)) %>%
    group_by(i, npos, cidxs) %>% nest() %>% rename(cnts=data) %>% ungroup() %>%
    select(i,npos,cidxs,cnts) %>%
    rename(tidx = i) %>% mutate(tidx=glue("t{tidx}"))
#}}}
rtss2 = rtss %>% inner_join(tic, by=c('tidx')) %>%
    select(tidx,tpm,IQR,entropy,support,tid,peakType,gid,npos,everything())

tss2 = rtss2 %>% filter(peakType %in% peakTypes)

#{{{ gTSS - merge replicates
ti = findOverlaps(SE_ctss_m, rowRanges(SE_gtss)) %>% as_tibble() %>%
    select(cidx=1, i=2)
tic = assay(SE_ctss_m, 'counts') %>% as_tibble() %>%
    mutate(cidx = 1:n()) %>%
    inner_join(ti, by='cidx') %>%
    gather(grp, cnt, -i, -cidx) %>%
    arrange(i, grp, cidx) %>%
    group_by(i, grp) %>% summarise(cnts=list(cnt), cidxs=list(cidx)) %>%
    ungroup() %>%
    mutate(grp = factor(grp, levels=tcond$grp)) %>%
    arrange(i, grp) %>%
    mutate(npos = map_int(cnts, length)) %>%
    group_by(i, npos, cidxs) %>% nest() %>% rename(cnts=data) %>% ungroup() %>%
    select(i,npos,cidxs,cnts) %>%
    rename(gidx = i) %>% mutate(gidx=glue("g{gidx}"))
#}}}
gtss2 = gtss %>% inner_join(tic, by=c('gidx')) %>%
    select(gidx,tpm,IQR,entropy,gid,n_tss,npos,everything())

r3 = list(rtss=rtss2, tss=tss2, gtss=gtss2, SE_ctss_m=SE_ctss_m, ctss_m=ctss_m,
    tidxs=r2$tidxs, gidxs=r2$gidxs)
fo = glue('{dirw}/03.rep.merged.rds')
saveRDS(r3, fo)
#}}}

# characterize shifting scores by running cg.job.03.R

#{{{ shape characterization
fi = glue('{dirw}/03.rep.merged.rds')
r3 = readRDS(fi)
rtss=r3$rtss; tss=r3$tss; gtss=r3$gtss

#{{{ plot shape IQR vs. entropy 
calc_SI <- function(t_cnts, nsam=6) {
    #{{{ calc shape index using tag count matrix
    y = colSums(matrix(unlist(t_cnts %>% pull(cnts)), byrow=T, nrow=nsam))
    ty = tibble(i=1:length(y), cnt=y) %>%
        mutate(p = cnt/sum(y)) %>% mutate(ent = p*(log2(p)))
    si = 2 + sum(ty$ent)
    si
    #}}}
}

tp1 = rtss %>% select(x=entropy, y=IQR) %>% mutate(pnl = 'All TSS clusters')
tp2 = tss %>% select(x=entropy, y=IQR) %>% mutate(pnl = 'Genic TSS cluters')
tp3 = gtss %>% select(x=entropy, y=IQR) %>% mutate(pnl = 'Merged genic TSS cluters')
#{{{ plot
tp = rbind(tp1,tp2,tp3) %>% mutate(pnl=as_factor(pnl)) %>%
    mutate(x = pmin(7, x)) %>%
    mutate(y = pmin(300,y))
tpl = tp %>% group_by(pnl) %>% nest() %>% ungroup() %>%
    mutate(n = map_int(data, nrow)) %>%
    mutate(fit = map(data, ~ lm(y~x, data=.x))) %>%
    mutate(tidied = map(fit, glance))%>%
    unnest(tidied) %>%
    mutate(lab = glue("N={number(n,accuracy=1,big.mark=',')}<br>adjusted R<sup>2</sup> = {number(adj.r.squared,accuracy=.01)}"))
xlab = "Shape Entropy (SE)"
ylab = "Interquantile width (IQR)"
xbrks = seq(0,7, by=2)
ybrks = c(0,100,200,300)
#
p = ggplot(tp, aes(x=x,y=y)) +
    #geom_hline(yintercept=c(10,50), linetype='dashed', size=.5) +
    #geom_vline(xintercept=-1, linetype='dashed', size=.5) +
    geom_hex(bins=80) +
    geom_richtext(data=tpl, aes(x=0,y=300,label=lab,hjust=0,vjust=1), size=2.5) +
    scale_x_continuous(name=xlab,breaks=xbrks,expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name=ylab,breaks=ybrks,expand=expansion(mult=c(.02,.02))) +
    scale_fill_viridis(name='density', direction=-1) +
    facet_wrap(~pnl, nrow=1) +
    otheme(legend.pos='bottom.right', legend.title=T, legend.dir='v',
           legend.box='h', legend.vjust=.5, panel.spacing=.3,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) + o_margin(.3,.3,.3,.3)
#}}}
#}}}
fo = glue("{dirw}/10.IQR.entropy.pdf")
ggsave(p, file=fo, width=10, height=4)

#{{{ make shape categories
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
    mutate(shape.si = ifelse(SI <= -1, "SI <= -1", "SI > -1")) %>%
    mutate(shape.si = factor(shape.si)) %>%
    mutate(shape = factor(shape, levels=shapes)) %>%
    mutate(shape.s = factor(shape.s, levels=shapess))
gtss = gtss %>% rename(iqr=IQR) %>%
    mutate(shape = map_chr(iqr, iqr2shape)) %>%
    mutate(shape.s = map_chr(iqr, iqr2shape, opt='s')) %>%
    mutate(shape.si = ifelse(SI <= -1, "SI <= -1", "SI > -1")) %>%
    mutate(shape.si = factor(shape.si)) %>%
    mutate(shape = factor(shape, levels=shapes)) %>%
    mutate(shape.s = factor(shape.s, levels=shapess))
tss %>% count(shape) %>% mutate(p=n/sum(n))
gtss %>% count(shape) %>% mutate(p=n/sum(n))

tp = tss %>% count(shape,shape.si) %>%
    rename(tag1 = shape.si, tag2 = shape)
p = cmp_cnt1(tp, ytext=T, ypos='right', legend.title='Shape Index') +
    o_margin(.1,.3,.1,.3) +
    theme(legend.position='none')
fo = glue("{dirw}/10.si.iqr.prop.pdf")
ggsave(p, file=fo, width=3, height=5)
#}}}

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
fo = glue('{dirw}/11.ftype.panel.pdf')
ggsave(p1, file=fo, width=7, height=7)

#{{{ Allison and Maike's TSS shape
tp = tss %>% select(shape=shape.s, tpm.cond) %>%
    unnest(tpm.cond) %>%
    mutate(s = ifelse(support/nrep >= .2, "E", "x")) %>%
    filter(s == 'E') %>%
    count(cond, shape) %>%
    select(tag1=cond,tag2=shape,n)

p = cmp_proportion1(tp, ytext=T, ypos='right', legend.title='shape:') +
    o_margin(.1,.3,.1,.3) +
    theme(legend.position='none')
fo = glue("{dirw}/13.shape.batch.pdf")
ggsave(p, file=fo, width=5, height=5)
#}}}

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
    #filter(peakType != 'intergenic') %>%
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
fp = glue('{dirf}/sf03b.pdf')
fp = glue('{dirw}/10.tsne.pdf')
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
ti %>% filter(gt=='B73',tis=='shoot',s=='E') %>% distinct(i,tag) %>% count(tag)
ti %>% filter(gt=='B73',s=='E') %>% distinct(i,tag) %>% count(tag)
ti %>% filter(gt == 'B73', s=='E') %>% distinct(i,tag) %>% count(tag)
tib = ti %>% filter(gt == 'B73', s=='E', tag=='gene')
tib %>% count(cond) %>% pull(n) %>% mean()
tib %>% distinct(i)

tib = ti %>% filter(gt == 'B73', s=='E')
tia = ti %>% filter(gt == 'A632', s=='E') %>% distinct(i, tag) %>%
    filter(! i %in% tib$i) %>% count(tag) %>% print(n=4)

# alt-TSSs in non-shoot B73 tissues
tig = ti %>% filter(gt=='B73',tis=='shoot',s=='E',tag=='gene')
tit = ti %>% filter(gt=='B73',tag=='tss') %>%
    select(i,tis,s) %>% spread(tis, s)
tit %>%
    filter(shoot=='x', root=='E'|stem=='E'|husk=='E') %>%
    inner_join(g2t, by=c('i'='tidx')) %>%
    count(gidx %in% tig$i)
#}}}
#}}}

#{{{ export & share - all
fi = glue("{dirw}/01.tss.gtss.rds")
r1 = readRDS(fi)
tss=r1$tss; gtss=r1$gtss

#{{{ igv
diro = glue("{dird}/09_igv")
to = tss %>%
    mutate(start=start-1, tstart=domi-1, tend=domi, width=end-start,
           score=support, id = tidx,
           rgb = ifelse(is.na(gid), '0,128,128', '0,0,0'),
           blockCnt=1, blockSizes=sprintf("%d,",width), blockStarts="0,") %>%
    mutate(note1 = sprintf("width=%d", width)) %>%
    mutate(note2 = sprintf("IQR=%d", IQR)) %>%
    mutate(note3 = sprintf("peakType=%s", peakType)) %>%
    select(chrom, start, end, id, score, srd, tstart, tend, rgb,
           blockCnt, blockSizes, blockStarts,
           note1, note2, note3, gid) %>%
    arrange(chrom, start, end)

fo = glue('{dirw}/91.tss.bed')
write_tsv(to,fo, col_names=F)
system("bgzip -c 91.tss.bed > 91.tss.bed.gz")
system("tabix -p bed 91.tss.bed.gz")
#}}}

#{{{ TSS & gTSS
diro = glue("{dird}/91_share")

to = tss %>% select(-tpm.cond)
fo = glue('{diro}/01.tss.tsv')
write_tsv(to, fo, na='')

to = gtss %>% select(-tidx,-tpm.cond)
fo = glue('{diro}/01.gtss.tsv')
write_tsv(to, fo, na='')
#}}}

#{{{ B73
fi = glue("{dirw}/02.tss.gtss.cond.rds")
r2 = readRDS(fi)
tss=r2$tss.b; gtss=r2$gtss.b

# B73 TSSs
to = tss %>%
    mutate(start=start-1, tstart=domi-1, tend=domi, width=end-start,
           score=support, id = tidx,
           rgb = ifelse(is.na(gid), '0,128,128', '0,0,0'),
           blockCnt=1, blockSizes=sprintf("%d,",width), blockStarts="0,") %>%
    mutate(note1 = sprintf("width=%d", width)) %>%
    mutate(note2 = sprintf("IQR=%d", IQR)) %>%
    mutate(note3 = sprintf("peakType=%s", peakType)) %>%
    select(chrom, start, end, id, score, srd, tstart, tend, rgb,
           blockCnt, blockSizes, blockStarts,
           note1, note2, note3, gid) %>%
    arrange(chrom, start, end)

fo = glue('{dirw}/92.tss.B73.bed')
write_tsv(to, fo, col_names=F)
#}}}

#{{{ intergenic single
fi = glue("{dirw}/01.tss.gtss.rds")
r1 = readRDS(fi)
tss=r1$tss
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

to = tss %>% filter(peakType=='intergenic', shape.s=='single') %>%
    mutate(start=domi-5, end=domi+5) %>%
    select(chrom,start,end,tidx) %>% arrange(chrom,start)
fo = glue("{diro}/04.te.single.bed")
write_tsv(to, fo, col_names=F)

fh = "~/projects/s3/zhoup-igv-data/Zmays-B73/chromatin/meta.xlsx"
th = read_xlsx(fh)
fi = glue("{diro}/x.tab")
ti = read_tsv(fi, skip=3,col_names=F)

nsam = 30
ti2 = ti %>% mutate(i=1:n()) %>%
    gather(cidx, cnt, -i) %>% mutate(cidx=as.integer(str_sub(cidx,2)))
ti2s = tibble(sam = sprintf('s%02d', 1:nsam)) %>% crossing(j=1:101) %>%
    arrange(sam, j) %>% mutate(cidx=1:n())
ti3 = ti2 %>% inner_join(ti2s, by='cidx') %>%
    replace_na(list(cnt=0)) %>%
    group_by(sam, j) %>%
    summarise(cnt = sum(cnt)) %>% ungroup()

#{{{ meta plot for picked samples
tp = th %>% mutate(pnl.x=histone, pnl.y=glue("{tissue}_{rep}")) %>%
    mutate(pnl = glue("{histone}_{tissue}_{rep}")) %>%
    inner_join(ti3, by=c('sid'='sam'))
tps = tibble(idx=c(1,51,101), lab=c('-500','TSS','+500'))
p = ggplot(tp) +
    geom_line(aes(j, cnt)) +
    #geom_vline(xintercept=tps$idx, linetype='dashed', size=.5) +
    scale_x_continuous(breaks=tps$idx, labels=tps$lab, expand=expansion(mult=c(.05,.05))) +
    scale_y_continuous(expand=expansion(mult=c(.05,.05))) +
    facet_wrap(pnl~., ncol=5, scale='free') +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T, strip.style='light',
           xgrid=T,xtick=T, ytick=T, xtitle=F, xtext=T, ytext=T)
fo = glue("{diro}/x.pdf")
ggsave(p, file=fo, width=12, height=8)
#}}}
#}}}

#{{{ check ovlp w. AGPv4 annotated TSSs
fi = '~/projects/genome/data/Zmays_B73/50_annotation/10.tsv'
ti = read_tsv(fi)

ti2 = ti %>% group_by(tid) %>%
    summarise(chrom=chrom[1], s=min(start), e= max(end), srd=srd[1]) %>% ungroup() %>%
    mutate(tss = ifelse(srd=='-', e, s))

to2 = ti2 %>% mutate(s = tss-1, e = tss) %>% distinct(chrom,s,e) %>%
    arrange(chrom,s,e)
fo = glue("{dirw}/b73.tss.bed")
write_tsv(to2, fo, col_names=F)
#}}}
#}}}

