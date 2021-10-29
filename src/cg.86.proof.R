source("functions.R")
dirw = glue("{dird}/86_proof")

#{{{ extract first 150bp sequences
ti = gcfg$gene.loc %>%
    filter(etype=='exon') %>%
    group_by(gid,tid,chrom,srd) %>%
    summarise(start = min(start), end = max(end)) %>% ungroup() %>%
    mutate(tss = ifelse(srd=='-', end, start))
to = ti %>%
    mutate(rb = ifelse(srd=='-', tss-150+1, tss)) %>%
    mutate(re = ifelse(srd=='-', tss, tss+150-1)) %>%
    mutate(rb = rb - 1, score = '.') %>%
    select(chrom, rb, re, gid, score, srd)
fo = glue("{dirw}/01.bed")
write_tsv(to, fo, col_names=F)
# bedtools getfasta -fi $ref/10.fasta -bed 01.bed -fo 01.fa -nameOnly
#}}}

aln = 'hisat2'
aln = 'star'
#{{{ mapping stats
fi = glue("{dirw}/11_{aln}/bamstats.tsv")
tm = read_tsv(fi)
#
ti2 = tm %>% select(sid,unpair,unpair_map,unpair_map0, unpair_map_hq, unpair_map_hq0, unpair_unmap) %>%
    separate(sid, c('sid','gt'), sep='-') %>%
    mutate(gt = str_replace(gt, 'Zmays_', '')) %>%
    mutate(gt = str_replace(gt, '.stat', ''))
ti2 %>% mutate(x=unpair-unpair_map-unpair_unmap) %>% count(x)
types=c('unmap','multi_mm','uniq_mm','multi_pm','uniq_pm')
tp = ti2 %>% mutate(multi=unpair_map-unpair_map_hq) %>%
    mutate(multi_pm=unpair_map0-unpair_map_hq0, multi_mm=multi-multi_pm) %>%
    mutate(uniq_pm=unpair_map_hq0, uniq_mm=unpair_map_hq-uniq_pm) %>%
    select(sid, unmap=unpair_unmap, multi_pm,multi_mm,uniq_pm,uniq_mm) %>%
    inner_join(ths %>% select(sid,pnl,grp,grp.x,grp.xlab), by=c('sid')) %>%
    select(-sid) %>%
    gather(type, cnt, -pnl, -grp, -grp.x, -grp.xlab) %>% rename(x=pnl) %>%
    mutate(type=factor(type,levels=types)) %>%
    select(pnl=grp, tag1=grp.x, tag2=type, n=cnt) %>% mutate(n=n/1e6)

cols5 = pal_simpsons()(8)[c(3,1,2,4,5)]
p1 = cmp_proportion(tp, xangle=0, alph=.4, acc=.1,
    lab.size=2, oneline=F, nc=4, fills=cols5,
    expand.x=c(.02,.02), expand.y=c(.01,.04), legend.title='',
    legend.pos='top.center.out', legend.dir='h', legend.vjust=-1,
    margin= c(.5,.1,.1,.1)) 
fo = glue("{dirw}/12.map.rate.{aln}.pdf")
ggsave(p1, filename = fo, width = 10, height = 6)
#}}}




