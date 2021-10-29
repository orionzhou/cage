source("functions.R")
dirw = glue("{dird}/81_trim")
setwd(dirw)
tha %>% count(tis,gt,cond) %>% print(n=30)

#{{{ CAGE-seq: trimmed read size dist 
merge_cnt <- function(x1,x2,x3,x4, types=c('trimmed','too_short','too_long','untrimmed')) {
    #{{{
    tibble(type=types, x = list(x1,x2,x3,x4)) %>%
        unnest(x) %>%
        mutate(n=as.integer(Count)) %>%
        separate(Length, c('size','size2'), sep='-', fill='right') %>%
        mutate(size = as.numeric(size)) %>% select(type, size, n)
    #}}}
}
ti = th2 %>% #slice(1:2) %>%
    select(sid, pnl) %>%
    mutate(fi = glue("{dirw}/05_trimmed/{sid}.tsv")) %>%
    mutate(x1 = map(fi, read_tsv, col_types='cn')) %>%
    mutate(fi = glue("{dirw}/05_too_short/{sid}.tsv")) %>%
    mutate(x2 = map(fi, read_tsv, col_types='cn')) %>%
    mutate(fi = glue("{dirw}/05_too_long/{sid}.tsv")) %>%
    mutate(x3 = map(fi, read_tsv, col_types='cn')) %>%
    mutate(fi = glue("{dirw}/05_untrimmed/{sid}.tsv")) %>%
    mutate(x4 = map(fi, read_tsv, col_types='cn')) %>%
    mutate(x = pmap(list(x1,x2,x3,x4), merge_cnt)) %>%
    select(sid, pnl, x) %>% unnest(x) %>%
    mutate(type = ifelse(type =='untrimmed', 'too_long', type)) %>%
    group_by(sid, pnl, type, size) %>% summarise(n=sum(n)) %>% ungroup()

size_min = 18; size_max = 25
tis = ti %>% group_by(sid, type) %>%
    summarise(n = sum(n)) %>% ungroup() %>% spread(type, n)
tu = tha %>% select(sid, pnl, spots) %>% inner_join(tis, by=c('sid')) %>%
    mutate(num1 = glue("{number(too_long/1e6,accuracy=.1)}M")) %>%
    mutate(num2 = glue("{number(spots/1e6,accuracy=.1)}M")) %>%
    mutate(pct = glue("{percent(too_long/spots,accuracy=1)}")) %>%
    mutate(num3 = glue("{number(trimmed/1e6,accuracy=.1)}M")) %>%
    mutate(pct3 = glue("{percent(trimmed/spots,accuracy=1)}")) %>%
    mutate(lab = glue("{size_min}-{size_max}bp: {num3}/{num2} ({pct3})"))
tu %>% count(spots-trimmed-too_long-too_short)

xmax = max(ti$size); ymax=max(ti$n)
cols3 = pal_aaas()(4)[c(1,3,2)]
p = ggplot(ti) +
    #geom_segment(data=tps,aes(x=0,xend=xmax,y=0,yend=0), size=.5) +
    geom_histogram(aes(x=size,y=n,fill=type), stat='identity',position='stack') +
    geom_text(data=tu,aes(x=xmax,y=ymax,label=lab), size=2, vjust=1, hjust=1) +
    scale_x_continuous(breaks=c(20,40,60), expand=expansion(mult=c(.03,.03))) +
    scale_y_continuous(expand=expansion(mult=c(0,.01)),
                       labels=unit_format(unit = "M", scale=1e-6)) +
    scale_fill_manual(values=cols3)+
    facet_wrap(~pnl, ncol=9) +
    otheme(legend.pos='top.center.out', legend.vjust=-1, legend.dir='h',
           legend.box='h', panel.border=T,
           xtick=T, xtext=T, ytext=T, ytick=T)
fo = glue("{dirw}/10.pdf")
ggsave(p, file=fo, width=12, height=10)
#}}}

ths = thf %>% arrange(tis,gt,cond,rep) %>%
    mutate(grp = glue("{tis}.{gt}"), grp.xlab = glue("{cond}.{rep}")) %>%
    mutate(grp = as_factor(grp)) %>%
    group_by(grp) %>% mutate(grp.x = 1:n()) %>% ungroup()
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




