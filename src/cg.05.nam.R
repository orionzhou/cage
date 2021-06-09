source("functions.R")
dirw = glue("{dird}/05_nam")
setwd(dirw)
thf2 = thf %>% mutate(cond=glue("{cond.l}.{Replicate}")) %>% select(sid=SampleID,cond)

#{{{ mapping stats
fi = glue("{dirw}/50_final/bamstats.tsv")
ti = read_tsv(fi)

ti2 = ti %>% select(sid,unpair,unpair_map,unpair_map_hq, unpair_unmap) %>%
    separate(sid, c('sid','gt'), sep='-') %>%
    mutate(gt = str_replace(gt, 'Zmays_', '')) %>%
    mutate(gt = str_replace(gt, '.stat', ''))
ti2 %>% mutate(x=unpair-unpair_map-unpair_unmap) %>% count(x)
tp = ti2 %>% mutate(multi=unpair_map-unpair_map_hq) %>%
    select(sid,gt, unmap=unpair_unmap, multi,uniq=unpair_map_hq) %>%
    inner_join(thf2, by=c('sid')) %>% select(-sid) %>%
    gather(type, cnt, -cond,-gt)

p1 = ggplot(tp) +
    geom_histogram(aes(x=gt,y=cnt, fill=type), stat='identity', position='fill') +
    facet_wrap(.~cond, nrow=1) +
    coord_flip() +
    scale_x_discrete(expand=expansion(mult=c(.01,.01))) +
    scale_fill_simpsons() +
    otheme(legend.pos='top.center.out', legend.dir='h', legend.vjust=-1,
           xtext=F, ytext=T)
fo = glue("{dirw}/03.mapping.stats.pdf")
ggsave(p1, filename = fo, width = 10, height = 6)
#}}}


#{{{ ratio distribution 
fg = '~/projects/genome/data2/syntelog/xref.maize.v5.tsv'
tg = read_tsv(fg) %>% filter(type=='syntelog') %>% select(qry,tgt,gid1,gid2)

fi = glue("{dirw}/50_final/featurecounts.rds")
ti = readRDS(fi) %>%
    separate(sid, c('sid','gt'), sep='-') %>%
    mutate(gt = str_replace(gt, 'Zmays_', ''))
ti2 = ti %>% inner_join(thf2, by=c('sid')) %>% select(cond,gid, gt,cnt)

ti3a = ti2 %>% filter(gt=='B73v5')
#
x = ti2 %>% filter(gt!='B73v5', cnt>0) %>% rename(gid2=gid) %>%
    inner_join(tg, by=c('gt'='qry', 'gid2'='gid2')) %>%
    select(cond, gt, gid=gid1, cnt) %>%
    spread(gt, cnt)
gts = colnames(x)[-c(1:2)]
ti3b = x %>% rowwise() %>%
    mutate(cnt1 = max(c_across(B97:W22), na.rm=T),
           gt1 = gts[which.max(c_across(B97:W22))]) %>%
    ungroup() %>%
    select(cond, gid, gt1, cnt1)
#
ti3 = ti3a %>% rename(gt0=gt, cnt0=cnt) %>%
    inner_join(ti3b, by=c('cond','gid'))

tp = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1))
p1 = ggplot(tp) +
    geom_histogram(aes(x = ratio), bins=50) +
    facet_wrap(~cond, scale='free') +
    scale_x_continuous(expand=expansion(mult=c(.05,.05))) +
    otheme(legend.pos='none', legend.dir='v',
              xtitle=T, ytitle=T, xtext=T, ytext=T)
fo = glue("{dirw}/test.pdf")
ggsave(p1, filename = fo, width = 10, height = 8)
#}}}

y = tp %>% filter(ratio < .05) %>% count(gid) %>%
    filter(n>=8)
tp %>% filter(gid %in% y$gid) %>% count(gt1) %>% arrange(desc(n)) %>% print(n=30)
