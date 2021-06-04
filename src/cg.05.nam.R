source("functions.R")
dirw = glue("{dird}/05_nam")
setwd(dirw)

fg = '~/projects/genome/data2/syntelog/xref.maize.v5.tsv'
tg = read_tsv(fg) %>% select(qry,tgt,gid1,gid2)

read_fc <- function(fi, gt, tg) {
    #{{{
    if(file.exists(fi)) {
        ti = read_tsv(fi, skip=1) %>%
            select(gid=1,rc=7)
        if (gt == 'B73v5')
            ti
        else {
            ti %>% rename(gid2=gid) %>%
                inner_join(tg %>% filter(qry==gt) %>% select(gid1, gid2), by='gid2') %>%
                select(gid=gid1, rc)
        }
    }
    else
        NULL
    #}}}
}

fi = glue("{dirw}/50_final/featurecounts.rds")
ti = readRDS(fi) %>%
    separate(sid, c('sid','gt'), sep='-') %>%
    mutate(gt = str_replace(gt, 'Zmays_', ''))

thf2 = thf %>% mutate(cond=glue("{cond.l}.{Replicate}")) %>% select(sid=SampleID,cond)
ti2 = ti %>% inner_join(thf2, by=c('sid')) %>% select(-sid)

ti3a = ti2 %>% filter(gt=='B73v5')
ti3b = ti2 %>% filter(gt!='B73v5', cnt>0) %>% rename(gid2=gid) %>%
    inner_join(tg, by=c('gt'='qry', 'gid2'='gid2')) %>%
    select(cond, gt, gid=gid1, cnt) %>%
    arrange(cond,gid, desc(cnt)) #%>%
    group_by(cond, gid) %>% slice(1) %>% ungroup() %>%
    rename(gt1=gt, cnt1=cnt)

ti3 = ti3a %>% inner_join(ti3b, by=c('cond','gid')) %>%
    rename(gt0=gt, cnt0=cnt)

tp = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1))
p1 = ggplot(tp) +
    geom_histogram(aes(x = ratio), bins=50) +
    facet_wrap(~cond, scale='free') +
    scale_x_continuous(limits=c(-.1,1.1)) +
    otheme(legend.pos='none', legend.dir='v',
              xtitle=T, ytitle=T, xtext=T, ytext=T)
fo = glue("{dirw}/test.pdf")
ggsave(p1, filename = fo, width = 10, height = 8)


