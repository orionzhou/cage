source("functions.R")
dirw = glue("{dird}/05_nam")
setwd(dirw)

#{{{ preprocess read counts
fi = glue("{dirw}/50_final/featurecounts.uniq.rds")
tia = readRDS(fi) %>% rename(sid0 = sid) %>%
    mutate(sid = str_sub(sid0, 1, 3)) %>%
    mutate(gt = str_sub(sid0, 11, -6)) %>%
    select(sid, gt, gid, cnt)
tia %>% distinct(gt) %>% print(n=30)
fi = glue("{dirw}/50_final/featurecounts.multi.rds")
tib = readRDS(fi) %>% rename(sid0 = sid) %>%
    mutate(sid = str_sub(sid0, 1, 3)) %>%
    mutate(gt = str_sub(sid0, 11, -7)) %>%
    select(sid, gt, gid, cnt)
tib %>% distinct(gt) %>% print(n=30)

rc = tia %>% rename(uniq = cnt) %>%
    inner_join(tib %>% rename(multi=cnt), by=c('sid','gt','gid'))

fo = glue("{dirw}/00.rds")
saveRDS(rc, fo)
#}}}


#{{{ ratio distribution 
fg = '~/projects/genome/data2/syntelog/xref.maize.v5.tsv'
tg = read_tsv(fg) %>% filter(type=='syntelog') %>% select(qry,tgt,gid1,gid2)

fi = glue("{dirw}/00.rds")
rc = readRDS(fi)

ti2 = rc %>% select(sid,gid, gt, cnt=uniq)
ti3a = ti2 %>% filter(gt=='B73v5')
#
x = ti2 %>% filter(gt!='B73v5', cnt>0) %>% rename(gid2=gid) %>%
    inner_join(tg, by=c('gt'='qry', 'gid2'='gid2')) %>%
    select(sid, gt, gid=gid1, cnt) %>%
    spread(gt, cnt)
gts = colnames(x)[-c(1:2)]
ti3b = x %>% rowwise() %>%
    mutate(cnt1 = max(c_across(B97:W22), na.rm=T),
           gt1 = gts[which.max(c_across(B97:W22))]) %>%
    ungroup() %>%
    select(sid, gid, gt1, cnt1)
#
ti3 = ti3a %>% rename(gt0=gt, cnt0=cnt) %>%
    inner_join(ti3b, by=c('sid','gid'))

#{{{ per-rep ratio distri
tp = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    inner_join(thf2, by='sid') %>%
    filter(! str_detect(sname, "stem") & !str_detect(sname, 'husk')) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1))
p1 = ggplot(tp) +
    geom_histogram(aes(x = ratio), bins=50) +
    facet_wrap(~sname, scale='free_x') +
    scale_x_continuous(breaks=c(.25,.5,.75),expand=expansion(mult=c(.03,.03))) +
    otheme(legend.pos='none', legend.dir='v', panel.border=F,
              xtitle=T, ytitle=T, xtext=T, xtick=T, ygrid=T, ytext=T, ytick=T)
#}}}
fo = glue("{dirw}/11.rep.ratio.pdf")
ggsave(p1, filename = fo, width = 10, height = 8)

#{{{ stem & husk
tp = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    inner_join(thf2, by='sid') %>%
    filter(str_detect(sname, "stem") | str_detect(sname, 'husk')) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1))
p1 = ggplot(tp) +
    geom_histogram(aes(x = ratio), bins=50) +
    facet_wrap(~sname, scale='free_x') +
    scale_x_continuous(breaks=c(.25,.5,.75),expand=expansion(mult=c(.03,.03))) +
    otheme(legend.pos='none', legend.dir='v', panel.border=F,
              xtitle=T, ytitle=T, xtext=T, xtick=T, ygrid=T, ytext=T, ytick=T)
#}}}
fo = glue("{dirw}/11.rep.ratio.stem.husk.pdf")
ggsave(p1, filename = fo, width = 5, height = 3)

#{{{ output per-rep ratios
cnames = thf2 %>%
    #filter(! str_detect(sname, "stem") & !str_detect(sname, 'husk')) %>%
    select(sname) %>% crossing(suf=c('cnt0','cnt1','gt1','ratio')) %>%
    mutate(cname = glue("{sname}_{suf}")) %>% arrange(sname) %>% pull(cname)
tt = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    inner_join(thf2, by='sid') %>%
    #filter(! str_detect(sname, "stem") & !str_detect(sname, 'husk')) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1)) %>%
    select(gid, sname, cnt0, cnt1, gt1, ratio) %>%
    pivot_wider(names_from=sname, values_from=c(cnt0,cnt1,gt1,ratio),
        names_glue="{sname}_{.value}") %>%
    select(gid, eval(cnames))
#}}}
fo = glue("{dirw}/05.ratio.rep.tsv")
write_tsv(tt, fo, na='')
#}}}

#{{{ rep-merged ratio
ti2m = ti2 %>% inner_join(thf2, by='sid') %>%
    group_by(gid, grp, gt) %>% summarise(cnt=sum(cnt)) %>% ungroup()
ti3a = ti2m %>% filter(gt=='B73v5')
#
x = ti2m %>% filter(gt!='B73v5', cnt>0) %>% rename(gid2=gid) %>%
    inner_join(tg, by=c('gt'='qry', 'gid2'='gid2')) %>%
    select(grp, gt, gid=gid1, cnt) %>%
    spread(gt, cnt)
gts = colnames(x)[-c(1:2)]
ti3b = x %>% rowwise() %>%
    mutate(cnt1 = max(c_across(B97:W22), na.rm=T),
           gt1 = gts[which.max(c_across(B97:W22))]) %>%
    ungroup() %>%
    select(grp, gid, gt1, cnt1)
#
ti3 = ti3a %>% rename(gt0=gt, cnt0=cnt) %>%
    inner_join(ti3b, by=c('grp','gid'))

#{{{ per-rep ratio distri
tp = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    #filter(! str_detect(grp, "stem") & !str_detect(grp, 'husk')) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1))
p1 = ggplot(tp) +
    geom_histogram(aes(x = ratio), bins=50) +
    facet_wrap(~grp, scale='free_x') +
    scale_x_continuous(breaks=c(.25,.5,.75),expand=expansion(mult=c(.03,.03))) +
    otheme(legend.pos='none', legend.dir='v', panel.border=F,
              xtitle=T, ytitle=T, xtext=T, xtick=T, ygrid=T, ytext=T, ytick=T)
#}}}
fo = glue("{dirw}/12.ratio.pdf")
ggsave(p1, filename = fo, width = 6, height = 4)

#{{{ density plot
require(hexbin)
x = thf %>% select(grp=cond.s,gt=Genotype,tis=Tissue) %>% distinct(grp,gt,tis)
tx = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1)) %>%
    select(gid,grp,ratio) %>% spread(grp,ratio) %>%
    select(gid,Bs,As,Br,Ar) %>%
    gather(grp, ratio, -gid) %>%
    inner_join(x, by='grp') %>% select(-grp) %>%
    spread(gt, ratio) %>% filter(!is.na(B73), !is.na(A632))

tpl = tx %>% filter(B73>=.9, A632 <=.1) %>% count(tis)
tp = tx %>% rename(x=B73, y=A632)
xlab = 'B73 ratio'; ylab = 'A632 ratio'
p = ggplot(tp, aes(x=x,y=y)) +
    geom_hex(bins=30) +
    geom_rect(aes(xmin=.9,xmax=1,ymin=0,ymax=.1), color='black',size=.3, fill=NA) +
    geom_text(data=tpl,aes(x=.95,y=.05,label=n), color='black',size=2) +
    scale_x_continuous(name=xlab,expand=expansion(mult=c(.02,.02))) +
    scale_y_continuous(name=ylab,expand=expansion(mult=c(.02,.02))) +
    scale_fill_viridis(name='density', direction=-1) +
    #scale_fill_manual(name='', values=pal_npg()(3)[c(3,2,1)]) +
    facet_wrap(tis~., nrow=2) +
    otheme(legend.pos='right', legend.title=T, legend.dir='v',
           legend.box='v', panel.spacing=.3,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) + o_margin(.3,.3,.3,.3)
fo = glue("{dirw}/13.ratio.density.pdf")
ggsave(p, file=fo, width=5, height=8)
#}}}


#{{{ output rep-merged ratios
cnames = thf2 %>%
    distinct(grp)  %>% crossing(suf=c('cnt0','cnt1','gt1','ratio')) %>%
    mutate(cname = glue("{grp}_{suf}")) %>% arrange(grp) %>% pull(cname)
tt = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1)) %>%
    select(gid, grp, cnt0, cnt1, gt1, ratio) %>%
    pivot_wider(names_from=grp, values_from=c(cnt0,cnt1,gt1,ratio),
        names_glue="{grp}_{.value}") %>%
    select(gid, eval(cnames))
#}}}
fo = glue("{dirw}/06.ratio.tsv")
write_tsv(tt, fo, na='')
#}}}




