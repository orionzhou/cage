source("functions.R")
dirw = glue("{dird}/85_nanopore")
fh = glue("{dirw}/00.meta.tsv")
th = read_tsv(fh) %>% select(sid=1,tis=2,gt=3,rep=5) %>%
    mutate(lab = glue("{tis}.{gt}.{rep}")) %>% mutate(lab = as_factor(lab)) %>%
    mutate(grp = glue("{tis}.{gt}")) %>% mutate(grp = as_factor(grp))

aln = 'hisat2'
aln = 'star'
diri = glue("{dirw}_data/21_{aln}")
#{{{ mapping stats
fi = glue("{diri}/bamstats.tsv")
tm = read_tsv(fi)
#
ti2 = tm %>% select(sid,unpair,unpair_map,unpair_map0, unpair_map_hq, unpair_map_hq0, unpair_unmap) %>%
    separate(sid, c('sid','gt'), sep='-') %>%
    mutate(gt = str_replace(gt, 'Zmays_', '')) %>%
    mutate(gt = str_replace(gt, '.stat', ''))
ti2 %>% mutate(x=unpair-unpair_map-unpair_unmap) %>% count(x)
types=c('unmap','multi_mm','multi_pm','uniq_mm','uniq_pm')
tp = ti2 %>% mutate(multi=unpair_map-unpair_map_hq) %>%
    mutate(multi_pm=unpair_map0-unpair_map_hq0, multi_mm=multi-multi_pm) %>%
    mutate(uniq_pm=unpair_map_hq0, uniq_mm=unpair_map_hq-uniq_pm) %>%
    select(sid, unmap=unpair_unmap, multi_pm,multi_mm,uniq_pm,uniq_mm) %>%
    inner_join(th %>% select(sid,lab), by=c('sid')) %>%
    select(-sid) %>% gather(type, cnt, -lab) %>%
    mutate(type=factor(type,levels=types)) %>%
    select(tag1=lab, tag2=type, n=cnt) %>% mutate(n=n/1e6)

cols5 = pal_simpsons()(8)[c(3,1,2,4,5)]
p1 = cmp_proportion1(tp, xangle=0, alph=.4, acc=.1,
    lab.size=2, oneline=F,  fills=cols5, ypos='left',
    expand.x=c(.02,.02), expand.y=c(.01,.04), legend.title='',
    legend.pos='none', legend.dir='h', legend.vjust=-1,
    xtext=T, xtick=T, ytext=T, ytick=T,
    margin= c(.2,.1,.1,2.1)) 
fo = glue("{dirw}/22.map.rate.{aln}.pdf")
ggsave(p1, filename = fo, width = 6, height = 4)
#}}}


#{{{ preprocess read counts
fi = glue("{diri}/featurecounts.uniq.rds")
tia = readRDS(fi) %>% rename(sid0 = sid) %>%
    mutate(sid = str_sub(sid0, 1, 3)) %>%
    mutate(gt = str_sub(sid0, 11, -6)) %>%
    select(sid, gt, gid, cnt) %>% filter(gt != "B73")
tia %>% distinct(gt) %>% print(n=30)
fi = glue("{diri}/featurecounts.multi.rds")
tib = readRDS(fi) %>% rename(sid0 = sid) %>%
    mutate(sid = str_sub(sid0, 1, 3)) %>%
    mutate(gt = str_sub(sid0, 11, -7)) %>%
    select(sid, gt, gid, cnt) %>% filter(gt != "B73")
tib %>% distinct(gt) %>% print(n=30)

rc = tia %>% rename(uniq = cnt) %>%
    inner_join(tib %>% rename(multi=cnt), by=c('sid','gt','gid'))

fo = glue("{dirw}/00.rds")
saveRDS(rc, fo)
#}}}

#{{{ ratio distribution 
fg = '~/projects/genome/data2/syntelog/xref.maize.v5.tsv'
tg = read_tsv(fg) %>% filter(type=='syntelog') %>% select(qry,tgt,gid1,gid2) %>%
    mutate(qry=str_replace(qry,'^Zmays_','')) %>%
    mutate(tgt=str_replace(tgt,'^Zmays_',''))

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
    mutate(cnt1 = max(c_across(CML103:W22), na.rm=T),
           gt1 = gts[which.max(c_across(CML103:W22))]) %>%
    ungroup() %>%
    select(sid, gid, gt1, cnt1)
#
ti3 = ti3a %>% rename(gt0=gt, cnt0=cnt) %>%
    inner_join(ti3b, by=c('sid','gid'))

#{{{ per-rep ratio dist
tp = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    inner_join(th, by='sid') %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1))
p1 = ggplot(tp) +
    geom_histogram(aes(x = ratio), bins=50) +
    facet_wrap(~lab, scale='free_x', ncol=3) +
    scale_x_continuous(breaks=c(.25,.5,.75),expand=expansion(mult=c(.03,.03))) +
    scale_y_continuous(name='num. genes', expand=expansion(mult=c(.01,.02))) +
    otheme(legend.pos='none', legend.dir='v', panel.border=F,
              xtitle=T, ytitle=T, xtext=T, xtick=T, ygrid=T, ytext=T, ytick=T)
#}}}
fo = glue("{dirw}/31.rep.ratio.2.pdf")
ggsave(p1, filename = fo, width = 8, height = 6)

#{{{ write per-rep ratios
cnames = th %>% filter(row_number() <=2 ) %>%
    select(lab) %>% crossing(suf=c('cnt0','cnt1','gt1','ratio')) %>%
    mutate(cname = glue("{lab}_{suf}")) %>% arrange(lab) %>% pull(cname)
tt = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    inner_join(th, by='sid') %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1)) %>%
    select(gid, lab, cnt0, cnt1, gt1, ratio) %>%
    pivot_wider(names_from=lab, values_from=c(cnt0,cnt1,gt1,ratio),
        names_glue="{lab}_{.value}") %>%
    select(gid, eval(cnames))
#}}}
fo = glue("{dirw}/35.ratio.rep.2.tsv")
write_tsv(tt, fo, na='')
#}}}

#{{{ rep-merged ratio
ti2m = ti2 %>% inner_join(th %>% select(sid,grp), by='sid') %>%
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

#{{{ ratio distr
tp = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1))
p1 = ggplot(tp) +
    geom_histogram(aes(x = ratio), bins=50) +
    facet_wrap(~grp, scale='free_x', ncol=2) +
    scale_x_continuous(breaks=c(.25,.5,.75),expand=expansion(mult=c(.03,.03))) +
    scale_y_continuous(name='num. genes', expand=expansion(mult=c(.01,.02))) +
    otheme(legend.pos='none', legend.dir='v', panel.border=F,
              xtitle=T, ytitle=T, xtext=T, xtick=T, ygrid=T, ytext=T, ytick=T)
#}}}
fo = glue("{dirw}/36.ratio.pdf")
ggsave(p1, filename = fo, width = 6, height = 3)

#{{{ density plot
require(hexbin)
x = th %>% select(grp,gt) %>% distinct(grp,gt)
tx = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1)) %>%
    select(gid,grp,ratio) %>% spread(grp,ratio) %>%
    select(gid,shoot.B73,shoot.A632) %>%
    gather(grp, ratio, -gid) %>%
    inner_join(x, by='grp') %>% select(-grp) %>%
    spread(gt, ratio) %>% filter(!is.na(B73), !is.na(A632))

tpl = tx %>% filter(B73>=.9, A632 <=.1) %>% count()
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
    otheme(legend.pos='right', legend.title=T, legend.dir='v',
           legend.box='v', panel.spacing=.3,
           xtext=T, xtick=T, xtitle=T, ytitle=T, ytext=T, ytick=T,
           xgrid=T, ygrid=T) + o_margin(.3,.3,.3,.3)
fo = glue("{dirw}/37.ratio.density.pdf")
ggsave(p, file=fo, width=5, height=4)
#}}}

#{{{ write rep-merged ratios
cnames = th %>%
    distinct(grp) %>% crossing(suf=c('cnt0','cnt1','gt1','ratio')) %>%
    mutate(cname = glue("{grp}_{suf}")) %>% arrange(grp) %>% pull(cname)
tt = ti3 %>% filter(cnt0 >= 10 | cnt1 >= 10) %>%
    mutate(ratio = cnt0/pmax(cnt0, cnt1)) %>%
    select(gid, grp, cnt0, cnt1, gt1, ratio) %>%
    pivot_wider(names_from=grp, values_from=c(cnt0,cnt1,gt1,ratio),
        names_glue="{grp}_{.value}") %>%
    select(gid, eval(cnames))
#}}}
fo = glue("{dirw}/38.ratio.tsv")
write_tsv(tt, fo, na='')
#}}}


gids0 = tt %>% filter(shoot.B73.1_ratio < .01 | shoot.B73.2_ratio < .01) %>%
    pull(gid)

x2 = tg %>% count(qry, gid2)
x3 = x2 %>% filter(n>1)
x4 = tg %>% inner_join(x3, by=c('qry','gid2')) %>%
    distinct(gid1) %>% pull(gid1)

sum(gids0 %in% x4)
