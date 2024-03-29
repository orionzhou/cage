source("functions.R")
dirw = glue("{dird}/03_qc")
#
require(BiocParallel)
ncpu = 8
bpparam <- MulticoreParam(workers=ncpu)

fi = glue('{dirw}/03.rep.merged.rds')
r3 = readRDS(fi)
rtss=r3$rtss; tss=r3$tss; gtss=r3$gtss; SE_ctss_m=r3$SE_ctss_m; ctss_m=r3$ctss_m

#{{{ functions to calc shape metrics
require(kSamples)
#
freq2sample <- function(cnts) rep(1:length(cnts), cnts)
shape.test <- function(xs, method='qn') {
    #{{{
    xs2 = lapply(xs, freq2sample)
    if(method == 'ad') {
        ad.test(xs2)$ad %>% as_tibble() %>% pluck(3,1)
    } else if(method == 'qn') {
        qn.test(xs2)$qn %>% as_tibble() %>% pluck(1,2)
    } else if(method == 'ks') {
        stopifnot(length(xs2) == 2)
        ks.test(xs2[[1]], xs2[[2]])$p.value
    }
    #}}}
}
shifting_score <- function(xs) {
    #{{{
    stopifnot(length(xs) == 2)
    x1 = xs[[1]]; x2 = xs[[2]]
    if(sum(x1) > sum(x2)) { x0 = x1; x1 = x2; x2 = x0 }
    x1r = rev(x1); x2r = rev(x2)
    s1 = (max(cumsum(x1) - cumsum(x2))) / sum(x1)
    s2 = (max(cumsum(x1r) - cumsum(x2r))) / sum(x1r)
    max(s1, s2)
    #}}}
}
ks_dist <- function(xs) {
    #{{{
    xs2 = lapply(xs, freq2sample)
    stopifnot(length(xs2) == 2)
    ks.test(xs2[[1]], xs2[[2]])$statistic[['D']]
    #}}}
}
cmp_shift_score <- function(ti) {
    #{{{
    ti1a = ti %>% select(i, npos, tpm.grp) %>% unnest(tpm.grp)
    ti1b = ti %>% select(i, cnts) %>% unnest(cnts)
    ti2 = ti1a %>% inner_join(ti1b, by=c('i','grp'))
    #{{{ multi-sample test
    grps4 = c("shoot.B","shoot.A",'root.B','root.A')
    r2 = ti2 %>% filter(npos > 1) %>%
        filter(tpm >= 1, grp %in% grps4) %>%
        group_by(i) %>%
        summarise(ngrp = n(), cnts = list(cnts)) %>%
        ungroup() %>%
        filter(ngrp >= 2)
    #y <- bplapply(r2$cnts, shape.test, method='qn', BPPARAM = bpparam)
    #r.multi = r2 %>% mutate(pval.qn.raw = unlist(y)) %>% select(i,pval.qn.raw)
    #}}}
    #
    #{{{ 2-sample comparison
    r3 = ti2 %>% filter(npos > 1) %>%
        select(i, grp, tpm, cnts)
    imix <- function(tpmA, tpmB, cntsA, cntsB)
        tpmA / (tpmA+tpmB) * cntsA + tpmB / (tpmA+tpmB) * cntsB
    r3a = r3 %>% pivot_wider(names_from=grp, values_from=c(tpm,cnts)) %>%
        mutate(tpm_root.C = (tpm_root.A + tpm_root.A) / 2) %>%
        mutate(tpm_shoot.C = (tpm_shoot.B + tpm_shoot.A) / 2) %>%
        mutate(cnts_root.C = pmap(list(tpm_root.A, tpm_root.B, cnts_root.A, cnts_root.B), imix)) %>%
        mutate(cnts_shoot.C = pmap(list(tpm_shoot.A, tpm_shoot.B, cnts_shoot.A, cnts_shoot.B), imix)) %>%
        select(i, ends_with('root.C'), ends_with('shoot.C')) %>%
        pivot_longer(cols=!i, names_to=c(".value","grp"), names_pattern="(.*)_(.*)")
    r3 = r3 %>% bind_rows(r3a)
    r3b = cmps %>%
        inner_join(r3, by=c('grp1'='grp')) %>% rename(tpm1=tpm, cnts1=cnts) %>%
        inner_join(r3, by=c('i','grp2'='grp')) %>% rename(tpm2=tpm, cnts2=cnts) %>%
        filter(tpm1 >= 1, tpm2 >= 1) %>%
        mutate(cnts = map2(cnts1, cnts2, list))
    y1 <- bplapply(r3b$cnts, ks_dist, BPPARAM = bpparam)
    #y1 <- bplapply(r3b$cnts, shifting_score, BPPARAM = bpparam)
    #y2 <- bplapply(r3b$cnts, shape.test, method='ks', BPPARAM = bpparam)
    #y3 <- bplapply(r3b$cnts, shape.test, method='qn', BPPARAM = bpparam)
    #
    r.pair = r3b %>% select(i,grp1,grp2) %>%
        mutate(shift = unlist(y1)) %>%
        #mutate(pval.ks.raw = unlist(y2)) %>%
        #mutate(pval.qn.raw = unlist(y3)) %>%
        group_by(i) %>% nest() %>% rename(cmp = data)
    #}}}
    #ti2s = ti2 %>% group_by(i) %>%
        #summarise(ncond4 = sum(tpm>=1 & cond %in% conds4)) %>% ungroup()
    r = ti %>%# inner_join(ti2s, by='i') %>%
        #left_join(r.multi, by='i') %>%
        left_join(r.pair, by='i')
    r
    #}}}
}
#}}}

ti = rtss %>% rename(i=tidx)
rtss2 = cmp_shift_score(ti) %>%
    select(i,tpm,IQR,entropy,support,peakType,gid,npos,tpm.grp,cmp)

ti = gtss %>% rename(i=gidx)
gtss2 = cmp_shift_score(ti) %>% select(i,tpm,gid,n_tss,npos,cmp)

r5 = list(rtss=rtss2, gtss=gtss2)
fo = glue("{dirw}/05.shift.rds")
saveRDS(r5, fo)

#{{{ ### benchmark multidplyr and bpparallel
#require(multidplyr)
#cl <- new_cluster(ncpu)
#cluster_library(cl, "tidyverse")
#cluster_library(cl, "kSamples")
#cluster_copy(cl, 'qn.test4')
#cluster_copy(cl, 'freq2sample')

#r2 = ttss %>% filter(npos >= 2) %>% slice(1:5000)
#system.time(y <- lapply(r2$cnts, qn.test4))
#system.time( y <- r2 %>%
    #partition(cl) %>%
    #mutate(tsr = map_dbl(cnts, qn.test4)) %>%
    #select(i,tsr) %>%
    #collect()
#)
#system.time(y <- bplapply(r2$cnts, qn.test4, BPPARAM = bpparam))
#}}}

