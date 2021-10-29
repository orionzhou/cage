#require(CAGEfightR)
#require(GenomicFeatures)
source("functions.R")
dirw = glue("{dird}/06_snp")

#{{{ create TSS bed 
fi = glue("{dird}/03_qc/02.tss.gtss.B73.rds")
r2 = readRDS(fi)
SE_ctss=r2$SE_ctss; SE_rtss=r2$SE_rtss; SE_tss=r2$SE_tss; SE_gtss=r2$SE_gtss;
rtss=r2$rtss; tss=r2$tss; gtss=r2$gtss

tag_size=20
tb = rtss %>% select(chrom,start,end,tidx, srd) %>%
    mutate(start=start-1) %>%
    mutate(start = ifelse(srd=='-', start-tag_size, start)) %>%
    mutate(end = ifelse(srd=='-', end, end+tag_size)) %>% select(-srd)
fo = glue("{dirw}/01.tss.bed")
write_tsv(tb, fo, col_names=F)
#}}}

#{{{ characterize variants in TSSs 
fi1 = glue("{dirw}/04.ovlp.cov.bed")
ti1 = read_tsv(fi1, col_names=F) %>%
    select(chrom=1,start=2,end=3,tid=4,chrom2=5,start2=6,end2=7,GT=8,DP=9,GQ=10)
ti1a = ti1 %>% filter(!is.na(DP), !is.na(GQ)) %>% mutate(size=end-start) %>%
    group_by(tid,size) %>%
    summarise(base_LQ=sum(DP<10 & GQ<30)) %>% ungroup()


fi2 = glue("{dirw}/04.ovlp.vnt.bed")
ti2 = read_tsv(fi2, col_names=F) %>%
    select(chrom=1,start=2,end=3,tid=4,chrom2=5,start2=6,end2=7,ref=8,alt=9,GT=10,DP=11,GQ=12)
ti2a = ti2 %>% filter(chrom2=='.') %>% select(tid)

vnt = ti2 %>% filter(chrom2!='.') %>%
    mutate(DP=as.integer(DP), GQ=as.integer(GQ)) %>%
    mutate(vType=ifelse(nchar(ref)==1 & nchar(alt)==1, 'snp', 'indel')) %>%
    mutate(qual=ifelse(GQ >= 30 & GT %in% c('1/1','1|1'), 'HQ', 'LQ'))
ti2b = vnt %>% mutate(type=glue("{vType}_{qual}")) %>%
    select(tid, type) %>% count(tid, type) %>%
    spread(type, n)
ti2ab = ti2a %>% bind_rows(ti2b) %>%
    replace_na(list(snp_HQ=0, snp_LQ=0, indel_HQ=0, indel_LQ=0))

tv = ti1a %>% inner_join(ti2ab, by='tid') %>%
    mutate(stype = 'other') %>%
    mutate(stype=ifelse(base_LQ==0 & indel_HQ+indel_LQ+snp_HQ+snp_LQ==0, 'identical', stype)) %>%
    mutate(stype=ifelse(base_LQ==0 & indel_HQ+indel_LQ+snp_LQ==0 & snp_HQ==1, '1_snp', stype)) %>%
    mutate(stype=ifelse(base_LQ==0 & indel_HQ+indel_LQ+snp_LQ==0 & snp_HQ==2, '2_snp', stype))
tv %>% count(stype)
#}}}

fo = glue("{dirw}/10.tss.vnt.rds")
r = list(tv=tv, vnt=vnt)
saveRDS(r, fo)

#{{{ validate with IBD segments 
fi = glue("{dirw}/01.tss.bed")
tb = read_tsv(fi, col_names=c('chrom','start','end','tid'))
fi = glue("{dirw}/10.tss.vnt.rds")
tv = readRDS(fi)$tv

fi = glue("{dirw}/22.ibd.bed")
ti = read_tsv(fi, col_names=c('chrom','start','end','chrom2','start2','end2','tid','olen'))

tv %>% mutate(ibd = ifelse(tid %in% ti$tid, 'IBD', 'non-IBD')) %>%
    count(ibd, stype) %>%
    group_by(ibd) %>%
    mutate(lab=glue("{n} ({percent(n/sum(n))})")) %>% ungroup() %>%
    select(ibd, stype, lab) %>% spread(ibd, lab)

ti %>% inner_join(tv %>% select(tid, stype), by='tid') %>%
    count(chrom,start,end, stype) %>%
    group_by(chrom,start,end) %>%
    mutate(lab=glue("{n} ({percent(n/sum(n))})")) %>% ungroup() %>%
    mutate(size = (end-start)/1e6) %>%
    select(chrom,start,end,size, stype,lab) %>% spread(stype, lab) %>% print(n=42)

tv2 = tv %>% mutate(ibd = ifelse(tid %in% ti$tid, 'IBD', 'non-IBD'))
fo = glue("{dirw}/11.tss.vnt.rds")
r2 = list(tv=tv2, vnt=vnt)
saveRDS(r2, fo)
#}}}

