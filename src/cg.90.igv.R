source("functions.R")
require(jsonlite)
diro = glue('{dird}/tracks')
dir_s3 = "/datalus/weiyu/projects/s3/zhoup-nfo/zm.cg20a.2"
org = 'Zmays_B73v5'
make_cage_track <- function(sid, name, s3_pre) {
    #{{{
    url1 = glue("{s3_pre}/{sid}.plus.bw")
    url2 = glue("{s3_pre}/{sid}.minus.bw")
    tracks =  tibble(type='wig',format='bigwig',url=c(url1,url2),
                     color=c('red','green'))
    list(name=name, height=50, displayMode="SQUISHED", type='merged',
         tracks = tracks)
    #}}}
}

#{{{ copy bigwig and bam files
tho = tha %>% mutate(sname = glue("{grp}.{str_pad(grp.x,2,pad='0')}")) %>%
    select(grp, grp.x, sname, cond, sid)

tho %>%# slice(1:2) %>%
    mutate(fi=glue("{dir_s3}/34_ctss/{sid}-{org}.plus.bw")) %>%
    mutate(fo=glue("{diro}/01_bigwig/{sname}.plus.bw")) %>%
    mutate(j = map2_lgl(fi, fo, file.copy, overwrite=T)) %>%
    pull(j)
tho %>%# slice(1:2) %>%
    mutate(fi=glue("{dir_s3}/34_ctss/{sid}-{org}.minus.bw")) %>%
    mutate(fo=glue("{diro}/01_bigwig/{sname}.minus.bw")) %>%
    mutate(j = map2_lgl(fi, fo, file.copy, overwrite=T)) %>%
    pull(j)

tho %>%# slice(1:2) %>%
    mutate(fi=glue("{dir_s3}/20_bam/{sid}-{org}.sorted.bam")) %>%
    mutate(fo=glue("{diro}/01_bam/{sname}.bam")) %>%
    mutate(j = map2_lgl(fi, fo, file.copy, overwrite=T)) %>%
    pull(j)
tho %>%# slice(1:2) %>%
    mutate(fi=glue("{dir_s3}/20_bam/{sid}-{org}.sorted.bam.csi")) %>%
    mutate(fo=glue("{diro}/01_bam/{sname}.bam.csi")) %>%
    mutate(j = map2_lgl(fi, fo, file.copy, overwrite=T)) %>%
    pull(j)

to = tho %>% rename(group=grp, i=grp.x, sample=sname, condition=cond)
fo = glue("{diro}/00.meta.tsv")
write_tsv(to, fo)
#}}}

#{{{ merge bigwig and bam
size = glue("/datalus/genomes/{org}/15_intervals/01.chrom.sizes")
merge_bws <- function(fis, fo, size, minus=F) {
    #{{{
    cat(fo, "\n")
    fi = str_c(fis, sep=" ", collapse=" ")
    tag = ifelse(minus, "-threshold=-9999999999", '')
    cmd = glue("bigWigMerge {tag} {fi} tmp.bg")
    system(cmd)
    cmd = glue("bedGraphToBigWig tmp.bg {size} {fo}")
    system(cmd)
    system("rm tmp.bg")
    #}}}
}
merge_bams <- function(fis, fo) {
    #{{{
    cat(fo, "\n")
    fi = str_c(fis, sep=" ", collapse=" ")
    cmd = glue("sambamba merge -t 4 {fo} {fi}")
    system(cmd)
    #system(glue("sambamba index {fo}"))
    #}}}
}
thfo = tho %>% filter(sid %in% thf$sid) %>%
    mutate(bw1=glue("{diro}/01_bigwig/{sname}.plus.bw")) %>%
    mutate(bw2=glue("{diro}/01_bigwig/{sname}.minus.bw")) %>%
    mutate(bam=glue("{diro}/01_bam/{sname}.bam")) %>%
    group_by(grp) %>%
    summarise(bw1=list(bw1), bw2=list(bw2), bam=list(bam)) %>% ungroup() %>%
    mutate(obw1=glue("{diro}/05_bigwig_merged/{grp}.plus.bw")) %>%
    mutate(obw2=glue("{diro}/05_bigwig_merged/{grp}.minus.bw")) %>%
    mutate(obam=glue("{diro}/05_bam_merged/{grp}.bam"))

thfo %>% mutate(j=map2_int(bw1, obw1, merge_bws, size=!!size)) %>% pull(j)
thfo %>% mutate(j=map2_int(bw2, obw2, merge_bws, size=!!size, minus=T)) %>% pull(j)

thfo %>% mutate(j=map2_int(bam, obam, merge_bams)) %>% pull(j)
#}}}

#{{{ [old] IGV track json - raw bigwig
jts = thf %>%
    group_by(Genotype,Tissue) %>% mutate(i = 1:n()) %>% ungroup() %>%
    arrange(Genotype, Tissue) %>%
    mutate(name=glue("{Genotype} {Tissue} rep{i} ({SampleID})")) %>%
    mutate(jt = map2(SampleID, name, make_cage_track, s3_pre="https://s3.msi.umn.edu/zhoup-igv-data/Zmays-B73/cage/raw_bw"))

j = list()
j[['label']] = "CAGE-seq Raw"
j[['description']] = "CAGE-seq experiments"
j[['tracks']] = jts %>% pull(jt)
fo = file.path(diro, 'raw_bw.json')
write(toJSON(j, auto_unbox=T, pretty=T), fo)

to = thf %>%
    group_by(Genotype,Tissue) %>% mutate(i = 1:n()) %>% ungroup() %>%
    arrange(Genotype, Tissue) %>%
    mutate(name=glue("{Genotype} {Tissue} rep{i} ({SampleID})")) %>%
    mutate(fo1 = glue("https://s3.msi.umn.edu/zhoup-igv-data/Zmays-B73/cage/raw_bw/{SampleID}.plus.bw")) %>%
    mutate(fo2 = glue("https://s3.msi.umn.edu/zhoup-igv-data/Zmays-B73/cage/raw_bw/{SampleID}.minus.bw")) %>%
    select(name,Genotype,Tissue,Treatment,Replicate,fo1,fo2)
fo = glue("{dird}/91_share/00.bigwigs.tsv")
write_tsv(to, fo)
#}}}

#{{{ [old] merge replicates + IGV track json
bw_str <- function(sids, srd='plus')
    str_c(str_c(sids, srd, "bw", sep='.'), collapse=" ")
x = thf %>% group_by(Genotype,Tissue) %>%
    summarise(sids = list(SampleID)) %>% ungroup() %>%
    mutate(name = glue("{Genotype} {Tissue}")) %>%
    mutate(str1 = map_chr(sids, bw_str, srd='plus')) %>%
    mutate(str2 = map_chr(sids, bw_str, srd='minus'))
x %>% select(str1)

jts = x %>%
    mutate(name=glue("{Genotype}_{Tissue}")) %>%
    mutate(jt = map2(name, name, make_cage_track, s3_pre="https://s3.msi.umn.edu/zhoup-igv-data/Zmays-B73/cage/merged_bw")) %>%
    pull(jt)

j = list()
j[['label']] = "CAGE-seq merged"
j[['description']] = "CAGE-seq experiments (replicate merged)"
j[['tracks']] = jts
fo = file.path(diro, 'merged_bw.json')
write(toJSON(j, auto_unbox=T, pretty=T), fo)
#}}}


