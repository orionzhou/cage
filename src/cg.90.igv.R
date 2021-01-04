source("functions.R")
require(jsonlite)
diro = file.path(dird, 'tracks')
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

#{{{ IGV track json - raw bigwig
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
#}}}

#{{{ merge replicates + IGV track json
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


