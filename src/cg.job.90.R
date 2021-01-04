require(CAGEfightR)
require(GenomicFeatures)
source("functions.R")

#{{{ phantom5 - mouse & human
#{{{ mouse
dirw = glue("{dird}/90_refs/Phantom5_human")
url_pre='https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/basic/human.tissue.hCAGE'
genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
setwd(dirw)
#}}}
#{{{ human
#dirw = glue("{dird}/90_refs/Phantom5_mouse")
#url_pre="https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/basic/mouse.tissue.hCAGE"
#genomeName = 'BSgenome.Mmusculus.UCSC.mm10'
#setwd(dirw)
#}}}

diro = glue("{dirw}/04_ctss")
fh = glue("{dirw}/00.txt")
th = read_tsv(fh) %>% select(name=3,fname=25) %>%
    mutate(sid=glue("s{str_pad(1:n(), 3, side='left', pad='0')}")) %>%
    select(sid, name, fname)

fo = glue("{dirw}/01.meta.tsv")
write_tsv(th %>% select(-fname), fo)

th2 = th %>%
    mutate(fname = str_replace(fname, ".bam", ".ctss.bed.gz")) %>%
    mutate(fo1=glue("{diro}/{sid}.plus.bw")) %>%
    mutate(fo2=glue("{diro}/{sid}.minus.bw")) %>%
    mutate(x = pmap(list(fname,fo1,fo2), download_ctss_bed, url_pre=url_pre,
                    genome=genomeName))
#}}}

