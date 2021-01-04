require(CAGEfightR)
require(GenomicFeatures)
source("functions.R")
#txdb = load_txdb("Zmays_B73")
txdb = load_txdb("Athaliana")

#{{{ read human & mouse - only locations no IQR
fi = file.path(dird, '91_human_mouse', 'refTSS_v3.1_human_coordinate.hg38.bed')
fi = file.path(dird, '91_human_mouse', 'refTSS_v3.1_mouse_coordinate.mm10.bed')
ti = read_tsv(fi, col_names=F) %>% mutate(size=X3-X2, id=1:n()) %>%
    mutate(shape = ifelse(size==1, shapes[1], ifelse(size <= 10, shapes[2], shapes[3]))) %>%
    mutate(shape = factor(shape, levels=shapes))
#}}}

#{{{ Le2020
dirw = glue("{dird}/90_refs/Le2020")
dsg = tibble(Name='wt', Genotype="WT") %>% as.data.frame()

SE_ctss = read_cage_bws(dirw, dsg, minSupport=1)
head(assay(SE_ctss,"TPM"))
table(mcols(SE_ctss)$support)
#}}}

#{{{ Kurihara2018
dirw = glue("{dird}/90_refs/Kurihara2018")
dsg = read_tsv(glue("{dirw}/01.meta.tsv")) %>% rename(Name=SampleID)

diri = glue("{dirw}/04_ctss")
SE_ctss = read_cage_bws(diri, dsg, minSupport=2)
SE_tss = make_tss(SE_ctss, unexpressed=.2, minSamples=2)
tss = rowRanges(SE_tss) %>% as_tibble() %>%
    mutate(i = 1:n()) %>%
    select(i, chrom=seqnames,start,end,strand,pTPM=score,domi=thick.start,
           IQR, support, tid=txID, peakType=peakTxType, gid=geneID)

r1 = list(SE_ctss=SE_ctss, SE_tss=SE_tss, tss=tss)
fo = glue("{dirw}/10.tss.rds")
saveRDS(r1, fo)
#}}}

#{{{ phantom5 - mouse & human - [long time]
# run cg.job.90.R
#{{{ human
dirw = glue("{dird}/90_refs/Phantom5_human")
genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
txdb = makeTxDbFromUCSC(genome="hg38", tablename="refGene")
setwd(dirw)
#}}}

#{{{ mouse
dirw = glue("{dird}/90_refs/Phantom5_mouse")
genomeName = 'BSgenome.Mmusculus.UCSC.mm10'
txdb = makeTxDbFromUCSC(genome="mm10", tablename="refGene")
setwd(dirw)
#}}}

diri = glue("{dirw}/04_ctss")
fh = glue("{dirw}/01.meta.tsv")
th = read_tsv(fh)

#{{{ create CAGEfightR object [large mem]
dsg = th %>% rename(Name=sid)

SE_ctss = read_cage_bws(diri, dsg, minSupport=2)
SE_tss = make_tss(SE_ctss, txdb, pooledCutoff=1, unexpressed=1, minSamples=10)
tss = rowRanges(SE_tss) %>% as_tibble() %>%
    mutate(i = 1:n()) %>%
    select(i, chrom=seqnames,start,end,strand,pTPM=score,domi=thick.start,
           IQR, support, tid=txID, peakType=peakTxType, gid=geneID)

r1 = list(SE_tss=SE_tss, tss=tss)
fo = glue("{dirw}/10.tss.rds")
saveRDS(r1, fo)
#}}}
#}}}

