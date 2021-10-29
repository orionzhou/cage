source("functions.R")
#require(CAGEr)
#require(BSgenome.Zmays.B73v5)
#require(CAGEfightR)
gen = "BSgenome.Zmays.B73v5"
dirw = glue('{dird}/01_qc')




#{{{ create CAGEexp object [long time]
#fis = sprintf("%s/34_ctss/%s.bw", dirr, thf$SampleID)
fis = sprintf("%s/20_bam/%s.bam", dirr, thf$SampleID)
ce = CAGEexp(genomeName=gen, inputFiles=fis, inputFilesType="bam",
             sampleLabels = thf$lab)
getCTSS(ce, sequencingQualityThreshold = 10, mappingQualityThreshold = 20,
    removeFirstG = F, correctSystematicG = F, useMulticore = F, nrCores = NULL)
#
fo = file.path(dirw, '01.rds')
saveRDS(ce, fo)
#}}}

fi = file.path(dirw, '01.rds')
ce = readRDS(fi)

#{{{ CTSS count stats
tc = CTSStagCountDF(ce) %>% as_tibble()
tc2 = tc %>% mutate(tss=1:n()) %>% gather(sid, cnt, -tss) %>%
    filter(cnt > 0) %>%
    dplyr::count(sid)

tp = tc2 %>% dplyr::rename(lab=sid) %>% inner_join(th, by='lab') %>%
    mutate(tis_gt = str_c(Tissue, Genotype, sep='_')) %>%
    mutate(cond = str_c(Treatment, Replicate, sep="_")) %>%
    select(tis_gt, cond, n) %>%
    spread(tis_gt, n)
fo = file.path(dirw, '03.ctss.cnt.tsv')
write_tsv(tp, fo, na='')
#}}}

#{{{ sample correlation
tis = 'root'
tis = 'shoot'
gt = 'B73'
gt = 'LH143'
#
labs = th %>% filter(Tissue==tis, Genotype==gt) %>% pull(lab)
fo = sprintf("%s/13.corr.%s.%s.pdf", dirw, tis, gt)
pdf(fo, height=8, width=8)
corr.m <- plotCorrelation2(ce, samples = labs,
    tagCountThreshold = 5, applyThresholdBoth = F, method = "pearson")
dev.off()
#}}}

#{{{ test CAGEr import
fi = '/scratch.global/zhoux379/nf/work/rnaseq/cg18a/10/2941a9ca63561152b79267012c5425/s03.bam'
x = CAGEr:::import.bam.ctss(fi, filetype='bam', sequencingQualityThreshold=10,mappingQualityThreshold=20,removeFirstG=F,correctSystematicG=F,genome=coerceInBSgenome(gr))
x = CAGEr:::import.bam.ctss(fi, filetype='bam', sequencingQualityThreshold=10,mappingQualityThreshold=20,removeFirstG=F,correctSystematicG=F, genome=NULL)

fi = '~/Downloads/Follicle%20Associated%20Epithelium%2c%20pool2.CNhs13211.10262-104D1.mm10.nobarcode.ctss.bed.gz'
x = CAGEr:::import.bedCTSS(fi)

x2 = as.data.frame(x)
colnames(x2)[1]='chr'
x3 = as(x2, 'CAGEexp')
genomeName(x3) = 'BSgenome.Mmusculus.UCSC.mm10'
CAGEr:::exportCTSStoBedGraph(x3, values = "raw", format = "BigWig")
system(glue("mv score.CTSS.raw.plus.bw {fo1}"))
system(glue("mv score.CTSS.raw.minus.bw {fo2}"))
#}}}
