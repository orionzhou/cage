require(devtools)
load_all('~/git/rmaize')
require(ggpubr)
require(ggforce)
require(knitr)
require(kableExtra)
options(knitr.table.format = "latex")
dirg = '~/projects/genome/data'
dirp = '~/projects/cage'
dird = glue('{dirp}/data')
dirr = glue('{dird}/raw')
dirf = glue('{dird}/95_figures/plots')
gcfg = read_genome_conf()
#
fh = glue('{dird}/01.meta.tsv')
th0 = read_tsv(fh)
tgl = gcfg$gene %>% mutate(loc = glue("{chrom}:{start}-{end}")) %>% select(gid,loc)

#{{{ process th
tissues = c('shoot','root','stem','husk')
gts = c("B73",'LH143',"B73xLH143")
gts = c("B73",'A632',"B73xA632")
notes = c("normal",'cold_control','drought_control','2015','cold','drought','AMG','Maike')
tha = th0 %>% replace_na(list(Treatment='normal')) %>%
    mutate(Genotype = str_replace(Genotype, 'LH143','A632')) %>%
    mutate(Tissue=factor(Tissue, levels=tissues)) %>%
    mutate(Genotype=factor(Genotype, levels=gts)) %>%
    mutate(Treatment=factor(Treatment, levels=notes)) %>%
    arrange(Tissue, Genotype, Treatment) %>%
    mutate(lab = str_c(Tissue,Genotype,Treatment,Replicate, sep='_'))
th = tha %>%
    filter(!Treatment %in% c('2015','cold','drought'))
th %>% dplyr::count(Tissue, Genotype, Treatment) %>% print(n=20)
thf = th %>%
    filter(Tissue %in% c("stem", "husk") | Treatment %in% c("normal",'cold_control','drought_control')) %>%
    filter(!SampleID %in% c("s22",'s24','s25'))
thfs = thf %>% distinct(Genotype, Tissue) %>% arrange(Genotype,Tissue) %>%
    mutate(cond.s = c("Bs","Br",'Bstem','Bhusk','As','Ar','Hs','Hr')) %>%
    mutate(cond.l = str_c(Genotype, Tissue, sep=':')) %>%
    mutate(cond.s=fct_inorder(cond.s), cond.l=fct_inorder(cond.l))
thf = thf %>% inner_join(thfs, by=c('Genotype','Tissue')) %>%
    select(SampleID,Genotype,Tissue,Treatment,Replicate,lab,cond.s,cond.l)
shapes = c('single (1 bp)', 'steep  (2-10 bp)', 'broad  (>= 10 bp)')
shapess = c('single', 'steep', 'broad')
#}}}
#
read_cage_bws <- function(diri, dsg, minSupport=2) { # dsg is a df with 'Name'
    #{{{ read in, remove singletons
    fps = glue("{diri}/{dsg$Name}.plus.bw")
    fms = glue("{diri}/{dsg$Name}.minus.bw")
    rownames(dsg) = dsg$Name
    names(fps) = dsg$Name; names(fms) = dsg$Name
    #
    SE_ctss_raw = quantifyCTSSs(plusStrand=BigWigFileList(fps),
                             minusStrand=BigWigFileList(fms), design=dsg)
    assay(SE_ctss_raw,'counts') = abs(assay(SE_ctss_raw,'counts'))
    SE_ctss_raw %>%
        calcTotalTags(inputAssay="counts", outputColumn="totalTags") %>%
        calcTPM(inputAssay='counts', outputAssay='TPM', totalTags='totalTags') %>%
        calcPooled(inputAssay="TPM") %>%
        calcSupport(inputAssay="counts", outputColumn="support", unexpressed=0) %>%
        subset(support >= minSupport) %>%
        calcTPM(totalTags="totalTags") %>%
        calcPooled()
    #}}}
}
make_tss <- function(SE_ctss, txdb, pooledCutoff=.1, unexpressed=.1, minSamples=2) {
#{{{ TSS candidate
    tc = clusterUnidirectionally(SE_ctss, pooledCutoff=pooledCutoff, mergeDist=20)
    #tc %>% as_tibble %>% filter(seqnames=='B01', start >= 55500, end <= 55800)
    #tc %>% as_tibble %>% filter(seqnames=='B01', start >= 199000, end <= 199200)
    SE_tss_raw = SE_ctss %>% quantifyClusters(clusters=tc, inputAssay='counts') %>%
        calcTPM(totalTags='totalTags') %>%
        calcPooled()
    #
    SE_tss = SE_tss_raw %>%
        #subsetBySupport(inputAssay='TPM',unexpressed=.1,minSamples=2) %>%
        subsetBySupport(inputAssay='TPM',unexpressed=unexpressed,minSamples=minSamples) %>%
        calcTPM(totalTags='totalTags') %>%
        calcPooled() %>%
        calcShape(pooled=SE_ctss, outputColumn='IQR', shapeFunction=shapeIQR,
            lower=.1, upper=.9) %>%
        assignTxID(txModels=txdb, outputColumn='txID') %>%
        assignTxType(txModels=txdb, outputColumn='peakTxType', swap='thick') %>%
        assignGeneID(geneModels=txdb, outputColumn='geneID')
    SE_tss
    #}}}
}
subset_samples <- function(sids, ex) {
    #{{{
    #fctss1 = fctss[,sids] %>% calcPooled(inputAssay = 'TPM')
    x = ex[,sids] %>%
        calcPooled(inputAssay="TPM") %>%
        calcSupport(inputAssay="TPM", outputColumn="support", unexpressed=1)
        #subset(support >= 1) %>%
        #calcShape(pooled=fctss1, outputColumn='IQR', shapeFunction=shapeIQR,
                  #lower=.1, upper=.9)
    nrep = length(sids)
    rowRanges(x) %>% as_tibble() %>%
        mutate(tpm = score/nrep) %>%
        mutate(nrep = nrep) %>%
        mutate(i = 1:n()) %>%
        select(i, tpm, support, nrep)
    #}}}
}
#
download_ctss_bed <- function(fname, fo1, fo2, url_pre, genome) {
    #{{{
    system(glue("wget {url_pre}/{str_replace_all(fname,'%','%25')}"))
    #
    require(CAGEr)
    x = CAGEr:::import.bedCTSS(fname)
    x2 = as.data.frame(x)
    colnames(x2)[1]='chr'
    x3 = as(x2, 'CAGEexp')
    genomeName(x3) = genome
    CAGEr:::exportCTSStoBedGraph(x3, values = "raw", format = "BigWig")
    #
    system(sprintf("mv score.CTSS.raw.plus.bw %s", fo1))
    system(sprintf("mv score.CTSS.raw.minus.bw %s", fo2))
    #
    system(glue("rm {fname}"))
    #}}}
}


