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
fh = glue('{dird}/00.meta.tsv')
th0 = read_tsv(fh)
tgl = gcfg$gene %>% mutate(loc = glue("{chrom}:{start}-{end}")) %>%
    select(gid,ttype,loc)
peakTypes = c("promoter","proximal",'fiveUTR')

#{{{ process th
tiss = c('shoot','root','stem','husk')
gts = c("B73",'A632',"F1")
conds = c("r",'cold.ck','drought.ck','2015','cold','drought','AMG')
tha = th0 %>%
    rename(sid=SampleID, tis=Tissue, gt=Genotype, cond=Treatment, rep=Replicate) %>%
    mutate(gt = ifelse(gt=='B73xA632', 'F1', gt)) %>%
    mutate(tis = factor(tis, levels=tiss)) %>%
    mutate(gt = factor(gt, levels=gts)) %>%
    replace_na(list(cond='r')) %>%
    mutate(cond=str_replace(cond, "_control", ".ck")) %>%
    mutate(cond=ifelse(cond == 'Maike', 'r', cond)) %>%
    mutate(cond = factor(cond, levels=conds)) %>%
    select(sid,tis,gt,cond,rep,spots,avgLength) %>%
    arrange(tis, gt, cond, rep) %>%
    mutate(pnl = glue("{tis}.{gt}.{cond}.{rep}")) %>%
    mutate(pnl = as_factor(pnl)) %>%
    mutate(grp = glue("{tis}.{gt}"), grp.xlab = glue("{cond}.{rep}")) %>%
    mutate(grp = as_factor(grp)) %>%
    group_by(grp) %>% mutate(grp.x = 1:n()) %>% ungroup()
th = tha %>% filter(!cond %in% c('2015','cold','drought'))
th %>% dplyr::count(tis, gt, cond) %>% print(n=20)
gt_map = c("B73"='B','A632'='A','F1'='H')
thf = th %>%
    filter(tis %in% c("shoot", "root"), cond %in% c("r",'cold.ck','drought.ck')) %>%
    #filter(!sid %in% c("s22","s25")) %>%
    filter(!sid %in% c("s22")) %>%
    select(sid,tis,gt,spots,avgLength,note=grp.xlab) %>%
    mutate(gt2 = gt_map[gt]) %>% mutate(gt2=factor(gt2, levels=gt_map)) %>%
    group_by(tis,gt) %>% mutate(rep=1:n()) %>% ungroup() %>%
    mutate(grp = glue("{tis}.{gt2}")) %>% mutate(grp = as_factor(grp)) %>%
    mutate(cond = glue("{tis}.{gt2}.{rep}")) %>% mutate(cond = as_factor(cond)) %>%
    select(sid,tis,gt,gt2,cond,grp,rep,note) %>% print(n=40)
#
shapes = c('single (1 bp)', 'steep  (2-10 bp)', 'broad  (>= 10 bp)')
shapess = c('single', 'steep', 'broad')
#{{{ comparisons
cmps = tibble(cmp = c(
"root.B shoot.B",
"root.A shoot.A",
"root.B root.A",
"shoot.B shoot.A",
"root.H root.B",
"root.H root.A",
"root.H root.C",
"shoot.H shoot.B",
"shoot.H shoot.A",
"shoot.H shoot.C"
)) %>%
    separate(cmp, c("grp1",'grp2'), sep=' ', remove=F) %>%
    mutate(cmp = as_factor(cmp))# %>%
    #left_join(thf %>% select(cond1=cond.s,cond1.l=cond.l), by='cond1') %>%
    #left_join(thf %>% select(cond2=cond.s,cond2.l=cond.l), by='cond2') %>%
    #mutate(cmp.l=glue("{cond1.l} vs {cond2.l}")) %>%
    #mutate(cmp.l = as_factor(cmp.l))
#}}}
cmps5 = cmps$cmp[1:5]
cmps4 = cmps$cmp[c(1,7,5,6)]
cols_shape = pal_npg()(5)
cols_shift = brewer.pal(5, 'Accent')
cols_var = pal_d3()(5)
#}}}
#
read_cage_bws <- function(diri, dsg, minSupport=2, genome='Zmays_B73v5') { # dsg is a df with 'Name'
    #{{{ read in, remove singletons
    fps = glue("{diri}/{dsg$Name}-{genome}.plus.bw")
    fms = glue("{diri}/{dsg$Name}-{genome}.minus.bw")
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
        calcTotalTags() %>%
        calcTPM() %>%
        calcPooled()
    #
    SE_tss = SE_tss_raw %>%
        subsetBySupport(inputAssay='TPM',unexpressed=unexpressed,minSamples=minSamples) %>%
        calcTPM() %>%
        calcPooled() %>%
        calcShape(pooled=SE_ctss, outputColumn='IQR', shapeFunction=shapeIQR,
            lower=.1, upper=.9) %>%
        calcShape(pooled=SE_ctss, outputColumn='entropy', shapeFunction=shapeEntropy) %>%
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


