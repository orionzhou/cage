source("functions.R")
dirw = glue("{dird}/09_igv")
setwd(dirw)
source("~/source/git/trackplot/trackplot.R")

fis = c("B73_shoot.plus.bw","B73_shoot.minus.bw")
fis = glue("{dird}/tracks/merged_bw/{fis}")

gtf = "~/projects/genome/data/Zmays_B73/50_annotation/15.gtf"
mark = tibble(chr = 'B10',
              start=90664308,
              end=90665308,
              name='hs1'
              )
p = trackplot(
  bigWigs = fis,
  loci = "B10:90,663,308-90,665,449",
  draw_gene_track = T,
  gene_model = gtf, isGTF = TRUE,
  mark_regions = mark,
  custom_names = c("shoot +", "shoot -")
)
fo = glue("{dirw}/t.pdf")
ggsave(p, fo, width=8, height=6)
