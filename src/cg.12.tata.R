source("functions.R")
require(universalmotif)
dirw = glue("{dird}/12_tata")
setwd(dirw)

fi = 'tata.txt'
y = read_uniprobe(fi)
