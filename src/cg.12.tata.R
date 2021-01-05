source("functions.R")
require(universalmotif)
dirw = glue("{dird}/12_tata")
setwd(dirw)

fi = 'tata.txt'
mtfs = read_uniprobe(fi)
fo = 'tata.meme'
write_meme(mtfs, fo)
