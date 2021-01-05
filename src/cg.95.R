source('functions.R')
require(magick)
require(ggplotify)
setwd(dirf)

#{{{ Fig 1
a = as.grob(image_read_pdf("regions.pdf"))
b = readRDS('f1b.rds')
c = readRDS('f1c.rds') + o_margin(0,.3,0, 1.5) +
    theme(legend.position=c(.5,1), legend.justification=c(.5,1), legend.direction='vertical')
d = readRDS('f1d.rds') +
    o_margin(.3,.3,0,1)
ggarrange(a, c,
    ggarrange(b, d, nrow=1, ncol=2, widths=c(1,1,1), labels=LETTERS[3:4]),
    nrow=3, heights=c(1,1.5,1.5), labels = LETTERS[1:2]) %>%
ggexport(filename="f1.pdf", width=6, height=8)
#}}}

#{{{ f2
a = readRDS('f2a.rds')
b = readRDS('f2b.rds')
c = readRDS('f2c.rds')
ggarrange(a, b, c,
    nrow=1, ncol=3, widths=c(1,1,1.2), labels = LETTERS[1:3]) %>%
ggexport(filename="f2.pdf", width=9, height=3)
#}}}

#{{{ sf01
a = readRDS('sf02a.rds') + facet_wrap(~pnl,nrow=4,scale='free')
b = readRDS('sf02b.rds') + facet_wrap(~pnl,nrow=4) +
    theme(legend.position='none')
ggarrange(a, b, nrow=1, widths=c(1,1), heights=c(1,1),
      labels = LETTERS[1:2]) %>%
ggexport(filename="sf02.pdf", width=5, height=8)
#}}}

