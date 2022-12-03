rm(list=ls())
library(tidyverse)
library(stringi)
library(ggplot2)
Sys.setlocale("LC_TIME", "English")

###################################################################################################
# FTM surveillance
###################################################################################################

setwd("...")
load("calendar.RData")

# maf <- data.frame()
# for (w in calendar$colWeek) {
#   
#   load(paste0("w.muta.", w, ".RData")) # load weekly mutation data
#   load(paste0("w.meta.", w, ".RData")) # load weekly metadata
#   
#   w.muta.modify <- w.muta %>%
#     mutate(mutation = paste0(refvar, refpos, qvar)) %>%
#     select(sample, mutation)
#   
#   occur.prev <- data.frame(table(w.muta.modify$mutation))
#   names(occur.prev) <- c("mutation", "count")
#   occur.prev$n <- nrow(w.meta)
#   occur.prev$prev <- with(occur.prev, count / n * 100)
#   
#   maf <- rbind(maf, data.frame(colWeek = w,
#                                occur.prev))
#   rownames(maf) <- NULL
#   # save(occur.prev, file = paste0("./MAF/maf.", w, ".RData"))
#   
# }

# save(maf, file = paste0("/maf.all.RData"))
load("maf.all.RData")

muta.list <- unique(as.character(maf %>% filter(prev > 0.5) %>% pull(mutation))) #Only show mutations with MAF > 0.5%

dummy <- calendar %>%
  select(colWeek, n)
maf.all <- data.frame()
for (m in muta.list) {
  x <- subset(maf, mutation==m)
  x <- merge(dummy, x, by = c("colWeek", "n"), all.x = T) %>%
    mutate(mutation = m,
           count = ifelse(is.na(count), 0, count),
           prev = ifelse(is.na(prev), 0, prev))
  maf.all <- rbind(maf.all, x)
}

Sys.setlocale("LC_TIME", "English")
ts.label <- rep("", nrow(calendar))
ts.label[seq(1, nrow(calendar), 1)] <- format(calendar$Start[seq(1, nrow(calendar), 1)], "%b %d, %y")
# pdf("MAF.pdf", width = 14, height = 7)
ggplot(data = maf.all, aes(x = colWeek, y = prev, group = mutation)) +
  geom_line() +
  scale_x_discrete(labels = ts.label) +
  labs(x = "Sampleing week", y = "Mutation frequency (%)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid = element_line(colour = "grey88"),
        legend.key = element_rect(fill = NA),
        # axis.title = element_text(size = 40),
        # axis.text = element_text(size = 30),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
# dev.off()
