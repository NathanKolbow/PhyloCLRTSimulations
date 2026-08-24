library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggh4x)
library(cowplot)

############## Rough only first ##############
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df <- rbind(
		read.csv("dat-eps0.0.csv"),
		read.csv("dat-eps0.1.csv"),
		read.csv("dat-eps0.01.csv"),
		read.csv("dat-eps0.001.csv"),
		read.csv("dat-eps0.0001.csv"),
		read.csv("dat-eps0.5.csv"),
		read.csv("dat-eps1.0.csv")
 	) %>%
	mutate(
		ngt = paste0("ngt=", ngt),
		eps = paste0("eps=", eps),
		result = if_else(test == "CLIC", sign(result) * log(abs(result)), result)
	)
nrow(df)


alpha <- 0.05
df %>%
	group_by(test, model, eps, ngt, gamma, t) %>%
	summarise(
		type1rate = mean(if_else(test == "CLIC", result < 0, result <= alpha)),
		result = mean(result),
		.groups = "drop_last"
	) %>%
	mutate(cutoff = if_else(test == "CLIC", 0, alpha)) %>%
	ggplot(aes(x = gamma, y = type1rate, color = ngt)) +
	geom_point() +
	geom_line() +
	facet_nested(test + model ~ eps, scales="free")
