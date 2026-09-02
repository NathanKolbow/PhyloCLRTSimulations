setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(ggplot2)

df <- read.csv("dat-eps0.0.csv") %>%
	mutate(result = if_else(test == "CLIC", sign(result) * log10(abs(result) + 0.1), result))
df %>%
	mutate(ngt = factor(ngt)) %>%
	group_by(model, test, eps, gamma, ngt) %>%
	summarise(result = mean(result), .groups="drop") %>%
	ggplot(aes(x = gamma, y = result, color = ngt)) +
	geom_hline(yintercept=0, color="black", linetype="dashed", alpha=0.25) +
	geom_line() +
	facet_grid(test ~ model, scales="free")

rbind(
	filter(df, model == "quartet") %>%
		rename(terr=snaqabssumterr, gammaerr=snaqabssumgammaerr) %>%
		select(ngt, model, terr, gammaerr),
	filter(df, model != "quartet") %>%
		rename(terr=newabssumterr, gammaerr=newabssumgammaerr) %>%
		select(ngt, model, terr, gammaerr)
) %>%
	distinct() %>%
	pivot_longer(cols = c(terr, gammaerr), values_to = "value") %>%
	rename(errortype = name) %>%
	group_by(ngt, model, errortype) %>%
	summarise(
		ymean = mean(value),
		ymin = min(value),
		ymax = max(value),
		.groups = "drop"
	) %>%
	mutate(errortype = if_else(errortype == "terr", "Edge", "Gamma")) %>%
	mutate(ngt = factor(ngt)) %>%
	ggplot(aes(x = ngt, color = model)) +
	geom_errorbar(aes(y = ymean, ymin = ymin, ymax = ymax), position=position_dodge2()) +
	geom_point(aes(y = ymean), position=position_dodge2(0.9)) +
	facet_wrap(~errortype, scales="free")
