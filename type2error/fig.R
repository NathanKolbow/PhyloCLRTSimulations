setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(ggplot2)

df <- rbind(
	read.csv("dat-eps0.0.csv"),
	read.csv("dat-eps0.1.csv"),
	read.csv("dat-eps0.01.csv"),
	read.csv("dat-eps0.001.csv"),
	read.csv("dat-eps0.0001.csv"),
	read.csv("dat-eps1.0.csv")
) %>%
	mutate(
		# sign-preserving, monotone log scale (log10(|x| + 0.1) flips the sign of |x| < 0.9)
		result = if_else(test == "CLIC", sign(result) * log10(1 + abs(result)), result),
		ngt = factor(ngt)
	)

# Average results
plots <- list()
for(itertest in unique(df$test)) {
	plots[[length(plots)+1]] <- df %>%
		filter(test == itertest) %>%
		group_by(model, test, eps, gamma, ngt) %>%
		summarise(result = mean(result), .groups="drop") %>%
		ggplot(aes(x = gamma, y = result, color = ngt)) +
		geom_hline(yintercept=0, color="black", linetype="dashed", alpha=0.25) +
		geom_line() +
		facet_grid(model ~ eps, scales="free") +
		ggtitle(itertest) +
		theme_bw()
}
plots


# Parameter estimates
# NOTE: `newabssumterr` in data generated before `absparamerrors` is fixed is inflated by
# pendant edges (-1 vs. 1) and root placement; see bug-scripts/04-absparamerrors.jl
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
