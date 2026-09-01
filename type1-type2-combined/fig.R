library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggh4x)
library(cowplot)

# WARN: gamma is in all DFs but means different things in type-1 tests vs. type-2 tests

############## Rough only first ##############
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
type1df <- rbind(
		read.csv("../type1error/dat-eps0.0.csv"),
		read.csv("../type1error/dat-eps0.0001.csv"),
		read.csv("../type1error/dat-eps0.001.csv"),
		read.csv("../type1error/dat-eps0.01.csv"),
		read.csv("../type1error/dat-eps0.1.csv"),
		read.csv("../type1error/dat-eps0.5.csv"),
		read.csv("../type1error/dat-eps1.0.csv")
	) %>% mutate(truth = "H0")
type2df <- rbind(
		read.csv("../../CLRT/src/model_expansion/tests-perfect-data/dat-eps0.0.csv"),
		read.csv("../../CLRT/src/model_expansion/tests-perfect-data/dat-eps0.001.csv"),
		read.csv("../../CLRT/src/model_expansion/tests-perfect-data/dat-eps0.01.csv"),
		read.csv("../../CLRT/src/model_expansion/tests-perfect-data/dat-eps0.1.csv"),
		read.csv("../../CLRT/src/model_expansion/tests-perfect-data/dat-eps0.1.csv"),
		read.csv("../../CLRT/src/model_expansion/tests-perfect-data/dat-eps0.5.csv"),
		read.csv("../../CLRT/src/model_expansion/tests-perfect-data/dat-eps1.0.csv")
	) %>% mutate(truth = "H1")
nrow(type1df)
nrow(type2df)


# ROC
# TEST <- "cLR"

make_roc_data <- function(dat) {
	eps <- 1e-12
	dat_alphas <- sort(unique(c(-0.01, 0, dat$result, 1)))
	tibble(
			alpha = dat_alphas,
			fpr = unlist(lapply(dat_alphas, function(x) { mean(dat$result[dat$truth == "H0"] < x+eps) })),
			tpr = unlist(lapply(dat_alphas, function(x) { mean(dat$result[dat$truth == "H1"] < x+eps) }))
			# fpr = unlist(lapply(dat_alphas, function(x) { mean(dat$result[dat$truth == "H1"] > x-eps) })),	# false negative
			# tpr = unlist(lapply(dat_alphas, function(x) { mean(dat$result[dat$truth == "H0"] > x-eps) }))	# true negative
			# fpr = unlist(lapply(dat_alphas, function(x) { mean(dat$result[dat$truth == "H1"] < x+eps) })),	# POWER
			# tpr = unlist(lapply(dat_alphas, function(x) { mean(dat$result[dat$truth == "H0"] > x-eps) }))	# TYPE-1 NON-ERROR
		) %>%
		arrange(alpha) %>%
		uncount(if_else(row_number() == 1, 1, 2)) %>%
		mutate(tpr = lag(tpr, default = tpr[1]))
}

# filtered for now while we have little data
NGAMMA_GROUPS <- 4
roc_type2df <- type2df %>%
	mutate(gamma_group = cut_width(gamma, width = 0.5 / NGAMMA_GROUPS, boundary = 0))
roc_type1df <- type1df %>%
	group_by(model, test, eps, ngt) %>%
	uncount(NGAMMA_GROUPS) %>%
	mutate(gamma_group = levels(roc_type2df$gamma_group)[row_number() %% NGAMMA_GROUPS + 1])
roc_df <- rbind(roc_type2df, roc_type1df) %>%
	group_by(model, test, eps, ngt, gamma_group) %>%
	group_modify(~ make_roc_data(.x)) %>%
	mutate(
		ngtlabel = paste0(ngt, " loci"),
		epslabel = paste0("ϵ = ", eps)
	)
auc_df <- roc_df %>%
	group_by(test, model, eps, ngt, gamma_group, fpr) %>%
	summarise(tpr = max(tpr)) %>%
	mutate(
		fpr1 = lead(fpr, default=1),
		area = (fpr1 - fpr) * tpr
	) %>%
	group_by(test, model, eps, ngt, gamma_group) %>%
	summarise(AUC = sum(area)) %>%
	mutate(
		labely = 0.05 + 0.03 * which(gamma_group == levels(gamma_group)),
		text = sprintf("AUC: %.2f", round(AUC, digits=2))
	) %>%
	drop_na()

# All NGT and eps levels for one test
roc_df %>%
	drop_na() %>%
	ggplot(aes(x = fpr, y = tpr, color = gamma_group)) +
	geom_line(linewidth=1) +
	facet_nested(ngtlabel + eps ~ test + model) +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
	labs(
		x = "False Positive Rate",
		y = "True Positive Rate",
		color = expression(gamma * " under " * H[1])
	) +
	scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
	scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
	theme_bw()


# Labelled example
TEST <- "cLRI"
labeldf <- roc_df %>%
	filter(ngt == 100 & eps == 0.001 & test == TEST & model == "joint") %>%
	drop_na()

tibble_alpha_point <- function(dat, alphas) {
	tbl <- tibble()
	for(a in alphas) {
		mindiff <- min(abs(dat$alpha - a))
		imatch <- which(abs(dat$alpha - a) == mindiff)
		fpr <- max(dat$fpr[imatch])
		tbl <- rbind(tbl,
			tibble(
				alpha = a,
				fpr = fpr
			)
		)
	}
	tbl
}

annotation_df <- labeldf %>%
	group_by(gamma_group) %>%
	filter(gamma_group == levels(gamma_group)[1]) %>%
	group_modify(~ tibble_alpha_point(.x, c(0.01, 0.05))) %>%
	mutate(ymin = 0, ymax = 1)


ggplot(labeldf, aes(x = fpr, y = tpr, color = gamma_group)) +
	geom_line() +
	facet_grid(epslabel ~ ngtlabel) +
	geom_abline(slope = 1, intercept = 0, linetype = "solid", alpha = 1) +
	labs(
		x = "False Positive Rate",
		y = "True Positive Rate",
		color = expression(gamma * " under " * H[1])
	) +
	scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
	scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
	theme_bw() +
	# Specific alpha labelling
	geom_vline(data = annotation_df, aes(xintercept = fpr), linetype="longdash", alpha=0.5) +
	annotate("text", x = 0.11, y = 0.45, label = "α = 0.05", hjust=0) +
	# AUC labelling
	geom_text(data = filter(auc_df, ngt == 100 & eps == 0.001 & test == TEST), aes(x = 0.75, y = labely, color = gamma_group, label = text), size=5)


###### Select N (test, model, eps) combinations for each (ngt, gamma_group) combo ######
S <- auc_df %>%
	group_by(ngt, gamma_group) %>%
	filter(AUC == max(AUC))
Slabeldf <- auc_df %>%
	group_by(ngt, gamma_group) %>%
	filter(AUC == max(AUC)) %>%
	select(test, model, eps, ngt, gamma_group, AUC, text) %>%
	group_by(model, eps, ngt, gamma_group) %>%
	summarise(
		test = if_else(n() == 1, first(test), paste(test, collapse = ", ")),
		AUC = first(AUC),
		text = first(text)
	) %>%
	group_by(ngt, gamma_group) %>%
	filter(row_number() == 1) %>%
	mutate(
		ngtlabel = paste0(ngt, " loci")
	)
best_roc_df <- roc_df %>%
	group_by(ngt, gamma_group) %>%
	semi_join(S, by = c("test", "model", "eps", "ngt", "gamma_group")) %>%
	group_by(ngt, gamma_group) %>%
	arrange(desc(eps)) %>%
	mutate(
		test = if_else(length(unique(test)) == 1, paste0(first(test), "-", paste(unique(model), collapse=","), " (eps = ", paste(unique(eps), collapse=","), ")"), "Tie")
	) %>%
	ungroup() %>%
	arrange(tpr) %>%
	distinct(model, eps, ngt, gamma_group, test, fpr, .keep_all=TRUE)
prepend_df <- best_roc_df %>%
	group_by(ngt, gamma_group) %>%
	filter(row_number() == 1) %>%
	mutate(fpr = 0, tpr = 0)
append_df <- best_roc_df %>%
	filter(row_number() == n()) %>%
	mutate(fpr = 1)

proc_best <- rbind(prepend_df, best_roc_df, append_df) %>%
	ggplot(aes(x = fpr, y = tpr, color = test, fill = test)) +
	geom_line(linewidth=1) +
	geom_area(position = "identity", alpha = 0.25) +
	facet_grid(ngtlabel ~ gamma_group) +
	scale_x_continuous(expand = c(0.0125, 0.0125), limits = c(0, 1)) +
	scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
	theme_bw() +
	geom_text(data = Slabeldf, aes(x = 0.5, y = 0.5, label = text), fontface = "bold", color = "black", inherit.aes=FALSE, show.legend=FALSE, hjust=0.5) +
	labs(
		x = "False Positive Rate",
		y = "True Positive Rate",
		color = "Test",
		fill = "Test"
	)
proc_best

pdf("roc-best_tests.pdf", width=10, height=7)
proc_best
dev.off()


###### SAME AS ABOVE, BUT: show all ties individually ######
S <- auc_df %>%
	group_by(ngt, gamma_group) %>%
	filter(AUC == max(AUC))
Slabeldf <- auc_df %>%
	group_by(ngt, gamma_group) %>%
	filter(AUC == max(AUC)) %>%
	select(test, model, eps, ngt, gamma_group, AUC, text) %>%
	group_by(model, eps, ngt, gamma_group) %>%
	summarise(
		test = if_else(n() == 1, first(test), paste(test, collapse = ", ")),
		AUC = first(AUC),
		text = first(text)
	) %>%
	group_by(ngt, gamma_group) %>%
	filter(row_number() == 1) %>%
	mutate(
		ngtlabel = paste0(ngt, " loci")
	)
best_roc_df <- roc_df %>%
	group_by(ngt, gamma_group) %>%
	semi_join(S, by = c("test", "model", "eps", "ngt", "gamma_group")) %>%
	group_by(ngt, gamma_group) %>%
	arrange(desc(eps)) %>%
	mutate(
		test = if_else(length(unique(test)) == 1, paste0(first(test), "-", paste(unique(model), collapse=","), " (ϵ = ", paste(unique(eps), collapse=","), ")"), paste0("Tie: ", paste(unique(test), collapse=",")))
	) %>%
	ungroup() %>%
	arrange(tpr) %>%
	distinct(model, eps, ngt, gamma_group, test, fpr, .keep_all=TRUE)
prepend_df <- best_roc_df %>%
	group_by(ngt, gamma_group) %>%
	filter(row_number() == 1) %>%
	mutate(fpr = 0, tpr = 0)
append_df <- best_roc_df %>%
	filter(row_number() == n()) %>%
	mutate(fpr = 1)

rbind(prepend_df, best_roc_df, append_df) %>%
	ggplot(aes(x = fpr, y = tpr, color = test, fill = test)) +
	geom_line(linewidth=1) +
	geom_area(position = "identity", alpha = 0.25) +
	facet_grid(ngtlabel ~ gamma_group) +
	scale_x_continuous(expand = c(0.0125, 0.0125), limits = c(0, 1)) +
	scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
	theme_bw() +
	geom_text(data = Slabeldf, aes(x = 0.5, y = 0.5, label = text), fontface = "bold", color = "black", inherit.aes=FALSE, show.legend=FALSE, hjust=0.5) +
	labs(
		x = "False Positive Rate",
		y = "True Positive Rate",
		color = "Test",
		fill = "Test"
	)


###### Select best (model, eps) for each (ngt, gamma_group) combo for ONE TEST ######
TEST <- "cw"
S <- auc_df %>%
	filter(test == TEST) %>%
	group_by(ngt, gamma_group) %>%
	filter(AUC == max(AUC))
Slabeldf <- auc_df %>%
	filter(test == TEST) %>%
	group_by(ngt, gamma_group) %>%
	filter(AUC == max(AUC)) %>%
	select(test, model, eps, ngt, gamma_group, AUC, text) %>%
	group_by(model, eps, ngt, gamma_group) %>%
	summarise(
		test = if_else(n() == 1, first(test), paste(test, collapse = ", ")),
		AUC = first(AUC),
		text = first(text)
	) %>%
	group_by(ngt, gamma_group) %>%
	filter(row_number() == 1) %>%
	mutate(
		ngtlabel = paste0(ngt, " loci")
	)
best_roc_df <- roc_df %>%
	filter(test == TEST) %>%
	group_by(ngt, gamma_group) %>%
	semi_join(S, by = c("test", "model", "eps", "ngt", "gamma_group")) %>%
	group_by(ngt, gamma_group) %>%
	arrange(desc(eps)) %>%
	mutate(
		test = if_else(length(unique(test)) == 1, paste0(first(test), "-", paste(unique(model), collapse=","), " (ϵ = ", paste(unique(eps), collapse=","), ")"), "Tie")
	) %>%
	ungroup() %>%
	arrange(tpr) %>%
	distinct(model, eps, ngt, gamma_group, test, fpr, .keep_all=TRUE)
prepend_df <- best_roc_df %>%
	group_by(ngt, gamma_group) %>%
	filter(row_number() == 1) %>%
	mutate(fpr = 0, tpr = 0)
append_df <- best_roc_df %>%
	filter(row_number() == n()) %>%
	mutate(fpr = 1)

proc_bestcw <- rbind(best_roc_df, append_df) %>%
	ggplot(aes(x = fpr, y = tpr, color = test, fill = test)) +
	geom_line(linewidth=1) +
	geom_area(position = "identity", alpha = 0.25) +
	facet_grid(ngtlabel ~ gamma_group) +
	scale_x_continuous(expand = c(0.0125, 0.0125), limits = c(0, 1)) +
	scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
	theme_bw() +
	geom_text(data = Slabeldf, aes(x = 0.5, y = 0.5, label = text), fontface = "bold", color = "black", inherit.aes=FALSE, show.legend=FALSE, hjust=0.5) +
	labs(
		x = "False Positive Rate",
		y = "True Positive Rate",
		color = "Test",
		fill = "Test"
	)
proc_bestcw

pdf("roc-best_cw.pdf", width=10, height=7)
proc_bestcw
dev.off()
