#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
project <- normalizePath(getwd())
.libPaths(c(file.path(project, "r_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggrepel)
  library(svglite)
  library(ragg)
  library(grid)
})

src <- file.path(project, "figures", "source_data")
main_dir <- file.path(project, "output")
ext_dir <- main_dir
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(
  ink = "#272727", grey = "#767676", light = "#D8D8D8", pale = "#F3F4F6",
  blue = "#2F6DA4", blue2 = "#83A9CE", teal = "#319A91", teal2 = "#A8D8D1",
  red = "#C94E48", red2 = "#E9AAA6", orange = "#D98B35", violet = "#8064A2",
  green = "#4C956C", gold = "#D6A33A"
)

theme_nc <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = pal[["ink"]]),
      axis.ticks = element_line(linewidth = 0.30, colour = pal[["ink"]]),
      axis.text = element_text(colour = pal[["ink"]], size = base_size - 0.4),
      axis.title = element_text(colour = pal[["ink"]], size = base_size),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.6),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = pal[["ink"]]),
      plot.title = element_text(size = base_size + 0.8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = pal[["grey"]]),
      plot.tag = element_text(size = 8, face = "bold", colour = pal[["ink"]]),
      panel.grid = element_blank(),
      legend.key.height = unit(3.5, "mm"),
      plot.margin = margin(3, 4, 3, 3)
    )
}
theme_set(theme_nc())

read_src <- function(fig, name) read.delim(file.path(src, fig, paste0(name, ".tsv")), check.names = FALSE)
clean_fate <- function(x) {
  x <- sub("^p_", "", x)
  x <- sub("^external_", "", x)
  gsub("_", " ", x)
}
clean_model <- function(x) {
  recode(x,
    baseline_latent_ridge = "baseline latent",
    baseline_scvi_ridge = "baseline latent",
    baseline_plus_perturbation = "baseline + perturbation",
    perturbation_descriptor = "perturbation descriptor",
    biological_program_ridge = "biological programs",
    biological_axis_ridge = "biological axes",
    training_patient_mean = "training mean",
    training_mean = "training mean",
    .default = gsub("_", " ", x)
  )
}
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}
scale_y_reordered <- function(...) {
  scale_y_discrete(labels = function(x) gsub("___.*$", "", x), ...)
}

save_pub <- function(plot, filename, width_mm=183, height_mm=145, dpi=300) {
  aliases <- c(Figure1_study_design="Figure1 study design",Figure2_lineage_prediction="Figure2 lineage prediction",Figure3_generalization="Figure3 generalization",Figure4_mechanism_axes="Figure4 mechanism axes",Figure5_spatial_niches="Figure5 spatial niches",Figure6_patient_fate_embedding="Figure6 patient fate embedding")
  base <- basename(filename)
  selection <- strsplit(Sys.getenv("FIGURE_SELECTION", ""), ",", fixed = TRUE)[[1]]
  if (length(selection) && nzchar(selection[1]) && !base %in% selection) return(invisible(NULL))
  if(base %in% names(aliases)) filename <- file.path(dirname(filename), aliases[[base]])
  if (base %in% c("Figure4_mechanism_axes", "Figure5_spatial_niches", "Figure6_patient_fate_embedding") && capabilities("aqua")) {
    grDevices::quartz(type="pdf",file=paste0(filename,".pdf"),width=width_mm/25.4,height=height_mm/25.4,family="Helvetica")
  } else if (base %in% c("Figure4_mechanism_axes", "Figure5_spatial_niches", "Figure6_patient_fate_embedding") && capabilities("cairo")) {
    grDevices::cairo_pdf(paste0(filename,".pdf"),width=width_mm/25.4,height=height_mm/25.4,family="Helvetica")
  } else {
    grDevices::pdf(paste0(filename,".pdf"),width=width_mm/25.4,height=height_mm/25.4,family="Helvetica",useDingbats=FALSE)
  }
  print(plot);dev.off()
  ragg::agg_png(file.path(main_dir,paste0(basename(filename),".png")),width=width_mm/25.4,height=height_mm/25.4,units="in",res=200,background="white")
  print(plot);dev.off()
}

mean_ci <- function(data, groups, value = "value") {
  data %>% group_by(across(all_of(groups))) %>% summarise(
    mean = mean(.data[[value]], na.rm = TRUE),
    se = sd(.data[[value]], na.rm = TRUE) / sqrt(sum(is.finite(.data[[value]]))),
    n = sum(is.finite(.data[[value]])), .groups = "drop"
  ) %>% mutate(ci = mapply(function(nn, ss) if (nn > 1) qt(0.975, nn - 1) * ss else NA_real_, n, se))
}

tile_text_plot <- function(data, x, y, fill, label = NULL, low = pal[["red2"]], mid = "white", high = pal[["blue"]]) {
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_tile(colour = "white", linewidth = 0.45) +
    scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = 0, na.value = pal[["light"]])
  if (!is.null(label)) p <- p + geom_text(aes(label = .data[[label]]), size = 1.65)
  p + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.title = element_blank())
}

# Figure 1 -------------------------------------------------------------------
inv <- read_src("figure1", "dataset_scale_and_evidence")
inv$evidence_tier <- gsub("gold lineage supervision", "lineage-linked", inv$evidence_tier)
inv$evidence_tier <- gsub("silver transport supervision", "transport-derived", inv$evidence_tier)
inv$qc_or_target_numeric[grepl("BeatAML", inv$dataset)] <- 317
ood <- read_src("figure1", "ood_validation_contract")

inv_plot <- inv %>% filter(is.finite(input_numeric), !grepl("post-treatment", dataset)) %>%
  mutate(dataset = factor(dataset, levels = rev(dataset)), evidence_tier = gsub(" probability", "", evidence_tier))
p1b <- ggplot(inv_plot, aes(input_numeric, dataset, colour = evidence_tier)) +
  geom_segment(aes(x = pmax(qc_or_target_numeric, 1), xend = input_numeric, yend = dataset), colour = pal[["light"]], linewidth = 1.2) +
  geom_point(aes(x = pmax(qc_or_target_numeric, 1)), shape = 21, fill = "white", size = 2.1, stroke = .65) +
  geom_point(size = 1.9) + scale_x_log10(labels = label_number(scale_cut = cut_si(""))) +
  scale_colour_manual(values = c("lineage-linked" = pal[["teal"]], "transport-derived" = pal[["violet"]],
                                  "derived program" = pal[["blue"]], "genetic derived program" = pal[["blue2"]],
                                  "spatial exploratory" = pal[["orange"]], "clinical projection" = pal[["red"]])) +
  labs(x = "Input and analysed units (log scale)", y = NULL, colour = "Evidence tier") + theme(legend.position = "bottom")

ood <- ood %>% mutate(validation = factor(validation, levels = c("leave-cell-line", "leave-perturbation", "leave-task", "leave-dataset")))
p1c <- ggplot(ood, aes(x = validation, y = 1, fill = validation)) +
  geom_tile(width = .88, height = .52, colour = "white", linewidth = .5) +
  geom_text(aes(label = held_out_unit), size = 1.8, lineheight = .9) +
  scale_fill_manual(values = c(pal[["blue2"]], pal[["teal2"]], "#DDD2EA", pal[["red2"]])) +
  scale_y_continuous(NULL, breaks = NULL) + labs(x = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(face = "bold", angle = 22, hjust = 1, size = 6), legend.position = "none")

tier_counts <- inv_plot %>% count(evidence_tier)
p1d <- ggplot(tier_counts, aes(n, reorder(evidence_tier, n), fill = evidence_tier)) +
  geom_col(width = .62) + geom_text(aes(label = n), hjust = -0.4, size = 2) +
  scale_fill_manual(values = c("lineage-linked" = pal[["teal"]], "transport-derived" = pal[["violet"]],
                               "derived program" = pal[["blue"]], "genetic derived program" = pal[["blue2"]],
                               "spatial exploratory" = pal[["orange"]], "clinical projection" = pal[["red"]])) +
  scale_x_continuous(expand = expansion(mult = c(0, .18)), breaks = 0:4) + labs(x = "Data sources", y = NULL) + theme(legend.position = "none")

p1b <- p1b + theme(legend.position = "none")
fig1 <- p1b + p1c + p1d + plot_layout(design = "AB
AC", widths = c(1, 1.08), heights = c(.8, 1)) +
  plot_annotation(tag_levels = "a")
save_pub(fig1, file.path(main_dir, "Figure1_study_design"), 183, 118)

# Figure 2 -------------------------------------------------------------------
branch <- read_src("figure2", "larry_branch_fold_metrics") %>% mutate(model = clean_model(model))
scalar <- read_src("figure2", "larry_scalar_fold_metrics") %>% mutate(model = clean_model(model), target = clean_fate(target))
lineage_scale <- read_src("figure2", "larry_lineage_scale")

branch_macro <- branch %>% filter(branch == "macro", metric %in% c("brier", "top_branch_accuracy"))
branch_sum <- mean_ci(branch_macro, c("experiment", "model", "metric"))
p2b <- ggplot(branch_sum %>% filter(metric == "brier"), aes(model, mean, fill = model)) +
  geom_col(width = .62) + geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = .15, linewidth = .35) +
  geom_point(data = branch_macro %>% filter(metric == "brier"), aes(y = value), position = position_jitter(width = .08), size = .65, inherit.aes = TRUE) +
  facet_wrap(~experiment) + scale_fill_manual(values = c("baseline latent" = pal[["blue"]], "training mean" = pal[["light"]])) +
  labs(x = NULL, y = "Macro Brier score") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
p2c <- ggplot(branch_sum %>% filter(metric == "top_branch_accuracy"), aes(model, mean, fill = model)) +
  geom_col(width = .62) + geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), width = .15, linewidth = .35) +
  geom_point(data = branch_macro %>% filter(metric == "top_branch_accuracy"), aes(y = value), position = position_jitter(width = .08), size = .65) +
  facet_wrap(~experiment) + scale_fill_manual(values = c("baseline latent" = pal[["blue"]], "training mean" = pal[["light"]])) +
  labs(x = NULL, y = "Top-branch accuracy") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

scalar_corr <- scalar %>% filter(metric == "spearman", model == "baseline latent")
p2d <- ggplot(scalar_corr, aes(value, target, colour = experiment)) +
  geom_vline(xintercept = 0, linetype = 2, colour = pal[["grey"]], linewidth = .3) +
  geom_point(alpha = .55, size = .75, position = position_jitter(height = .09)) +
  stat_summary(fun = mean, geom = "point", size = 2.0, shape = 18, position = position_dodge(width = .25)) +
  scale_colour_manual(values = c("in vitro" = pal[["blue"]], cytokine = pal[["teal"]])) +
  labs(x = "Held-clone Spearman correlation", y = NULL, colour = NULL) + theme(legend.position = "bottom")

scalar_mae <- scalar %>% filter(metric == "mae") %>% select(experiment, fold, target, model, value) %>%
  pivot_wider(names_from = model, values_from = value) %>% mutate(improvement = `training mean` - `baseline latent`)
p2e <- ggplot(scalar_mae, aes(improvement, target, colour = experiment)) +
  geom_vline(xintercept = 0, linetype = 2, colour = pal[["grey"]], linewidth = .3) +
  geom_point(alpha = .55, size = .75, position = position_jitter(height = .09)) +
  stat_summary(fun = mean, geom = "point", size = 2.0, shape = 18, position = position_dodge(width = .25)) +
  scale_colour_manual(values = c("in vitro" = pal[["blue"]], cytokine = pal[["teal"]])) +
  labs(x = "MAE improvement over training mean", y = NULL) + guides(colour = "none") + theme(legend.position = "none")

fig2 <- (p2b | p2c) / (p2d | p2e) + plot_layout(heights = c(1, 1.08)) +
  plot_annotation(tag_levels = "a")
save_pub(fig2, file.path(main_dir, "Figure2_lineage_prediction"), 183, 126)

# Figure 3 -------------------------------------------------------------------
mix <- read_src("figure3", "mixseq_generalization_summary") %>% mutate(
  fate = clean_fate(target), validation = gsub("_", "-", validation), model_label = clean_model(model)
)
mix_heat <- mix %>% filter(model == "baseline_scvi_ridge", metric %in% c("brier", "spearman")) %>%
  mutate(display = ifelse(metric == "brier", improvement_over_training_mean, estimate),
         metric_label = ifelse(metric == "brier", "Brier improvement", "Spearman"), label = sprintf("%.3f", display))
p3a <- ggplot(mix_heat, aes(fate, interaction(validation, metric_label, sep = " | "), fill = display)) +
  geom_tile(colour = "white", linewidth = .5) + geom_text(aes(label = label), size = 1.65) +
  scale_fill_gradient2(low = pal[["red"]], mid = "white", high = pal[["blue"]], midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Performance") + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

sens <- read_src("figure3", "mixseq_target_sensitivity") %>% filter(metric == "brier") %>%
  select(variant, validation, target, model, mean_value) %>% pivot_wider(names_from = model, values_from = mean_value) %>%
  mutate(improvement = training_mean - baseline_scvi_ridge, fate = clean_fate(target), validation = gsub("_", "-", validation))
sens_summary <- sens %>% group_by(validation, fate) %>% summarise(wins = sum(improvement > 0), settings = n(), .groups = "drop") %>%
  mutate(label = paste0(wins, "/", settings), fraction = wins / settings)
p3b <- ggplot(sens_summary, aes(fate, validation, fill = fraction)) + geom_tile(colour = "white", linewidth = .5) +
  geom_text(aes(label = label), size = 1.8) + scale_fill_gradient(low = "white", high = pal[["blue"]], limits = c(0,1)) +
  labs(x = NULL, y = NULL, fill = "Brier wins\nacross settings") + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

ablation_plot <- function(label, colour) {
  d <- read_src("figure3", paste0("replogle_", label, "_ablation")) %>%
    filter(metric %in% c("brier", "spearman")) %>% mutate(fate = clean_fate(target), validation = gsub("_", "-", validation), metric = tools::toTitleCase(metric))
  ggplot(d, aes(baseline_state_increment, fate, colour = validation)) +
    geom_vline(xintercept = 0, linetype = 2, colour = pal[["grey"]], linewidth = .3) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = .18, linewidth = .4) + geom_point(size = 1.4) +
    facet_wrap(~metric, scales = "free_x") + scale_x_continuous(breaks=function(x) pretty(x,n=2),labels=function(x) format(x,scientific=TRUE,digits=1)) + scale_colour_manual(values = c("leave-perturbation" = colour, "leave-gem-group" = pal[["red"]])) +
    labs(x = "Baseline-state increment", y = NULL, colour = NULL, title = toupper(label)) + theme(legend.position = "bottom",axis.text.x=element_text(size=5.6,angle=20,hjust=1))
}
p3c <- ablation_plot("k562", pal[["blue"]])
p3d <- ablation_plot("rpe1", pal[["teal"]])
p3d <- p3d + guides(colour = "none")

held <- bind_rows(
  read_src("figure3", "replogle_k562_ablation") %>% mutate(cell_line = "K562"),
  read_src("figure3", "replogle_rpe1_ablation") %>% mutate(cell_line = "RPE1")
) %>% filter(validation == "leave_gem_group", metric == "spearman") %>% mutate(fate = clean_fate(target))
p3e <- ggplot(held, aes(fate, baseline_state_increment, fill = cell_line)) + geom_hline(yintercept = 0, linetype = 2, colour = pal[["grey"]]) +
  geom_col(position = position_dodge(width = .72), width = .62) +
  geom_errorbar(aes(ymin=ci_low,ymax=ci_high),position=position_dodge(width=.72),width=.15,linewidth=.35) +
  scale_fill_manual(values = c(K562 = pal[["blue2"]], RPE1 = pal[["teal2"]])) +
  labs(x = NULL, y = "Held-GEM correlation increment", fill = NULL) + theme(legend.position = "bottom")

fig3 <- (p3a | p3b) / (p3c | p3d) / p3e + plot_layout(heights = c(.9, 1.05, .82)) +
  plot_annotation(tag_levels = "a")
save_pub(fig3, file.path(main_dir, "Figure3_generalization"), 183, 155)

# Figure 4 -------------------------------------------------------------------
figure4_theme <- theme_nc(7.5) + theme(
  plot.title = element_text(size = 8.2, face = "bold", margin = margin(b = 5)),
  axis.title = element_text(size = 7.5), axis.text = element_text(size = 7),
  legend.text = element_text(size = 6.8), legend.key.width = unit(4, "mm"),
  plot.margin = margin(4, 5, 3, 3)
)

wot <- read_src("figure4", "wot_replicate_metrics") %>% mutate(model_label = clean_model(model))
wot_macro <- wot %>% filter(branch == "macro", metric == "brier") %>%
  mutate(model_label = factor(model_label, levels = c("training mean", "biological programs", "baseline latent")))
wot_summary <- wot_macro %>% group_by(model_label) %>%
  summarise(value = mean(value), low = min(value), high = max(value), .groups = "drop")
p4a <- ggplot(wot_macro, aes(value, model_label, colour = model_label)) +
  geom_segment(data = wot_summary, aes(x = low, xend = high, yend = model_label), linewidth = .65, alpha = .55) +
  geom_point(size = 1.5, alpha = .6, position = position_jitter(height = .10, width = 0, seed = 401)) +
  geom_point(data = wot_summary, shape = 18, size = 2.5) +
  scale_colour_manual(values = c("baseline latent" = pal[["blue"]], "biological programs" = pal[["teal"]], "training mean" = pal[["grey"]])) +
  scale_y_discrete(labels = c("baseline latent" = "Baseline state", "biological programs" = "Biological programs", "training mean" = "Training mean")) +
  scale_x_log10(breaks = c(5e-6, 1e-5, 2e-5, 5e-5), labels = c("5", "10", "20", "50"), expand = expansion(mult = c(.10, .10))) +
  labs(x = expression("Macro Brier score ("*10^{-6}*"; log scale)"), y = NULL, title = "Transport-derived branch targets") +
  figure4_theme + theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank())

cross_raw <- read_src("figure4", "cross_task_dataset_metrics") %>%
  filter(metric == "spearman", model != "training_mean") %>%
  mutate(transfer = factor(validation, levels = c("leave_dataset_cross_platform_stress", "leave_dataset_larry_primary", "leave_task")),
         outcome = recode(target, relative_entropy = "Entropy", relative_rare_fate = "Rare fate"),
         y = as.numeric(transfer) + ifelse(outcome == "Entropy", .13, -.13))
cross <- cross_raw %>% group_by(transfer, outcome, y) %>% summarise(value = mean(value), .groups = "drop")
p4b <- ggplot(cross_raw, aes(value, y, colour = outcome)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#B5BAC0", linewidth = .35) +
  geom_point(size = 1.1, alpha = .38, position = position_jitter(height = .035, width = 0, seed = 402)) +
  geom_point(data = cross, aes(shape = outcome), size = 2.2) +
  scale_shape_manual(values = c(Entropy = 16, `Rare fate` = 17)) +
  scale_colour_manual(values = c(Entropy = pal[["blue"]], `Rare fate` = pal[["teal"]])) +
  scale_y_continuous(breaks = 1:3, labels = c("Cross-platform", "LARRY dataset", "Cytokine task"), limits = c(.55, 3.45), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-.1, .2, .1), limits = c(-.17, .25)) +
  labs(x = "Held-group Spearman correlation", y = NULL, colour = NULL, shape = NULL, title = "Cross-context transfer") +
  figure4_theme + theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())

axes <- read_src("figure4", "top_stable_biological_axes") %>%
  filter(target %in% c("p_apoptosis", "p_proliferation", "p_stress")) %>%
  group_by(target) %>% slice_max(abs(mean_coefficient), n = 5, with_ties = FALSE) %>% ungroup() %>%
  mutate(feature = sub("^[^_]+__", "", feature), feature = recode(feature, HEME_METABOLISM = "Heme metabolism"),
         outcome = factor(target, levels = c("p_apoptosis", "p_proliferation", "p_stress"), labels = c("Apoptosis", "Proliferation", "Stress")),
         feature_order = reorder_within(feature, mean_coefficient, outcome),
         coefficient = mean_coefficient * 1000, coefficient_sd = sd_coefficient * 1000,
         direction = ifelse(mean_coefficient > 0, "Positive", "Negative"))
coef_limit <- ceiling(max(abs(axes$coefficient) + axes$coefficient_sd) / 5) * 5
p4c <- ggplot(axes, aes(coefficient, feature_order, colour = direction)) +
  geom_vline(xintercept = 0, colour = "#C9CDD1", linewidth = .4) +
  geom_segment(aes(x = 0, xend = coefficient, yend = feature_order), linewidth = .5, alpha = .35) +
  geom_errorbar(aes(xmin = coefficient - coefficient_sd, xmax = coefficient + coefficient_sd), orientation = "y", width = .15, linewidth = .45) +
  geom_point(size = 2) + facet_wrap(~outcome, scales = "free_y", nrow = 1) + scale_y_reordered() +
  scale_colour_manual(values = c(Positive = pal[["blue"]], Negative = pal[["red"]]), breaks = c("Positive", "Negative")) +
  guides(colour = guide_legend(override.aes = list(shape = 16, linewidth = 0, alpha = 1))) +
  scale_x_continuous(limits = c(-coef_limit, coef_limit), breaks = c(-10, 0, 10), expand = expansion(mult = .06)) +
  labs(x = expression("Conditional coefficient ("*10^{-3}*")"), y = NULL, colour = NULL, title = "Outcome-associated transcriptional features") +
  figure4_theme + theme(legend.position = "bottom", legend.key.height = unit(2.5, "mm"), legend.margin = margin(t = 0, b = 0), strip.text = element_text(size = 7.5, face = "bold"),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.spacing.x = unit(7, "mm"))

null_dir <- file.path(project, "results")
null_raw <- read.delim(file.path(null_dir, "axis_permutation.tsv")) %>% filter(replicate >= 0)
stopifnot(all(table(null_raw$target) == 499))
null <- read.delim(file.path(null_dir, "axis_permutation_summary.tsv"))
fate_labels <- c(p_apoptosis = "Apoptosis", p_proliferation = "Proliferation", p_stress = "Stress", external_resistance_quantile = "Resistance")
fate_levels <- rev(names(fate_labels))
null_raw$fate <- factor(null_raw$target, levels = fate_levels, labels = fate_labels[fate_levels])
null <- null %>% mutate(fate = factor(target, levels = fate_levels, labels = fate_labels[fate_levels]), q_label = sprintf("q = %.2f", qvalue))
p4d <- ggplot(null_raw, aes(stable_axes, fate)) +
  geom_violin(orientation = "y", fill = "#E5E8EC", colour = "#C5CAD1", linewidth = .3, width = .72, scale = "width", trim = TRUE) +
  geom_errorbar(data = null, aes(x = null_mean, xmin = null_low, xmax = null_high), orientation = "y", width = .14, linewidth = .55, colour = "#89929C") +
  geom_point(data = null, aes(x = null_mean, colour = "Null mean"), size = 1.8) +
  geom_point(data = null, aes(x = observed, colour = "Observed"), shape = 18, size = 2.7) +
  geom_text(data = null, aes(x = 352, label = q_label), colour = pal[["ink"]], size = 2.45, hjust = 1) +
  scale_colour_manual(values = c(Observed = pal[["blue"]], `Null mean` = "#747F8B"), breaks = c("Observed", "Null mean")) +
  scale_x_continuous(limits = c(235, 354), breaks = c(250, 275, 300, 325, 350), expand = expansion(mult = c(0, .015))) +
  labs(x = "Sign-stable axes among 640 features", y = NULL, colour = NULL, title = "Sign stability against the permutation null") +
  figure4_theme + theme(legend.position = "bottom", axis.line.y = element_blank(), axis.ticks.y = element_blank())

fig4 <- p4a + p4b + p4c + p4d + plot_layout(design = "AABB\nCCCC\nDDDD", widths = c(1, 1, 1.2, 1.2), heights = c(.92, 1.14, .98)) +
  plot_annotation(tag_levels = "a")
save_pub(fig4, file.path(main_dir, "Figure4_mechanism_axes"), 183, 149)

# Figure 5 -------------------------------------------------------------------
spcells <- read_src("figure5", "representative_spatial_cells")
spatial_x <- range(spcells$x)
spatial_y <- range(spcells$y)
bar_x <- spatial_x[1] + 80
bar_y <- spatial_y[1] - 130
p5a <- ggplot(spcells, aes(x, y, colour = intrinsic_resistance_ecosystem)) +
  geom_point(size = .45, alpha = .95, stroke = 0) +
  scale_colour_gradient2(low = "#3D739F", mid = "#D5DADF", high = "#B7504C", midpoint = 0,
    limits = quantile(spcells$intrinsic_resistance_ecosystem, c(.02, .98), na.rm = TRUE), oob = squish,
    breaks = c(-.5, 0, .5, 1), guide = guide_colourbar(title.position = "top", title.hjust = .5,
      barwidth = unit(27, "mm"), barheight = unit(2.7, "mm"))) +
  annotate("segment", x = bar_x, xend = bar_x + 500, y = bar_y, yend = bar_y, linewidth = .9, colour = pal[["ink"]]) +
  annotate("text", x = bar_x + 250, y = bar_y - 100, label = "500 \u00b5m", size = 2.6, colour = pal[["ink"]]) +
  coord_fixed(xlim = spatial_x + c(-70, 70), ylim = spatial_y + c(-280, 60), expand = FALSE) +
  labs(x = NULL, y = NULL, colour = "Resistance-ecosystem score", title = "Baseline region (su001)") +
  theme_void(base_size = 7.5, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 8.1, face = "bold", hjust = .5, margin = margin(b = 4)),
    legend.position = "bottom", legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
    plot.margin = margin(4, 5, 2, 3))

moran <- read_src("figure5", "spatial_autocorrelation") %>% mutate(feature = clean_fate(sub("^intrinsic_", "", feature)), label = sprintf("%.2f", moran_i))
p5b <- ggplot(moran, aes(feature, sample, fill = moran_i)) + geom_tile(colour = "white", linewidth = .35) + geom_text(aes(label = label), size = 1.35) +
  scale_fill_gradient(low = "white", high = pal[["blue"]], limits = c(0, max(moran$moran_i))) + labs(x = NULL, y = NULL, fill = "Moran's I") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6))

niches <- read_src("figure5", "common_spatial_niches") %>% pivot_longer(c(apoptosis_pressure, proliferation, differentiation_proxy, stress, resistance_ecosystem), names_to = "axis", values_to = "score") %>%
  mutate(axis = clean_fate(axis), niche = paste0("N", spatial_niche), label = sprintf("%.2f", score))
p5c <- ggplot(niches, aes(axis, niche, fill = score)) + geom_tile(colour = "white", linewidth = .4) + geom_text(aes(label = label), size = 1.45) +
  scale_fill_gradient2(low = pal[["blue"]], mid = "white", high = pal[["red"]], midpoint = 0) + labs(x = NULL, y = NULL, fill = "Mean relative\naxis score") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

paired <- read_src("figure5", "paired_future_targets") %>% select(patient, starts_with("delta_absolute_")) %>%
  pivot_longer(-patient, names_to = "axis", values_to = "delta") %>% mutate(axis = clean_fate(sub("^delta_absolute_", "", axis)))
p5d <- ggplot(paired, aes(delta, axis, colour = patient)) + geom_vline(xintercept = 0, linetype = 2, colour = pal[["grey"]], linewidth = .3) +
  geom_point(size = 1.7) + scale_colour_manual(values = c(su001 = pal[["blue"]], su002 = pal[["orange"]], su007 = pal[["teal"]])) +
  labs(x = "Post minus baseline absolute score", y = NULL, colour = "Patient") + theme(legend.position = "right")

joint <- read_src("figure5", "joint_model_metrics") %>% filter(metric == "mae") %>% mutate(model = clean_model(model))
p5e <- ggplot(joint, aes(reorder(model, value), value, fill = model)) + geom_col(width = .62) + geom_text(aes(label = sprintf("%.1f", value)), hjust = -0.12, size = 1.8) +
  coord_flip() + scale_fill_manual(values = c("intrinsic" = pal[["blue2"]], "joint" = pal[["violet"]], "spatial niche" = pal[["teal2"]], "training mean" = pal[["light"]])) +
  scale_y_continuous(expand = expansion(mult = c(0,.14))) + labs(x = NULL, y = "LOPO MAE (n = 4 patients)") + theme(legend.position = "none")

sqc <- read_src("figure5", "spatial_qc_sensitivity") %>% filter(comparison %in% c("moran_concordance_vs_primary", "patient_absolute_axis_concordance_vs_primary")) %>%
  mutate(comparison = recode(comparison, moran_concordance_vs_primary = "Moran", patient_absolute_axis_concordance_vs_primary = "Patient axes"))
p5f <- ggplot(sqc, aes(label, value, group = comparison, colour = comparison)) + geom_hline(yintercept = 1, colour = pal[["light"]]) + geom_line(linewidth = .55) + geom_point(size = 1.5) +
  scale_colour_manual(values = c(Moran = pal[["blue"]], `Patient axes` = pal[["teal"]])) + coord_cartesian(ylim = c(.95,1.005)) +
  labs(x = "Minimum transcripts (primary = 5)", y = "Concordance with primary", colour = NULL) + theme(legend.position = "bottom")

p5a <- wrap_elements(full = p5a)
p5d <- p5d + theme(legend.position = "bottom", legend.direction = "horizontal")
p5f <- p5f + labs(y = "Concordance")
fig5 <- p5a + p5b + p5c + p5d + p5e + p5f +
  plot_layout(design = "AABB\nAACC\nDDDD\nEEFF", widths = c(1, 1, 1.05, 1.05), heights = c(1.05, 1, .78, .85)) +
  plot_annotation(tag_levels = "a")
save_pub(fig5, file.path(main_dir, "Figure5_spatial_niches"), 183, 179)

# Extended Data 2 ------------------------------------------------------------
ed2 <- read_src("extended_data2", "mixseq_complete_sensitivity") %>% filter(metric=="brier") %>%
  select(variant, min_cells, mad_multiplier, validation, target, model, mean_value) %>% pivot_wider(names_from=model, values_from=mean_value) %>%
  mutate(improvement=training_mean-baseline_scvi_ridge, fate=clean_fate(target), setting=paste0("n",min_cells," / ",mad_multiplier," MAD"), validation=gsub("_","-",validation))
p_ed2 <- ggplot(ed2,aes(fate,setting,fill=improvement)) + geom_tile(colour="white",linewidth=.35) + facet_wrap(~validation,nrow=1) +
  scale_fill_gradient2(low=pal[["red"]],mid="white",high=pal[["blue"]],midpoint=0) + labs(x=NULL,y="Target definition",fill="Brier improvement") +
  theme(axis.text.x=element_text(angle=35,hjust=1))
save_pub(p_ed2, file.path(ext_dir,"ExtendedData2_mixseq_sensitivity"), 183, 105)

# Extended Data 3 ------------------------------------------------------------
ed3 <- read_src("extended_data3", "cross_task_dataset_metrics") %>% filter(model!="training_mean") %>% mutate(target=clean_fate(target),validation=gsub("_"," ",validation))
ed3sum <- ed3 %>% group_by(validation,target,model,metric) %>% summarise(value=mean(value,na.rm=TRUE),.groups="drop") %>% filter(metric %in% c("spearman","mse"))
p_ed3 <- ggplot(ed3sum,aes(target,interaction(validation,metric,sep=" | "),fill=value)) + geom_tile(colour="white",linewidth=.35) + geom_text(aes(label=sprintf("%.2f",value)),size=1.35) +
  facet_wrap(~model,nrow=1) + scale_fill_gradient2(low=pal[["red"]],mid="white",high=pal[["blue"]],midpoint=0) + labs(x=NULL,y=NULL,fill="Mean") + theme(axis.text.x=element_text(angle=40,hjust=1),axis.text.y=element_text(size=6))
save_pub(p_ed3, file.path(ext_dir,"ExtendedData3_cross_transfer"), 183, 115)

# Extended Data 4 ------------------------------------------------------------
ed4 <- bind_rows(read_src("extended_data4","replogle_k562_fold_metrics") %>% mutate(cell_line="K562"), read_src("extended_data4","replogle_rpe1_fold_metrics") %>% mutate(cell_line="RPE1")) %>%
  mutate(fate=clean_fate(target),model=clean_model(model),validation=gsub("_","-",validation))
ed4sum <- mean_ci(ed4,c("cell_line","validation","fate","model","metric")) %>% filter(is.finite(mean))
p_ed4 <- ggplot(ed4sum,aes(model,mean,colour=model)) + geom_point(size=1.2) + geom_errorbar(aes(ymin=mean-ci,ymax=mean+ci),width=.15) +
  facet_grid(cell_line+validation~metric+fate,scales="free_y") + scale_colour_manual(values=c("baseline + perturbation"=pal[["blue"]],"perturbation descriptor"=pal[["teal"]],"training mean"=pal[["light"]])) +
  labs(x=NULL,y="Fold mean ± 95% CI",colour=NULL) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="bottom",strip.text=element_text(size=6))
save_pub(p_ed4, file.path(ext_dir,"ExtendedData4_replogle_full"), 183, 145)

# Extended Data 5 ------------------------------------------------------------
cov <- read_src("extended_data5","panel_module_coverage") %>% distinct(module_type,module,genes_specified,genes_present) %>% mutate(fraction=genes_present/genes_specified,module=clean_fate(module))
p_ed5a <- ggplot(cov,aes(fraction,reorder(module,fraction),fill=module_type))+geom_col(width=.65)+scale_x_continuous(labels=percent,limits=c(0,1))+scale_fill_manual(values=c(cell_type=pal[["blue2"]],fate=pal[["teal2"]]))+labs(x="Targeted-panel marker coverage",y=NULL,fill=NULL)+theme(legend.position="bottom")
sq <- read_src("extended_data5","sample_qc") %>% mutate(retained=cells_qc/cells_input)
p_ed5b <- ggplot(sq,aes(retained,reorder(sample,retained),fill=patient))+geom_col(width=.65)+scale_fill_manual(values=c(su001=pal[["blue"]],su002=pal[["orange"]],su006=pal[["violet"]],su007=pal[["teal"]]))+scale_x_continuous(labels=percent,limits=c(0,1))+labs(x="Cell retention",y=NULL,fill="Patient")+theme(legend.position="bottom")
sens5 <- read_src("extended_data5","spatial_qc_sensitivity") %>% filter(comparison=="retained_cells") %>% mutate(label=factor(label,levels=c("min5","min10","min20")))
p_ed5c <- ggplot(sens5,aes(label,value,group=1))+geom_line(colour=pal[["blue"]],linewidth=.5)+geom_point(colour=pal[["blue"]],size=1.4)+labs(x="Minimum transcript threshold",y="Retained cells")
edfig5 <- p_ed5a+p_ed5b+p_ed5c+plot_layout(widths=c(1.3,1,1))+plot_annotation(tag_levels="a")
save_pub(edfig5,file.path(ext_dir,"ExtendedData5_spatial_audit"),183,105)

# Clinical comparisons replace response-supervised cluster displays.
dat <- function(name) read.delim(file.path(project,"results",paste0(name,".tsv")),check.names=FALSE)
cohort_name <- function(x) recode(x,ispy2="I-SPY2",beataml="BeatAML")
model_name <- function(x) recode(x,PFE="PFE",clinical="Clinical",clinical_plus_PFE="Clinical + PFE",five_gene="Five genes",clinical_plus_five_gene="Clinical + five genes")
clinical_reference <- dat("clinical_bootstrap")
primary_model_map <- c(PFE="PFE", five_gene="five_gene", hr_her2_plus_arm="clinical", hr_her2_plus_arm_plus_PFE="clinical_plus_PFE", hr_her2_plus_arm_plus_five_gene="clinical_plus_five_gene")
clinical_adjusted <- dat("ispy2_arm_adjustment_pooled_metrics") %>% filter(model %in% names(primary_model_map)) %>%
  mutate(cohort="ispy2", model=unname(primary_model_map[model])) %>% select(all_of(names(clinical_reference)))
clinical <- bind_rows(clinical_adjusted, clinical_reference %>% filter(cohort=="beataml")) %>%
  filter(metric=="auroc") %>% mutate(cohort=cohort_name(cohort),model=model_name(model))
p6b <- ggplot(clinical,aes(estimate,model,colour=cohort))+geom_vline(xintercept=.5,linetype=2,colour=pal[["light"]])+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=.18,position=position_dodge(width=.45),linewidth=.4)+geom_point(size=1.4,position=position_dodge(width=.45))+scale_colour_manual(values=c(`I-SPY2`=pal[["blue"]],BeatAML=pal[["teal"]]))+labs(x="Pooled out-of-fold AUROC (95% CI)",y=NULL,colour=NULL)+theme(legend.position="bottom",axis.text.y=element_text(size=6.4))
ann <- dat("clinical_annotation_response") %>% mutate(cohort=cohort_name(cohort),label=gsub("−","-",label)) %>% filter(n>=10)
p6c <- ggplot(ann,aes(label,fraction,fill=cohort))+geom_col(width=.65)+geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.15,linewidth=.35)+geom_text(aes(label=paste0("n=",n)),y=.98,size=2)+facet_wrap(~cohort,scales="free_x")+scale_fill_manual(values=c(`I-SPY2`=pal[["blue2"]],BeatAML=pal[["teal2"]]))+scale_y_continuous(labels=percent,limits=c(0,1.04))+labs(x=NULL,y="Observed response fraction (95% CI)")+theme(legend.position="none",axis.text.x=element_text(angle=35,hjust=1,size=6))
cv <- dat("clinical_fold_metrics") %>% filter(model=="PFE") %>% mutate(cohort=cohort_name(cohort),metric=recode(metric,auroc="AUROC",average_precision="AP",balanced_accuracy="Bal. accuracy"))
chance <- data.frame(cohort=c("I-SPY2","BeatAML"),metric="AP",reference=c(319/987,217/317))
p6d <- ggplot(cv,aes(cohort,value,colour=cohort))+geom_hline(data=chance,aes(yintercept=reference,colour=cohort),linetype=2,linewidth=.4)+geom_boxplot(width=.52,outlier.shape=NA,linewidth=.4)+geom_point(position=position_jitter(width=.08,seed=20261116),size=.75)+facet_wrap(~metric)+scale_colour_manual(values=c(`I-SPY2`=pal[["blue"]],BeatAML=pal[["teal"]]))+labs(x=NULL,y="Fivefold PFE response metric")+theme(legend.position="none",axis.text.x=element_text(angle=25,hjust=1,size=6))
surv <- dat("clinical_survival") %>% filter(term=="response_coordinate") %>% mutate(label=ifelse(model=="unadjusted",paste0("Unadjusted (n=",n,")"),paste0("Age + ELN adjusted (n=",n,")")))
p6e <- ggplot(surv,aes(estimate,label))+geom_vline(xintercept=1,linetype=2,colour=pal[["grey"]])+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=.16,linewidth=.45,colour=pal[["teal"]])+geom_point(size=1.8,colour=pal[["teal"]])+scale_x_log10(breaks=c(.6,.8,1,1.2))+labs(x="BeatAML survival hazard ratio per SD (95% CI)",y=NULL)
fig6 <- p6b+p6c+p6d+p6e+plot_layout(design="AAAA\nBBCC\nDDDD",heights=c(1.05,1.25,.65))+plot_annotation(tag_levels="a")
save_pub(fig6,file.path(main_dir,"Figure6_patient_fate_embedding"),183,146)

# Supplementary calibration: binned means describe fractional branch targets.
cal <- dat("calibration") %>% mutate(target=gsub("^p_","",gsub("_"," ",target)))
dist <- dat("target_distributions") %>% mutate(target=gsub("^p_","",gsub("_"," ",target)))
for(ds in unique(cal$dataset[grepl("LARRY",cal$dataset)])) {
 cc <- cal %>% filter(dataset==ds)
 pp <- ggplot(cc,aes(predicted,observed))+geom_abline(slope=1,intercept=0,linetype=2,colour=pal[["grey"]],linewidth=.3)+geom_line(colour=pal[["blue"]],linewidth=.4)+geom_point(aes(size=n),colour=pal[["blue"]],alpha=.8)+facet_wrap(~target,ncol=4)+coord_equal(xlim=c(0,1),ylim=c(0,1))+scale_size(range=c(.4,2),guide="none")+labs(x="Mean predicted descendant fraction",y="Mean observed descendant fraction",title=ds)+theme(strip.text=element_text(size=6),panel.spacing.x=unit(6,"mm"))
 save_pub(pp,file.path(ext_dir,ifelse(grepl("in_vitro",ds),"Supplementary5_larry_in_vitro_calibration","Supplementary6_larry_cytokine_calibration")),183,145)
}
pp <- ggplot(cal %>% filter(dataset=="MIX-seq"),aes(predicted,observed,colour=validation))+geom_abline(slope=1,intercept=0,linetype=2,colour=pal[["grey"]],linewidth=.3)+geom_line(linewidth=.4)+geom_point(size=1)+facet_wrap(~target,ncol=2)+coord_equal(xlim=c(0,1),ylim=c(0,1))+scale_colour_manual(values=c(leave_cell_line=pal[["blue"]],leave_perturbation=pal[["orange"]]),labels=c("Cell line held out","Perturbation held out"))+labs(x="Mean prediction",y="Mean observed target",colour=NULL)+theme(legend.position="bottom",panel.spacing.x=unit(6,"mm"))
save_pub(pp,file.path(ext_dir,"Supplementary7_mixseq_calibration"),160,150)
for(ds in unique(dist$dataset)) {
 dd <- dist %>% filter(dataset==ds, validation %in% c("clone","leave_cell_line"))
 pp <- ggplot(dd,aes(value))+geom_histogram(aes(y=after_stat(count/sum(count))),bins=20,fill=pal[["blue2"]],colour="white",linewidth=.15)+facet_wrap(~target,scales="free",ncol=4)+labs(x="Observed target",y="Fraction of all targets",title=ds)+theme(strip.text=element_text(size=5.8),panel.spacing.x=unit(6,"mm"))
 # Each target distribution is normalized separately within its facet.
 pp <- ggplot(dd,aes(value))+geom_histogram(bins=20,fill=pal[["blue2"]],colour="white",linewidth=.15)+facet_wrap(~target,scales="free",ncol=4)+labs(x="Observed target",y="Number of targets",title=ds)+theme(strip.text=element_text(size=5.8),panel.spacing.x=unit(6,"mm"))
 fname <- if(ds=="MIX-seq") "Supplementary10_mixseq_distributions" else if(grepl("in_vitro",ds)) "Supplementary8_larry_in_vitro_distributions" else "Supplementary9_larry_cytokine_distributions"
 save_pub(pp,file.path(ext_dir,fname),183,if(ds=="MIX-seq") 75 else 160)
}
save_pub(p4e,file.path(ext_dir,"Supplementary11_velocity_eligibility"),130,85)
cat("Revised manuscript figures exported.\n")

extra <- read.delim(file.path(project,"results/additional_calibration.tsv"))
extra_dist <- read.delim(file.path(project,"results/additional_distributions.tsv"))
curve <- function(d) ggplot(d,aes(predicted,observed,colour=validation))+geom_abline(slope=1,intercept=0,linetype=2,colour=pal[["grey"]],linewidth=.3)+geom_line(linewidth=.4)+geom_point(size=.7)+coord_equal(xlim=c(0,1),ylim=c(0,1))+scale_colour_manual(values=c(replicate=pal[["blue"]],leave_cell_line=pal[["blue"]],leave_perturbation=pal[["blue"]],leave_gem_group=pal[["orange"]]))+labs(x="Mean prediction",y="Mean observed fraction",colour=NULL)+theme(legend.position="bottom",panel.spacing.x=unit(6,"mm"))
wc <- curve(extra %>% filter(dataset=="Waddington-OT"))+facet_wrap(~target,ncol=3)+guides(colour="none")
wd <- ggplot(extra_dist %>% filter(dataset=="Waddington-OT"),aes(value))+geom_histogram(bins=20,fill=pal[["blue2"]],colour="white",linewidth=.15)+facet_wrap(~target,ncol=3,scales="free")+labs(x="Observed transport target",y="Targets")
save_pub(wc/wd+plot_annotation(tag_levels="a"),file.path(ext_dir,"Supplementary12_transport_calibration"),183,210)
gc <- curve(extra %>% filter(dataset!="Waddington-OT"))+facet_grid(dataset~target)+scale_colour_manual(values=c(leave_perturbation=pal[["blue"]],leave_gem_group=pal[["orange"]]),labels=c(leave_perturbation="Perturbation held out",leave_gem_group="GEM held out"))
gd <- ggplot(extra_dist %>% filter(dataset!="Waddington-OT"),aes(value))+geom_histogram(bins=20,fill=pal[["blue2"]],colour="white",linewidth=.15)+facet_grid(dataset~target,scales="free_y")+labs(x="Observed program fraction",y="Pairs")
save_pub(gc/gd+plot_annotation(tag_levels="a"),file.path(ext_dir,"Supplementary13_genetic_calibration"),183,205)
