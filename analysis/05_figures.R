## 05_figures.R ---------------------------------------------------------------
## Four figures. Each one is meant to answer a question a reviewer would ask,
## not to decorate the report.
##
##   1. What do the trajectories actually look like?
##   2. Should we believe there are four clusters?
##   3. What distinguishes one cluster from another?
##   4. What does the model say, with uncertainty attached?
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(nlme)
})

cl   <- readRDS("data/clusters.rds")
milk <- readRDS("data/prepared.rds")
md   <- readRDS("data/model.rds")

feat <- cl$features
dat  <- merge(milk, feat[, c("cow", "cluster")], by = "cow")
dat$cluster <- factor(dat$cluster)

PAL <- c("#1b6ca8", "#c1582b", "#2e7d4f", "#7b4b94", "#b8932f",
         "#496a81", "#8c3a3a", "#3f7f7f")
pal <- PAL[seq_len(nlevels(dat$cluster))]

theme_set(
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      strip.text       = element_text(face = "bold", hjust = 0),
      plot.title       = element_text(face = "bold", size = 13),
      plot.subtitle    = element_text(colour = "grey35", size = 10),
      plot.caption     = element_text(colour = "grey45", size = 8, hjust = 0),
      legend.position  = "bottom"
    )
)

lab_cluster <- function(x) paste("Cluster", x)

## --- Figure 1: the trajectories -----------------------------------------------
## Every cow in grey behind its own cluster, so the reader can see the spread
## the cluster mean is hiding. A plot of cluster means alone always looks more
## convincing than the data deserves.

cl_mean <- aggregate(protein ~ cluster + week, data = dat, FUN = mean)

## The background layer must NOT carry the cluster column, otherwise
## facet_wrap splits it as well and each panel shows only its own cows in
## grey — which defeats the entire point of drawing a reference herd.
bg <- dat[, c("cow", "week", "protein")]

p1 <- ggplot() +
  geom_line(data = bg, aes(week, protein, group = cow),
            colour = "grey86", linewidth = 0.25) +
  geom_line(data = dat, aes(week, protein, group = cow, colour = cluster),
            linewidth = 0.35, alpha = 0.6) +
  geom_line(data = cl_mean, aes(week, protein, colour = cluster),
            linewidth = 1.3) +
  facet_wrap(~ cluster, labeller = as_labeller(lab_cluster)) +
  scale_colour_manual(values = pal, guide = "none") +
  labs(
    title    = "Milk protein trajectories by cluster",
    subtitle = "Every cow in the herd shown in grey; cows in the panel's cluster in colour; cluster mean in bold",
    x = "Week of lactation", y = "Protein (%)",
    caption  = "Data: nlme::Milk (79 cows, weeks 1-19, three diets)"
  )

ggsave("figures/trajectory_clusters.png", p1, width = 9, height = 6.5, dpi = 200)

## --- Figure 2: should we believe it? ------------------------------------------

sel  <- cl$selection
stab <- cl$stability
nb   <- read.csv("outputs/null_benchmark.csv")

d_sil <- rbind(
  data.frame(k = sel$k, value = sel$silhouette, what = "Observed"),
  data.frame(k = sel$k, value = nb$null_median,  what = "Null (shuffled data)")
)

p2a <- ggplot(d_sil, aes(k, value, colour = what, linetype = what)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = c("Observed" = pal[1],
                                 "Null (shuffled data)" = "grey55")) +
  scale_linetype_manual(values = c("Observed" = "solid",
                                   "Null (shuffled data)" = "dashed")) +
  scale_x_continuous(breaks = sel$k) +
  labs(title = "Average silhouette width",
       subtitle = "Height matters less than the gap between the two lines",
       x = "k", y = "Silhouette", colour = NULL, linetype = NULL)

p2b <- ggplot(sel, aes(k, wss)) +
  geom_line(linewidth = 0.9, colour = pal[1]) +
  geom_point(size = 2, colour = pal[1]) +
  scale_x_continuous(breaks = sel$k) +
  labs(title = "Within-cluster sum of squares",
       subtitle = "The elbow method, shown for completeness — it always slopes down",
       x = "k", y = "WSS")

p2c <- ggplot(stab, aes(factor(cluster), jaccard)) +
  geom_col(aes(fill = jaccard >= 0.75), width = 0.65) +
  geom_hline(yintercept = 0.60, linetype = "dashed", colour = "grey35") +
  geom_hline(yintercept = 0.75, linetype = "dotted", colour = "grey55") +
  geom_text(aes(label = sprintf("%.2f", jaccard)), vjust = -0.5, size = 3.2) +
  scale_fill_manual(values = c("TRUE" = pal[3], "FALSE" = "grey65"),
                    guide = "none") +
  ylim(0, 1) +
  labs(title = "Bootstrap cluster stability",
       subtitle = "Mean Jaccard over 1000 resamples. Dashed 0.60, dotted 0.75 (Hennig, 2007)",
       x = "Cluster", y = "Jaccard")

ggsave("figures/cluster_validation_silhouette.png", p2a, width = 6, height = 4, dpi = 200)
ggsave("figures/cluster_validation_elbow.png",      p2b, width = 6, height = 4, dpi = 200)
ggsave("figures/cluster_validation_stability.png",  p2c, width = 6, height = 4, dpi = 200)

## --- Figure 3: what separates the clusters ------------------------------------

Xz <- as.data.frame(cl$X)
Xz$cluster <- feat$cluster
zmeans <- aggregate(. ~ cluster, data = Xz, FUN = mean)

long <- reshape(zmeans, direction = "long",
                varying = setdiff(names(zmeans), "cluster"),
                v.names = "z", timevar = "feature",
                times = setdiff(names(zmeans), "cluster"),
                idvar = "cluster")

readable <- c(
  level_early  = "Starting level",
  slope_linear = "Rate of change (mid-lactation)",
  curvature    = "Curvature",
  turning_week = "Week of turning point",
  dev_sd       = "Deviation size",
  dev_acf1     = "Deviation persistence",
  dev_skew     = "Deviation skew"
)
long$feature <- factor(readable[long$feature], levels = rev(readable))

p3 <- ggplot(long, aes(factor(cluster), feature, fill = z)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%+.1f", z)), size = 3.1,
            colour = ifelse(abs(long$z) > 0.8, "white", "grey20")) +
  scale_fill_gradient2(low = "#1b6ca8", mid = "white", high = "#c1582b",
                       midpoint = 0, name = "SD from\nherd mean") +
  labs(title = "What separates the clusters",
       subtitle = "Cluster means in standard deviations from the herd average",
       x = "Cluster", y = NULL)

ggsave("figures/cluster_profiles.png", p3, width = 7.5, height = 5, dpi = 200)

## --- Figure 4: the model ------------------------------------------------------

grid <- md$grid

p4 <- ggplot(grid, aes(week, fit, colour = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pal, labels = lab_cluster, name = NULL) +
  scale_fill_manual(values = pal, labels = lab_cluster, name = NULL) +
  labs(
    title    = "Model-implied trajectories with 95% confidence bands",
    subtitle = "Mixed-effects model, random intercept per cow, AR(1) residuals, diet held at reference",
    x = "Week of lactation", y = "Protein (%)",
    caption  = paste(
      "Descriptive only. Clusters were derived from these trajectories, so separation here is guaranteed by construction.",
      "The non-circular test is the cluster-by-diet association in 04_model.R.", sep = "\n")
  )

ggsave("figures/model_predictions.png", p4, width = 8, height = 5.5, dpi = 200)

cat("Wrote 6 figures to figures/\n")
print(list.files("figures"))
