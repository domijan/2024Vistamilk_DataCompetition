
library(tidyverse)
library(broom)
library(GGally)
library(gridExtra)
library(RColorBrewer)
substructure <- read.csv("surface.csv")
substructure |> glimpse()
# substructure <- substructure |> rename(`sample number` = Sample)
substructure <- substructure |> rename(sample = Sample)

specTrain <- readRDS("specTrain.rds")
specTest <- readRDS("specTest.rds")
medoiddat <- readRDS("medoiddat.rds")
medoidtest <- readRDS("medoidtest.rds")

medoiddat <- medoiddat  |> select(-`lactose content`)
# medoidtest |> glimpse()
# specTest |> glimpse()
#
# specTrain |> dim()
# specTest |> dim()

medoiddat <- medoiddat |> left_join(substructure)
medoidtest <- medoidtest |> left_join(substructure)
medoids <- rbind(medoiddat, medoidtest)
medoids |> glimpse()
medoids <- medoids |> mutate(label = gsub('sample', '', sample))
medoids <- medoids |> mutate(substruct = interaction(Ml, Surface))
lactose <- read.csv("data/CompleteReferences.csv")
medoids <- medoids |> left_join(lactose)



pc <- medoids |> select(`412.71`: `3907.24`) |>
  prcomp(center = TRUE,
         scale. = TRUE)
p1 <- pc |>
  augment(medoids) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = (Surface), label = label)) +
  geom_text()+
  ylab("PC2") +
  xlab("PC1") +
theme(legend.position = "none")

p2 <- pc |>
  augment(medoids) |>
  ggplot(aes(.fittedPC1, .fittedPC3, color =(Surface), label = label)) +
  geom_text()+
  ylab("PC3") +
  xlab("PC1") +
theme(legend.position = "none")

grid.arrange(p1, p2, nrow = 1)


mds_meds <- medoids |> select(`412.71`: `3907.24`) |> dist() |> cmdscale()

# Plot the results
plot(mds_meds[, 1], mds_meds[, 2],
     type = "n", xlab = "MDS Dimension 1",
     ylab = "MDS Dimension 2")

points(mds_meds[, 1], mds_meds[, 2],
       pch = 21, bg = as.factor(medoids$Surface))
text(mds_meds[, 1], mds_meds[, 2],
     labels = substr(medoids$label, 1, 2),
     pos = 3, cex = 0.8)
# pc <- medoiddat |> select(`412.71`: `3907.24`) |>
  # prcomp(center = TRUE,
  #        scale. = TRUE)

# medoiddat <- medoiddat |> mutate(label = gsub('sample', '', sample))

# pc2 <- pc |> predict(medoidtest)
p1 <- pc |>
  augment(medoiddat) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = (Surface))) +
  geom_point()+
  ylab("PC2") +
  xlab("PC1")
  # theme(legend.position = "none")

p2 <- pc |>
  augment(medoiddat) |>
  ggplot(aes(.fittedPC1, .fittedPC3, color =(Surface))) +
  geom_point()+
  ylab("PC3") +
  xlab("PC1")
  # theme(legend.position = "none")

grid.arrange(p1, p2, nrow = 1)
p1 <- pc2 |>
  cbind(medoidtest) |>
  as_tibble() |>
  ggplot(aes(PC1, PC2, color = (Surface))) +
  geom_point()+
  ylab("PC2") +
  xlab("PC1") +
  theme(legend.position = "none")

p2 <- pc2 |>
  cbind(medoidtest) |>
  as_tibble() |>
  ggplot(aes(PC1, PC3, color =(Surface))) +
  geom_point()+
  ylab("PC3") +
  xlab("PC1") +
  theme(legend.position = "none")

grid.arrange(p1, p2, nrow = 1)


# ==================================================

specTrain<- specTrain |> left_join(substructure)
specTest<- specTest |> left_join(substructure)

pc <- specTrain |> select(`412.71`: `3907.24`) |>
  prcomp(center = TRUE,
         scale. = TRUE)

pc |>
  augment(specTrain) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = as.factor(Surface))) +
  geom_point(size = 0.5, alpha = 0.5)

pc |>
  augment(specTrain) |>
  ggplot(aes(.fittedPC1, .fittedPC3, color = as.factor(Surface))) +
  geom_point(size = 0.5, alpha = 0.5)



pc |>
  augment(specTrain) |>
  select(c(.fittedPC1: .fittedPC6, Surface)) |>
  ggpairs(columns = 1:6, aes(colour = as.factor(Surface)),
          lower = list(continuous = wrap("points", alpha = 0.5, size = 0.2)),
          upper = list(continuous = wrap("points", alpha = 0.5, size = 0.2)))


pc |>
  augment(specTrain) |>
  ggplot(aes(.fittedPC1, .fittedPC2, color = as.factor(Ml))) +
  geom_point(size = 0.5, alpha = 0.5)

pc |>
  augment(specTrain) |>
  ggplot(aes(.fittedPC1, .fittedPC3, color = as.factor(Ml))) +
  geom_point(size = 0.5, alpha = 0.5)



pc |>
  augment(specTrain) |>
  select(c(.fittedPC1: .fittedPC6, Ml)) |>
  ggpairs(columns = 1:6, aes(colour = as.factor(Ml)),
          lower = list(continuous = wrap("points", alpha = 0.5, size = 0.2)),
          upper = list(continuous = wrap("points", alpha = 0.5, size = 0.2)))

dat.tr <- specTrain
dat.tr <- dat.tr |> select(-c(`sample number`))
# dat.tr <- dat.tr |> select(-c(`sample number`, Ml, Surface, Replicate))
dat.tr$`lactose content` <- as.factor(dat.tr$`lactose content`)
dat.tr <- dat.tr |> rename(Diet = `lactose content`)
fit_lda <- MASS::lda(Diet ~ .,
                     dat.tr|> select(`412.71`: `3907.24`, Diet))
lda_meds  <-  as_tibble(as.matrix(medoids |> select(`412.71`: `3907.24`)) %*% fit_lda$scaling)
lda_meds <- lda_meds |> mutate(Surface=medoids$Surface,
                               label = medoids$label,
                               Ml = medoids$Ml,
                               substruct = medoids$substruct,
                               lactose = medoids$lactose)
table(lda_meds$Ml, lda_meds$lactose)
p4 <- lda_meds %>%
  ggplot(aes(x = LD1, y = LD2, col = as.factor(lactose), label = label)) +
  geom_text() +
  theme(legend.position = "none")

p5 <- lda_meds %>%
  ggplot(aes(x = LD3, y = LD4, col = as.factor(lactose), label = label)) +
  geom_text() +
  theme(legend.position = "none")

p6 <- lda_meds %>%
  ggplot(aes(x = LD5, y = LD6, col = as.factor(lactose), label = label)) +
  geom_text() +
  theme(legend.position = "none")

p7 <- lda_meds %>%
  ggplot(aes(x = LD2, y = LD7, col = as.factor(lactose), label = label)) +
  geom_text() +
  theme(legend.position = "none")
grid.arrange(p4, p5, p6, p7, nrow = 2)
#
# p5 <- lda_te2 |>
#   ggplot(aes(x = LD1, y = LD2, col = as.factor(PredictionsLDA))) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 30)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nprediction")
#
# p5b <- lda_te2 |>
#   ggplot(aes(x = LD1, y = LD2, col = lactose.content)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 30)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nlabels")
#
# p6 <- lda_tr %>%
#   ggplot(aes(x = LD3, y = LD4, col = Surface)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-15, 10)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Training set\nlabels")
#
#
# p7 <- lda_te %>%
#   mutate(lactose.content = as.factor(lactose.content))|>
#   filter(`lactose.content` %in% c(24,60)) |>
#   ggplot(aes(x = LD3, y = LD4, col = as.factor(PredictionsLDA))) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-15, 10)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nprediction")
#
# p7b <- lda_te %>%
#   mutate(lactose.content = as.factor(lactose.content))|>
#   filter(`lactose.content` %in% c(24,60)) |>
#   ggplot(aes(x = LD3, y = LD4, col = lactose.content)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-15, 10)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nlabels")
#
# p8 <- lda_tr %>%
#   ggplot(aes(x = LD5, y = LD6, col = Surface)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 6)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Training set\nlabels")
#
# p9 <- lda_te %>%
#   mutate(lactose.content = as.factor(lactose.content))|>
#   filter(`lactose.content` %in% c(24,60)) |>
#   ggplot(aes(x = LD5, y = LD6, col = as.factor(PredictionsLDA))) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 6)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nprediction")
#
# p9b <- lda_te %>%
#   mutate(lactose.content = as.factor(lactose.content))|>
#   filter(`lactose.content` %in% c(24,60)) |>
#   ggplot(aes(x = LD5, y = LD6, col = lactose.content)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 6)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nlabels")
#
# grid.arrange(p4, p5, p5b,p6, p7,p7b,  p8, p9, p9b, nrow = 3)


# lda_tr2 <-  as_tibble(as.matrix(dat.tr %>% select(-Diet)) %*% fit_lda$scaling)
# lda_tr2 <- lda_tr2 %>% mutate(Diet = dat.tr$Diet, `sample number` = specTrain$`sample number`, Surface=specTrain$Surface)
# saveRDS(lda_tr2, "lda_tr2.rds")
#
# lda_tr <- readRDS("lda_tr2.rds")
# lda_te <- readRDS("lda_te.rds")
# lda_te <- lda_te |> rename(`sample number` = samp)
# lda_te <- lda_te |> left_join(substructure)
myColors <- brewer.pal(8, "Spectral")
# lda_tr |> glimpse()
# lda_te |> glimpse()
# names(myColors) <- levels(lda_tr$Surface)
custom_colors <- scale_colour_manual(name = "Species Names", values = myColors)
#
#
# p4 <- lda_tr %>%
#   ggplot(aes(x = LD1, y = LD2, col = Surface)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 30)) +
#   theme(legend.position = "none")+
#   custom_colors +
#   ggtitle("Training set\nlabels")
#
# p5 <- lda_te2 |>
#   ggplot(aes(x = LD1, y = LD2, col = as.factor(PredictionsLDA))) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 30)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nprediction")
#
# p5b <- lda_te2 |>
#   ggplot(aes(x = LD1, y = LD2, col = lactose.content)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 30)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nlabels")
#
# p6 <- lda_tr %>%
#   ggplot(aes(x = LD3, y = LD4, col = Surface)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-15, 10)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Training set\nlabels")
#
#
# p7 <- lda_te %>%
#   mutate(lactose.content = as.factor(lactose.content))|>
#   filter(`lactose.content` %in% c(24,60)) |>
#   ggplot(aes(x = LD3, y = LD4, col = as.factor(PredictionsLDA))) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-15, 10)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nprediction")
#
# p7b <- lda_te %>%
#   mutate(lactose.content = as.factor(lactose.content))|>
#   filter(`lactose.content` %in% c(24,60)) |>
#   ggplot(aes(x = LD3, y = LD4, col = lactose.content)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-15, 10)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nlabels")
#
# p8 <- lda_tr %>%
#   ggplot(aes(x = LD5, y = LD6, col = Surface)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 6)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Training set\nlabels")
#
# p9 <- lda_te %>%
#   mutate(lactose.content = as.factor(lactose.content))|>
#   filter(`lactose.content` %in% c(24,60)) |>
#   ggplot(aes(x = LD5, y = LD6, col = as.factor(PredictionsLDA))) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 6)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nprediction")
#
# p9b <- lda_te %>%
#   mutate(lactose.content = as.factor(lactose.content))|>
#   filter(`lactose.content` %in% c(24,60)) |>
#   ggplot(aes(x = LD5, y = LD6, col = lactose.content)) +
#   geom_point(alpha = 0.1) +
#   ylim(c(-10, 10)) +
#   xlim(c(-10, 6)) +
#   theme(legend.position = "none")+
#   custom_colors+
#   ggtitle("Testing set\nlabels")
#
# grid.arrange(p4, p5, p5b,p6, p7,p7b,  p8, p9, p9b, nrow = 3)
