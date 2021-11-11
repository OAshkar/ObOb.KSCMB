library(tidyverse)
library(ComplexHeatmap)

#####################################################################################################
# Model Structure
structComp <- scale(read_csv("./MatResults/structComp.csv", col_names = F))
rownames_structComp <- read_csv("./MatResults/structComp_rowcolnames.csv", col_names = F)$X1 %>%
  as.vector()


column_ha <- HeatmapAnnotation(
  Diet = gsub(x = colnames_structComp, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\2"),
  mouse = gsub(x = colnames_structComp, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\1"),
  tissue = gsub(x = colnames_structComp, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\3")
)
png("./paper/fig/structcomp_cluster.png", units = "in", width = 20, height = 12, res = 400)
Heatmap(
  t(scale(t(structComp))),
  row_labels = rownames_structComp,
  column_labels = rownames_structComp,
  bottom_annotation =  column_ha
)
dev.off()

#####################################################################################################
#####################################################################################################
# Subsystems coverage
subCoverage <- scale(read_csv("MatResults/subCoverage.csv", col_names = F))
rownams_subCoverage <- read_csv("./MatResults/subCoverage_rownames.csv", col_names = F)$X1 %>%
  as.vector()
colnames_subCoverage <- read_csv("./MatResults/subCoverage_colnames.csv",
  col_names = F
)[1, ] %>%
  as.vector() %>%
  as.character()

column_ha <- HeatmapAnnotation(
  Diet = gsub(x = colnames_subCoverage, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\2"),
  mouse = gsub(x = colnames_subCoverage, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\1"),
  tissue = gsub(x = colnames_subCoverage, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\3")
)
png("./paper/fig/subsystemcoverage.jpg", units = "in", width = 20, height = 12, res = 400)
Heatmap(
  t(scale(t(subCoverage))),
  row_labels = rownams_subCoverage,
  column_labels = colnames_subCoverage,
  bottom_annotation =  column_ha
)
dev.off()
#####################################################################################################
# Metabolic Tasks Evaluation
metFunction <- read_csv("./MatResults/MetabolicFunction.csv", col_names = F) %>%
  mutate(
    Task = read_csv("./MatResults/MetabolicFunction_rownames.csv", col_names = F)$X1 %>%
      as.vector()
  )
x <- read_csv("./MatResults/MetabolicFunction_colnames.csv",
  col_names = F
)[1, ] %>%
  as.vector() %>%
  as.character()
colnames(metFunction) <- c(x, "Task")

metFunction2 <- metFunction %>%
  pivot_longer(names_to = "Sample", values_to = "Presence", cols = -Task) %>%
  mutate(
    Diet = gsub(.$Sample, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\2"),
    mouse = gsub(.$Sample, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\1"),
    tissue = gsub(.$Sample, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\3"),
    Presence = na_if(.$Presence, 0)
  )


png("./paper/fig/funcCompMatrix.jpg", units = "in", width = 20, height = 12, res = 400)
ggplot(metFunction2, aes(x = Sample, y = Task, size = Presence, color = Diet)) +
  geom_point(show.legend = T, drop = T) +
  facet_grid(~tissue, scales = "free", space = "free") +
  scale_y_discrete(position = "right") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = 18)
  ) +
  scale_size_continuous(guide = FALSE)
dev.off()

# elrm of task function
library(elrm)

# dataframe with 4 columns. func, n, diet, geno
# need to repeat the command to zero

funcDat <- metFunction %>%
  pivot_longer(names_to = "Sample", values_to = "Presence", cols = -Task) %>%
  mutate(
    Diet = gsub(.$Sample, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\2"),
    mouse = gsub(.$Sample, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\1"),
    tissue = gsub(.$Sample, pattern = ".(WT|ob.ob).(HFD|ND).(.*)", replacement = "\\3"),
    n = 1
  )




# n is the no of trial. Make col = 1
# y is the outcome

test <- funcData %>%
  group_by(tissue, mouse, diet) %>%
  filter(funcDat, Task == "Anaerobic rephosphorylation of ATP")

crash.elrm <- elrm(
  formula = Presence / n ~ Diet + mouse + tissue,
  interest = ~mouse,
  r = 4, iter = 20000,
  dataset = test,
  burnIn = 100
)
summary(crash.elrm)

mod <- glm(Presence ~ Diet * mouse * tissue,
  family = binomial,
  data = test
)
summary(mod)