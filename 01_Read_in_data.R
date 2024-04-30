library(readxl)
library(tidyverse)
library(here)

path <- "data/train_dataset.xlsx"
path2 <- "data/reference values.xlsx"
path3 <- "data/Reference test dataset.xlsx"
sheets <- excel_sheets(path)
labelsTest <- readxl::read_excel(path3, sheet = excel_sheets(path3))

specTrain <- NULL
for(i in 1:length(sheets)){
  spec <- readxl::read_excel(path, sheet = sheets[i])
  spec <- spec |>
    mutate(`sample number` = rep(sheets[i], nrow(spec)))
  specTrain <- rbind(specTrain, spec)
}

labels <- readxl::read_excel(path2, sheet = excel_sheets(path2))
specTrain |> glimpse()
specTrain <- specTrain |> left_join(labels)

# labels$`lactose content` |> table()
tbl <- sort(unique(labels$`lactose content`))

for (i in 1:length(tbl)){
  p1 <- specTrain |>
    mutate(ID = 1:nrow(specTrain)) |>
    # mutate(`sample number` = fct_reorder(as.factor(`sample number`), `lactose content`)) |>
    pivot_longer(-c(ID, `sample number`, `lactose content`)) |>
    mutate(name = as.numeric(name)) |>
    filter(`lactose content`==tbl[i]) |>
    ggplot(aes(x = name, y = value, group = ID, col = `sample number`)) +
    geom_line()+
    facet_wrap(~`sample number`)
  print(p1)
}

# for(i in 1:length(sheets)){
# spec <- readxl::read_excel(path, sheet = sheets[i])
#
# p1 <- spec |>
#   mutate(ID = 1:nrow(spec)) |>
#   pivot_longer(-ID) |>
#   mutate(name = as.numeric(name)) |>
#   ggplot(aes(x = name, y = value, group = ID, col = ID)) +
#   geom_line()
# print(p1)
# }
#
