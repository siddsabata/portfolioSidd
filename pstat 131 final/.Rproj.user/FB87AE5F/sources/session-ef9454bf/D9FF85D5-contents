library(janitor)
library(tidyverse)
library(tidymodels)
library(ggplot2)
library(corrplot)
library(ggthemes)
library(kableExtra)
library(themis)
library(ranger)
library(xgboost)
tidymodels_prefer()

pdb_data <- read.csv("data/pdb_data_no_dups.csv") 

protein_data <- pdb_data %>% 
  dplyr::filter(!str_detect(classification, "DNA") & !str_detect(classification, "RNA"))

protein_data_top10 <- protein_data %>% 
  group_by(classification) %>% 
  summarise(count=n()) %>% 
  arrange(-count) %>% 
  slice_head(n=10) %>% 
  select(classification)

# creating working dataframe
protein_data <- protein_data %>% 
  semi_join(protein_data_top10, by='classification') %>% 
  dplyr::filter(experimentalTechnique == 'X-RAY DIFFRACTION') %>% 
  select(-c(crystallizationMethod, pdbxDetails, crystallizationTempK, publicationYear, 
            experimentalTechnique, macromoleculeType, phValue)) %>% 
  na.omit() %>% 
  mutate(classification = factor(classification)) %>% 
  clean_names()

set.seed(1234)
# create training and test datasets 
protein_split <- initial_split(protein_data, prop = 0.8, strata = classification)
protein_train <- training(protein_split)
protein_test <- testing(protein_split)

# k-fold cross validation 
cv_folds <- vfold_cv(protein_train, v = 10, strata = classification)

# create recipe 
protein_recipe <- recipe(classification ~ residue_count + resolution + 
                           structure_molecular_weight + density_matthews + 
                           density_percent_sol, data = protein_train) %>% 
  step_normalize(all_predictors()) %>% 
  step_downsample(classification, under_ratio = 2) %>% 
  step_pca(residue_count, structure_molecular_weight, 
           num_comp = 1, prefix = "rc_smw_component")

bt_mod <- boost_tree(mtry = tune(), 
                     trees = tune(), 
                     learn_rate = tune()) %>%
  set_engine("xgboost") %>% 
  set_mode("classification")

bt_grid <- grid_regular(mtry(range = c(1, 4)), 
                        trees(range = c(200, 600)),
                        learn_rate(range = c(-10, -1)),
                        levels = 10)

bt_wf <- workflow() %>% 
  add_model(bt_mod) %>% 
  add_recipe(protein_recipe)

bt_results <- tune_grid(
  bt_wf,
  resamples = cv_folds,
  grid = bt_grid
)

saveRDS(bt_results, "bt_tuning_results.rds")