# TCGA Kidney Tumor Classification
# Revised R Script - Self-contained and executable

# ---------------------------
# 1. Package Installation & Loading
# ---------------------------
packages <- c("tidyverse", "caret", "randomForest", "nnet", "e1071", "kernlab", "ggfortify", "biomaRt")
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}
lapply(packages, install_if_missing)

# ---------------------------
# 2. Load and Merge Expression Data
# NOTE:
# This script assumes the following files exist in the ./data folder:
# - expr_kirc_df.rds
# - expr_kirp_df.rds
# - expr_kich_df.rds
# Please download them from 'https://github.com/narisarit/edx-capstone-validation.git' and place them in the `data/` folder before running this script.
# These datasets were originally downloaded from the GDC Data Portal by Dong Hyun Kang.
# They contain transcriptomic TPM data from TCGA kidney cancer samples, pathologically classified as KIRC, KIRP, and KICH.
# ---------------------------
expr_kirc_df <- readRDS("data/expr_kirc_df.rds")
expr_kirp_df <- readRDS("data/expr_kirp_df.rds")
expr_kich_df <- readRDS("data/expr_kich_df.rds")

expr_all_wide <- bind_rows(expr_kirc_df, expr_kirp_df, expr_kich_df) %>%
  relocate(sample_id, tumor_type)

# ---------------------------
# 3. PCA-Based Outlier Removal (PC1)
# ---------------------------
# Remove zero-variance columns (e.g., sample_id, tumor_type mistakenly treated as numeric)
numeric_expr <- expr_all_wide %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(where(~ sd(.) > 0))

# Perform PCA
pca_all <- prcomp(numeric_expr, scale. = TRUE)

# Plot initial PC1 vs PC2 to visualize outliers
pca_all_df <- as_tibble(pca_all$x) %>%
  mutate(sample_id = expr_all_wide$sample_id, tumor_type = expr_all_wide$tumor_type)

ggplot(pca_all_df, aes(x = PC1, y = PC2, color = tumor_type)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Initial PCA: PC1 vs PC2 before outlier removal")

# Identify PC1 outliers
pc1_scores <- pca_all$x[, 1]
thresh_hi <- mean(pc1_scores) + 3 * sd(pc1_scores)
thresh_lo <- mean(pc1_scores) - 3 * sd(pc1_scores)
extreme_pc1 <- expr_all_wide$sample_id[pc1_scores > thresh_hi | pc1_scores < thresh_lo]

# Filter clean data
expr_clean <- expr_all_wide %>% filter(!sample_id %in% extreme_pc1)

numeric_expr_clean <- expr_clean %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(where(~ sd(.) > 0))

# ---------------------------
# 4. PCA and Top Genes from PC3/PC4
# Visualize why PC3 and PC4 were chosen by comparing variance separation
pca_clean <- prcomp(numeric_expr_clean, scale. = TRUE)

pca_clean_df <- as_tibble(pca_clean$x) %>%
  mutate(tumor_type = expr_clean$tumor_type)

# Plot PC1 vs PC2
print(
  ggplot(pca_clean_df, aes(x = PC1, y = PC2, color = tumor_type)) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(title = "PCA after outlier removal: PC1 vs PC2")
)

# Plot PC2 vs PC3
print(
  ggplot(pca_clean_df, aes(x = PC2, y = PC3, color = tumor_type)) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(title = "PCA after outlier removal: PC2 vs PC3")
)

# Plot PC3 vs PC4
print(
  ggplot(pca_clean_df, aes(x = PC3, y = PC4, color = tumor_type)) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(title = "PCA after outlier removal: PC3 vs PC4")
)

# select top genes
rotation <- pca_clean$rotation
top_pc3_genes <- sort(abs(rotation[, 3]), decreasing = TRUE)[1:10]
top_pc4_genes <- sort(abs(rotation[, 4]), decreasing = TRUE)[1:10]
selected_genes <- unique(c(names(top_pc3_genes), names(top_pc4_genes)))

# ---------------------------
# 5. Prepare Modeling Data
# ---------------------------
model_data <- expr_clean %>%
  dplyr::select(all_of(selected_genes), tumor_type) %>%
  mutate(across(all_of(selected_genes), ~log2(. + 1)))

set.seed(42)
split_idx <- createDataPartition(model_data$tumor_type, p = 0.8, list = FALSE)
train <- model_data[split_idx, ]
test <- model_data[-split_idx, ]
train$tumor_type <- factor(train$tumor_type)
test$tumor_type <- factor(test$tumor_type)

# ---------------------------
# 6. Model Training: RF, Logit, SVM (with kernlab)
# ---------------------------
# Random Forest
rf_model <- randomForest(tumor_type ~ ., data = train)
pred_rf <- predict(rf_model, test)
conf_rf <- confusionMatrix(pred_rf, test$tumor_type)

# Logistic Regression
logit_model <- multinom(tumor_type ~ ., data = train)
pred_logit <- predict(logit_model, test)
conf_logit <- confusionMatrix(pred_logit, test$tumor_type)

# SVM via kernlab
train_matrix <- as.matrix(train[, setdiff(names(train), "tumor_type")])
train_labels <- train$tumor_type

svm_model <- ksvm(x = train_matrix, y = train_labels, type = "C-svc", kernel = "vanilladot", scaled = TRUE)
test_matrix <- as.matrix(test[, setdiff(names(test), "tumor_type")])
pred_svm <- predict(svm_model, newdata = test_matrix)
conf_svm <- confusionMatrix(pred_svm, test$tumor_type)

# ---------------------------
# 6.1 Top Gene Comparison Across Models
# ---------------------------
# Compare top 10 important genes from each model
rf_importance <- importance(rf_model)[, 1]
top10_rf <- sort(rf_importance, decreasing = TRUE)[1:10]

logit_importance <- apply(abs(summary(logit_model)$coefficients), 2, max)
top10_logit <- sort(logit_importance, decreasing = TRUE)[1:10]

W <- colSums(coef(svm_model)[[1]] * svm_model@xmatrix[[1]])
names(W) <- colnames(train_matrix)
top10_svm <- sort(abs(W), decreasing = TRUE)[1:10]

# Intersect top genes
top10_overlap <- Reduce(intersect, list(names(top10_rf), names(top10_logit), names(top10_svm)))
top10_overlap  # should contain the 3 selected genes

# ---------------------------
# 7. Minimal 3-Gene Signature Model
# The following three genes were selected because they appeared in the top-ranked importance lists of all three models:
# Random Forest, Logistic Regression, and SVM.
# - ENSG00000167646.14 → DNAAF3
# - ENSG00000169727.12 → GPS1
# - ENSG00000204237.5  → OXLD1
# This overlap was determined after evaluating model-based variable importance scores.
# ---------------------------
top3_genes <- c("ENSG00000167646.14", "ENSG00000169727.12", "ENSG00000204237.5")
model_top3 <- expr_clean %>%
  dplyr::select(all_of(top3_genes), tumor_type) %>%
  mutate(across(all_of(top3_genes), ~log2(. + 1)))

split_idx3 <- createDataPartition(model_top3$tumor_type, p = 0.8, list = FALSE)
train3 <- model_top3[split_idx3, ]
test3  <- model_top3[-split_idx3, ]
train3$tumor_type <- factor(train3$tumor_type)
test3$tumor_type  <- factor(test3$tumor_type)

# Random Forest on top 3 genes
rf3 <- randomForest(tumor_type ~ ., data = train3)
pred_rf3 <- predict(rf3, test3)
conf_rf3 <- confusionMatrix(pred_rf3, test3$tumor_type)

# Logistic Regression on top 3 genes
logit3 <- multinom(tumor_type ~ ., data = train3)
pred_logit3 <- predict(logit3, test3)
conf_logit3 <- confusionMatrix(pred_logit3, test3$tumor_type)

# SVM via kernlab on top 3 genes
train3_matrix <- as.matrix(train3[, setdiff(names(train3), "tumor_type")])
test3_matrix <- as.matrix(test3[, setdiff(names(test3), "tumor_type")])
svm3 <- ksvm(x = train3_matrix, y = train3$tumor_type, type = "C-svc", kernel = "vanilladot", scaled = TRUE)
pred_svm3 <- predict(svm3, newdata = test3_matrix)
conf_svm3 <- confusionMatrix(pred_svm3, test3$tumor_type)

# ---------------------------
# 8. Summary of Results
# ---------------------------
results <- list(
  Full = list(RandomForest = conf_rf, Logistic = conf_logit, SVM = conf_svm),
  Top3 = list(RandomForest = conf_rf3, Logistic = conf_logit3, SVM = conf_svm3)
)

sapply(results$Full, function(x) x$overall["Accuracy"])
sapply(results$Top3, function(x) x$overall["Accuracy"])


# View summary accuracy

# ---------------------------
# 9. References
# ---------------------------
# - TCGA (The Cancer Genome Atlas): https://www.cancer.gov/ccg/research/genome-sequencing/tcga
# - GDC Data Portal: https://portal.gdc.cancer.gov/
# - Bioconductor biomaRt package: https://bioconductor.org/packages/release/bioc/html/biomaRt.html
# - caret R package: https://topepo.github.io/caret/
# - kernlab R package: https://cran.r-project.org/web/packages/kernlab/index.html
# - randomForest R package: https://cran.r-project.org/web/packages/randomForest/
# - nnet R package: https://cran.r-project.org/web/packages/nnet/

