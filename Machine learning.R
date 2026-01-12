# ==============================================================================
# 1. Environment Setup and Library Loading
# ==============================================================================

# Use pacman for efficient package management
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  caret, DALEX, ggplot2, randomForest, kernlab, pROC, xgboost, 
  fs, shapviz, kernelshap, gridExtra, tibble
)

# Create output directory
dir_create("Results")

# Set seed for reproducibility
set.seed(20060414)

# Set working directory (Update this to your local path)
# setwd("C:/Users/whisk/Desktop/NPtotal/Machine learning")

# ==============================================================================
# 2. Data Preparation
# ==============================================================================

# Load datasets
expression_matrix <- read.csv("GSE100927.csv", row.names = 1)
core_genes        <- read.csv("hub gene.csv")[, 2]
sample_info       <- read.csv("group.CSV", row.names = 1)

# Extract and clean expression data for hub genes
data <- expression_matrix[core_genes, , drop = FALSE]
data <- na.omit(data)
data <- as.data.frame(t(data))

# Assign classification labels (ensure "group" exists in your CSV)
data$Type <- as.factor(sample_info[rownames(data), "group"])

# Data Partitioning: 70% Training, 30% Testing
inTrain <- createDataPartition(y = data$Type, p = 0.7, list = FALSE)
train   <- data[inTrain, ]
test    <- data[-inTrain, ]

# ==============================================================================
# 3. Model Training (Iterative Approach)
# ==============================================================================

model_settings <- data.frame(
  AlgorithmName  = c("Random Forest", "SVM", "Neural Network", "LASSO"),
  Implementation = c("rf", "svmRadial", "nnet", "glmnet")
)

# Define 5-fold repeated cross-validation
cv_control <- trainControl(
  method          = "repeatedcv",
  number          = 5,
  savePredictions = TRUE,
  classProbs      = TRUE # Required for ROC and probability analysis
)

trainedModels <- list()

for (idx in seq_len(nrow(model_settings))) {
  algoName <- model_settings$AlgorithmName[idx]
  algoImpl <- model_settings$Implementation[idx]
  
  message(paste0("Training: ", algoName, "..."))
  
  if (algoImpl == "svmRadial") {
    model <- train(Type ~ ., data = train, method = algoImpl,
                   prob.model = TRUE, trControl = cv_control)
  } else if (algoImpl == "nnet") {
    model <- train(Type ~ ., data = train, method = algoImpl,
                   trControl = cv_control, trace = FALSE)
  } else {
    model <- train(Type ~ ., data = train, method = algoImpl,
                   trControl = cv_control)
  }
  
  trainedModels[[algoName]] <- model
}

# ==============================================================================
# 4. Model Evaluation (DALEX)
# ==============================================================================

# Prediction wrapper for DALEX
p_fun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, 2]
}

# Binary outcome for validation (Assumes "Control" is the baseline)
yTest <- ifelse(test$Type == "Control", 0, 1)
explainers  <- list()
model_perfs <- list()

for (algoName in names(trainedModels)) {
  explainer <- explain(
    trainedModels[[algoName]],
    label            = algoName,
    data             = test[, -which(names(test) == "Type")],
    y                = yTest,
    predict_function = p_fun,
    verbose          = FALSE
  )
  explainers[[algoName]]  <- explainer
  model_perfs[[algoName]] <- model_performance(explainer)
}

# 4.1 Residual Analysis & Boxplots
pdf("Results/Model_Performance_Comparison.pdf", width = 12, height = 6)
p1 <- plot(model_perfs[[1]], model_perfs[[2]], model_perfs[[3]], model_perfs[[4]])
p2 <- plot(model_perfs[[1]], model_perfs[[2]], model_perfs[[3]], model_perfs[[4]], geom = "boxplot")
grid.arrange(p1, p2, ncol = 2)
dev.off()

# 4.2 ROC Curve Comparison
pdf("Results/ROC_Comparison.pdf", width = 6, height = 6)
colors <- c("chocolate", "aquamarine3", "darkolivegreen3", "dodgerblue3")
roc_curves <- list()
auc_labels <- c()

for (i in seq_along(trainedModels)) {
  algoName <- names(trainedModels)[i]
  pred     <- predict(trainedModels[[algoName]], newdata = test, type = "prob")
  roc_obj  <- roc(yTest, as.numeric(pred[, 2]))
  roc_curves[[algoName]] <- roc_obj
  auc_labels <- c(auc_labels, paste0(algoName, " (AUC: ", sprintf("%.3f", roc_obj$auc), ")"))
  
  plot(roc_obj, add = (i > 1), col = colors[i], legacy.axes = TRUE, 
       lwd = 2, main = "ROC Curve Comparison")
}
abline(a = 0, b = 1, lty = 3, col = "gray")
legend("bottomright", legend = auc_labels, col = colors, lwd = 2, bty = "n")
dev.off()

# ==============================================================================
# 5. Global Feature Importance (DALEX)
# ==============================================================================

pdf("Results/DALEX_Feature_Importance.pdf", width = 10, height = 10)
importance_plots <- list()

for (algoName in names(trainedModels)) {
  imp <- variable_importance(explainers[[algoName]], loss_function = loss_root_mean_square)
  
  # Export Top 12 Hub Genes
  write.table(head(imp[order(imp$dropout_loss, decreasing = TRUE), ], 12), 
              file = paste0("Results/Top_Genes_", gsub(" ", "_", algoName), ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
# (Note: Use standard plot() or ggplot to visualize the 'imp' objects here)
dev.off()

# ==============================================================================
# 6. SHAP Analysis (Best Model)
# ==============================================================================

# Select best model based on AUC
auc_values      <- sapply(roc_curves, function(x) x$auc)
best_model_name <- names(which.max(auc_values))
best_model      <- trainedModels[[best_model_name]]

message(paste("Best Model Selected:", best_model_name))

# Prepare SHAP values (Using KernelSHAP)
X_train <- train[, -which(names(train) == "Type")]
shap_values <- kernelshap(
  best_model,
  X        = X_train[1:min(70, nrow(X_train)), ], 
  bg_X     = X_train,
  pred_fun = function(obj, newdata) as.numeric(predict(obj, newdata = newdata, type = "prob")[, 2])
)
shap_vis <- shapviz(shap_values)

# 6.1 SHAP Feature Importance
pdf("Results/SHAP_Feature_Importance.pdf", width = 8, height = 6)
imp_data <- enframe(colMeans(abs(shap_vis$S)), name = "Feature", value = "Importance") %>%
  arrange(desc(Importance)) %>% head(15)

ggplot(imp_data, aes(x = reorder(Feature, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C", name = "Importance") +
  coord_flip() + theme_minimal() +
  labs(title = paste("Global SHAP Importance (", best_model_name, ")"),
       x = "Features", y = "Mean |SHAP Value|")
dev.off()

# 6.2 SHAP Beeswarm Plot
pdf("Results/SHAP_Beeswarm.pdf", width = 9, height = 7)
sv_importance(shap_vis, kind = "bee") + 
  theme_minimal() +
  labs(title = paste("SHAP Value Distribution:", best_model_name))
dev.off()

# 6.3 Local Explanations (Waterfall and Force Plot)
pdf("Results/SHAP_Individual_Explanations.pdf", width = 12, height = 8)
p_waterfall <- sv_waterfall(shap_vis, row_id = 1) + labs(title = "Waterfall Plot (Sample 1)")
p_force     <- sv_force(shap_vis, row_id = 1) + labs(title = "Force Plot (Sample 1)")
grid.arrange(p_waterfall, p_force, ncol = 2)
dev.off()