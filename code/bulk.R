rm(list = ls())

library(dplyr)  
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(tibble)



# 0,1,
seeds <- c(42)

for (seed in seeds) {
  set.seed(seed)
  
  load('/work/DongZ/project/CosMx/CRC/data/bulk/CRC.Rdata')
  genes <- readr::read_rds('/work/DongZ/project/CosMx/CRC/data/bulk/LR_genes.rds')
  
  
  bulk <- t(data)
  bulk <- as.data.frame(bulk[,intersect(genes,colnames(bulk))]) %>% rownames_to_column()
  colnames(bulk)[1] <- 'barcode'
  
  merged_data <- merge(bulk, MMR, by = "barcode")
  merged_data$subtype <- gsub("MSI-H|MSI-L", "MSI", merged_data$subtype)
  table(merged_data$subtype)
  

  X <- merged_data[, intersect(genes,colnames(bulk))]
  # 标签
  y <- as.factor(merged_data$subtype)  

  folds <- createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)

  results_list <- list()
  sample_pred_list <- list()  # 用于存储样本级预测结果
  

  rf_grid <- expand.grid(mtry = seq(1, ncol(X), by = 1))
  xgb_grid <- expand.grid(
    nrounds = c(50, 100),    
    max_depth = c(3, 5),      
    eta = c(0.01,0.05,0.1),     
    gamma = c(0, 0.1),         
    colsample_bytree = 0.8,   
    min_child_weight = 1,    
    subsample = 0.5         
  )
  glmnet_grid <- expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 2, length = 8))
  lasso_grid <- expand.grid(alpha = 1,lambda = 10^seq(-3, 2, length = 10))

  svm_grid <- expand.grid(C = 10^seq(-2, 2, length = 5)) # SVM线性核
  knn_grid <- expand.grid(k = seq(1, 21, by = 2)) # kNN (奇数k值)
  gbm_grid <- expand.grid(
    n.trees = c(100, 200, 500),       
    interaction.depth = c(2, 3, 5),  
    shrinkage = c(0.01, 0.05, 0.1),  
    n.minobsinnode = c(1, 3, 5)        
  )
  
  # 循环处理每个fold
  for (i in seq_along(folds)) {
    cat("\n\nProcessing Fold", i, "/", length(folds), "\n")
    
    validIndex <- folds[[i]]
    X_valid_raw <- X[validIndex, ]
    X_train_raw <- X[-validIndex, ]
    y_valid <- y[validIndex]
    y_train <- y[-validIndex]
    
    preProc <- preProcess(X_train_raw, method = c("center", "scale"))
    X_train_proc <- predict(preProc, X_train_raw)
    X_valid_proc <- predict(preProc, X_valid_raw)
    
    ctrl <- trainControl(
      method = "repeatedcv",
      number = 5,      # K值，通常取 5 或 10
      repeats = 3,     # 重复次数，建议 3 到 10 次
      classProbs = TRUE,
      verboseIter = TRUE,
      summaryFunction = twoClassSummary # 必须使用，因为是二分类问题
    )
    
    models <- list(
      rf = caret::train(X_train_proc, y_train,
                        method = "rf",
                        trControl = ctrl,
                        tuneGrid = rf_grid,
                        metric = "ROC"),
      xgb = caret::train(X_train_proc, y_train,
                         method = "xgbTree",
                         trControl = ctrl,
                         tuneGrid = xgb_grid,
                         metric = "ROC",
                         nthread = 50),
      glm = caret::train(X_train_proc, y_train,
                        method = "glm",
                        trControl = ctrl,
                        metric = "ROC",
                        family = binomial()),  # 明确指定二项分布（二分类）
      lasso = caret::train(X_train_proc, y_train,
                          method = "glmnet",
                          trControl = ctrl,
                          tuneGrid = lasso_grid,  # 使用专属LASSO网格
                          metric = "ROC"),
      svm = caret::train(X_train_proc, y_train,
                         method = "svmLinear",
                         trControl = ctrl,
                         tuneGrid = svm_grid,
                         metric = "ROC",
                         probability = TRUE),
      knn = caret::train(X_train_proc, y_train,
                         method = "knn",
                         trControl = ctrl,
                         tuneGrid = knn_grid,
                         metric = "ROC"),
      gbm = caret::train(X_train_proc, y_train,
                         method = "gbm",
                         trControl = ctrl,
                         tuneGrid = gbm_grid,
                         metric = "ROC",
                         verbose = FALSE)
    )

    fold_results <- list()
    
    for (model_name in names(models)) {
      model <- models[[model_name]]

      pred_prob <- predict(model, X_valid_proc, type = "prob")
      pred_class <- predict(model, X_valid_proc)

      cm <- confusionMatrix(pred_class, y_valid, positive = "MSS")

      roc_obj <- roc(
        response = y_valid,
        predictor = pred_prob[,"MSS"],
        levels = c("MSI", "MSS"),
        direction = "<"   # MSS 类别将被视为正类，而 MSI 类别将被视为负类
      )
      
      metrics <- data.frame(
        Fold = i,
        Model = model_name,
        ACC = as.numeric(cm$overall["Accuracy"]),
        Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
        Specificity = as.numeric(cm$byClass["Specificity"]),
        F1 = as.numeric(cm$byClass["F1"]),
        AUC = as.numeric(auc(roc_obj)),  # 显式转换为数值
        stringsAsFactors = FALSE
      )
      
      fold_results[[model_name]] <- metrics
 
      sample_pred_df <- data.frame(
        Fold = i,
        Model = model_name,
        Sample = rownames(X_valid_proc),
        TrueLabel = as.character(y_valid),
        PredProb = pred_prob[, "MSS"],
        stringsAsFactors = FALSE
      )
      
      sample_pred_list[[length(sample_pred_list) + 1]] <- sample_pred_df
    }

    results_list[[i]] <- bind_rows(fold_results)
  }

  final_results <- bind_rows(results_list)
  sample_pred_all <- bind_rows(sample_pred_list)

  summary_stats <- final_results %>%
    group_by(Model) %>%
    summarise(
      Mean_ACC = mean(ACC),
      Mean_F1 = mean(F1),
      Mean_AUC = mean(AUC),
      .groups = 'drop'
    )
  
  write.csv(final_results, paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_bulk/result/cv_metrics_detail_", seed, ".csv"), row.names = FALSE)
  write.csv(summary_stats, paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_bulk/result/cv_metrics_summary_", seed, ".csv"), row.names = FALSE)
  write.csv(sample_pred_all, paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_bulk/result/sample_predictions_all_folds_", seed, ".csv"), row.names = FALSE)

  model_names <- unique(sample_pred_all$Model)
  roc_combined <- list()
  
  for (model in model_names) {
    model_data <- sample_pred_all %>% filter(Model == model)
    
    roc_combined[[model]] <- roc(
      response = model_data$TrueLabel,
      predictor = model_data$PredProb,
      levels = c("MSI", "MSS"),
      direction = "<"
    )
  }

  mean_auc <- sapply(roc_combined, auc)
  
  plot_data <- data.frame()
  for (model in names(roc_combined)) {
    roc_obj <- roc_combined[[model]]
    model_data <- data.frame(
      Model = model,
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, model_data)
  }
  
  
  my_palette <- c(
    "glm" = "#E58F8E",   # 红色
    "xgb" = "#ecb46c",   # 橙色
    "rf" = "#559d4f",    # 绿色
    "svm" = "#6A5ACD",   # 紫罗兰色 (SlateBlue)
    "knn" = "#20B2AA",   # 青色 (LightSeaGreen)
    "lasso" = "#9370DB",    # 紫色 (MediumPurple)
    "gbm" = "#FF6347"    # 番茄红 (Tomato)
  )
  
  pdf(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_bulk/ROC_", seed, ".pdf"), 
      width = 8, height = 7)
  
  plot.roc(roc_combined$glm, 
           col = my_palette["glm"], 
           ylim = c(0, 1), 
           legacy.axes = TRUE,
           print.auc = TRUE, 
           print.auc.y = 0.7,
           auc.polygon.col = paste0(my_palette["glm"], "20"),
           lwd = 2)
  plot.roc(roc_combined$xgb, 
           add = TRUE, 
           col = my_palette["xgb"],
           print.auc = TRUE, 
           print.auc.y = 0.6,
           auc.polygon.col = paste0(my_palette["xgb"], "20"),
           lwd = 2)
  plot.roc(roc_combined$rf, 
           add = TRUE, 
           col = my_palette["rf"],
           print.auc = TRUE, 
           print.auc.y = 0.5,
           auc.polygon.col = paste0(my_palette["rf"], "20"),
           lwd = 2)
  plot.roc(roc_combined$svm, 
           add = TRUE, 
           col = my_palette["svm"],
           print.auc = TRUE, 
           print.auc.y = 0.4,
           auc.polygon.col = paste0(my_palette["svm"], "20"),
           lwd = 2)
  plot.roc(roc_combined$knn, 
           add = TRUE, 
           col = my_palette["knn"],
           print.auc = TRUE, 
           print.auc.y = 0.3,
           auc.polygon.col = paste0(my_palette["knn"], "20"),
           lwd = 2)
  plot.roc(roc_combined$lasso, 
           add = TRUE, 
           col = my_palette["lasso"],
           print.auc = TRUE, 
           print.auc.y = 0.2,
           auc.polygon.col = paste0(my_palette["lasso"], "20"),
           lwd = 2)
  plot.roc(roc_combined$gbm, 
           add = TRUE, 
           col = my_palette["gbm"],
           print.auc = TRUE, 
           print.auc.y = 0.1,
           auc.polygon.col = paste0(my_palette["gbm"], "20"),
           lwd = 2)
  legend("bottomright", 
         legend = c(paste0("glm (AUC = ", round(mean_auc["glm"], 3), ")"),
                    paste0("xgb (AUC = ", round(mean_auc["xgb"], 3), ")"),
                    paste0("rf (AUC = ", round(mean_auc["rf"], 3), ")"),
                    paste0("svm (AUC = ", round(mean_auc["svm"], 3), ")"),
                    paste0("knn (AUC = ", round(mean_auc["knn"], 3), ")"),
                    paste0("lasso (AUC = ", round(mean_auc["lasso"], 3), ")"),
                    paste0("gbm (AUC = ", round(mean_auc["gbm"], 3), ")")),
         col = my_palette, 
         lwd = 2, 
         cex = 0.9,
         bty = "n")

  dev.off()
  
  
  model_colors <- c(
    "glm" = "#E58F8E",   # 红色
    "xgb" = "#ecb46c",   # 橙色
    "rf" = "#559d4f",    # 绿色
    "svm" = "#6A5ACD",   # 紫罗兰色 (SlateBlue)
    "knn" = "#20B2AA",   # 青色 (LightSeaGreen)
    "lasso" = "#9370DB",    # 紫色 (MediumPurple)
    "gbm" = "#FF6347"    # 番茄红 (Tomato)
  )
  
  # 绘制AUC结果
  auc_plot <- ggplot(final_results, aes(x = Model, y = AUC, fill = Model)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1, size = 2) +
    scale_fill_manual(values = model_colors) +
    labs(title = "Model Performance Across 5 Folds",
         subtitle = "AUC Distribution",
         y = "AUC") +
    theme_minimal()
  
  # 绘制ACC结果
  acc_plot <- ggplot(final_results, aes(x = Model, y = ACC, fill = Model)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1, size = 2) +
    scale_fill_manual(values = model_colors) +
    labs(y = "Accuracy") +
    theme_minimal()
  
  # 绘制F1结果
  f1_plot <- ggplot(final_results, aes(x = Model, y = F1, fill = Model)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1, size = 2) +
    scale_fill_manual(values = model_colors) +
    labs(y = "F1 Score") +
    theme_minimal()
  
  # 保存图表
  ggsave(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_bulk/AUC_", seed, ".pdf"), auc_plot, width = 6, height = 5)
  ggsave(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_bulk/ACC_", seed, ".pdf"), acc_plot, width = 6, height = 5)
  ggsave(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_bulk/F1_", seed, ".pdf"), f1_plot, width = 6, height = 5)
  
}




