# nohup Rscript pseudobulk_sc.R > pseudobulk_sc.log 2>&1 &

library(dplyr)  
library(magrittr)
library(rlang)
library(pROC)
library(dplyr)  
library(randomForest)
library(caret)
library(tibble)

#### sc_pseudobulk ####

meta <- readr::read_rds('/work/xyding/2021/2021CRC/scRNAseq/processdata/Merge/MergeMetaDataNewMainTypes20241016.rds.gz')
meta <- meta[meta$orig.ident != 'OurData', ]

meta <- meta[!is.na(meta$MMR_Status), ]
meta <- meta[meta$Tissue == 'Tumor', ]
meta <- meta[,c('Patient', 'MMR_Status', 'Sample')]
table(meta$Sample, meta$MMR_Status)


# sc_1 <- readRDS("/work/DongZ/project/CosMx/CRC/data/sc_ref.rds")
sc <- readRDS("/work/xyding/2021/2021CRC/scRNAseq/processdata/Merge/WholeTissueList.rds.gz")

list1 <- colnames(sc$TCells)
list2 <- colnames(sc$BCells)
list3 <- colnames(sc$MyeloidCells)
list4 <- colnames(sc$Fibroblasts)
list5 <- colnames(sc$EndothelialCells)
list6 <- colnames(sc$EpithelialCellsT)
list7 <- colnames(sc$EpithelialCellsN)

list_all <- Reduce(union, list(list1, list2, list3, list4, list5, list6, list7))
list_unique <- unique(list_all)

Wholedata <- merge(x=sc[[1]],y=sc[2:length(sc)],merge.data=T)

patients <- unique(meta$Sample)
patients_whole <- unique(Wholedata$Sample)
all(patients%in%patients_whole)

Wholedata <- subset(Wholedata, Sample %in% patients)
Wholedata <- subset(Wholedata, Tissue == 'Tumor')

Wholedata$MMR <- meta$MMR_Status[match(Wholedata$Patient, meta$Patient)]
table(meta$Sample, meta$MMR_Status)
table(Wholedata$Sample, Wholedata$MMR)


### 生成拟bulk 数据 ###
bulk_data <- AggregateExpression(Wholedata, return.seurat = T, slot = "counts", assays = "RNA",
                                 group.by = c("Sample", "MMR")# 分别填写细胞类型、样本变量、分组变量的slot名称
)
bulk_data$MMR <- rownames(bulk_data@meta.data) %>%sub(".*_", "", .)


saveRDS(bulk_data, file = '/work/DongZ/project/CosMx/CRC/data/bulk/pseudobulk_sc_new.rds')

# bulk_data <- readRDS('/work/DongZ/project/CosMx/CRC/data/bulk/pseudobulk_sc.rds')
# 
# 
# 
# MMR <- as.data.frame(bulk_data$MMR) %>% set_names('subtype')  %>% rownames_to_column()
# 
# genes <- readr::read_rds('/work/DongZ/project/CosMx/CRC/data/bulk/LR_genes.rds')
# 
# data <- as.matrix(bulk_data@assays$RNA@counts)
# bulk <- t(data)
# bulk <- as.data.frame(bulk[,intersect(genes,colnames(bulk))]) %>% rownames_to_column()
# 
# 
# merged_data <- merge(bulk, MMR, by = "rowname")
# table(merged_data$subtype)
# 
# 
# # 基因表达
# X <- merged_data[, intersect(genes,colnames(bulk))]
# # 标签
# y <- as.factor(merged_data$subtype)  



# 0,1,

seeds <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,42,123,1234,520,1314,2025,777,999,
           100,200,300,400,500,600,700,800,900,1000,666)

seeds <-c(42)
for (seed in seeds) {
  set.seed(seed)
  
  bulk_data <- readRDS('/work/DongZ/project/CosMx/CRC/data/bulk/pseudobulk_sc_new.rds')
  
  MMR <- as.data.frame(bulk_data$MMR) %>% set_names('subtype')  %>% rownames_to_column()
  
  genes <- readr::read_rds('/work/DongZ/project/CosMx/CRC/data/bulk/LR_genes.rds')
  
  data <- as.matrix(bulk_data@assays$RNA@counts)
  bulk <- t(data)
  bulk <- as.data.frame(bulk[,intersect(genes,colnames(bulk))]) %>% rownames_to_column()
  
  
  merged_data <- merge(bulk, MMR, by = "rowname")
  table(merged_data$subtype)
  
  
  # 基因表达
  X <- merged_data[, intersect(genes,colnames(bulk))]
  # 标签
  y <- as.factor(merged_data$subtype)  
  
  
  # 创建分层5折交叉验证索引
  folds <- createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)
  
  # 初始化结果收集器
  results_list <- list()
  sample_pred_list <- list()  # 用于存储样本级预测结果
  
  # 定义模型调参网格
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
  # 新增模型网格
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
    
    # 划分训练集和验证集
    validIndex <- folds[[i]]
    X_valid_raw <- X[validIndex, ]
    X_train_raw <- X[-validIndex, ]
    y_valid <- y[validIndex]
    y_train <- y[-validIndex]
    
    # 预处理（在分完fold之后再标准化，防止数据泄露）
    preProc <- preProcess(X_train_raw, method = c("center", "scale"))
    X_train_proc <- predict(preProc, X_train_raw)
    X_valid_proc <- predict(preProc, X_valid_raw)
    
    
    # 训练控制设置 - 使用ROC作为主要指标
    # ctrl <- trainControl(
    #   method = "cv",
    #   number = 5,
    #   verboseIter = TRUE,
    #   classProbs = TRUE,
    #   summaryFunction = twoClassSummary
    # )
    ctrl <- trainControl(
      method = "repeatedcv",
      number = 5,      # K值，通常取 5 或 10
      repeats = 3,     # 重复次数，建议 3 到 10 次
      classProbs = TRUE,
      verboseIter = TRUE,
      summaryFunction = twoClassSummary # 必须使用，因为是二分类问题
    )
    
    # 训练模型 - 所有模型都使用metric = "ROC"
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
      # glm = caret::train(X_train_proc, y_train,
      #                    method = "glmnet",
      #                    trControl = ctrl,
      #                    tuneGrid = glmnet_grid,
      #                    metric = "ROC"),
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
    
    # 在验证集上预测并收集结果
    fold_results <- list()
    
    for (model_name in names(models)) {
      model <- models[[model_name]]
      
      # 获取预测概率和类别
      pred_prob <- predict(model, X_valid_proc, type = "prob")
      pred_class <- predict(model, X_valid_proc)
      
      # 计算性能指标
      cm <- confusionMatrix(pred_class, y_valid, positive = "MSS")
      
      # 修复：显式指定ROC参数
      roc_obj <- roc(
        response = y_valid,
        predictor = pred_prob[,"MSS"],
        levels = c("MSI", "MSS"),
        direction = "<"   # MSS 类别将被视为正类，而 MSI 类别将被视为负类
      )
      
      # 收集指标 - 确保AUC存储为数值
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
      
      # 收集样本级预测结果用于合并ROC曲线
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
    
    # 合并当前fold的结果
    results_list[[i]] <- bind_rows(fold_results)
  }
  
  # 合并所有fold的结果
  final_results <- bind_rows(results_list)
  sample_pred_all <- bind_rows(sample_pred_list)
  
  # 计算汇总统计
  summary_stats <- final_results %>%
    group_by(Model) %>%
    summarise(
      Mean_ACC = mean(ACC),
      Mean_F1 = mean(F1),
      Mean_AUC = mean(AUC),
      .groups = 'drop'
    )
  
  write.csv(final_results, paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_sc/result/cv_metrics_detail_", seed, ".csv"), row.names = FALSE)
  write.csv(summary_stats, paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_sc/result/cv_metrics_summary_", seed, ".csv"), row.names = FALSE)
  write.csv(sample_pred_all, paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_sc/result/sample_predictions_all_folds_", seed, ".csv"), row.names = FALSE)
  
  
  # 计算每个模型的合并ROC曲线
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
  
  # 计算平均AUC
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
  
  pdf(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_sc/ROC_", seed, ".pdf"), 
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
  
  # 关闭图形设备
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
  ggsave(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_sc/AUC_", seed, ".pdf"), auc_plot, width = 6, height = 5)
  ggsave(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_sc/ACC_", seed, ".pdf"), acc_plot, width = 6, height = 5)
  ggsave(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure6/ROC_sc/F1_", seed, ".pdf"), f1_plot, width = 6, height = 5)
  
  
  
}









# 找到最佳的seed重新绘图
seed <- 999

final_results <- read_csv(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure5/ROC_sc/result/cv_metrics_detail_", seed, ".csv"))
summary_stats <- read_csv(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure5/ROC_sc/result/cv_metrics_summary_", seed, ".csv"))
sample_pred_all <- read_csv(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure5/ROC_sc/result/sample_predictions_all_folds_", seed, ".csv"))

# 计算每个模型的合并ROC曲线
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

pdf("/work/DongZ/project/CosMx/CRC/figures_final/figure5/ROC_sc.pdf",
    width = 5, height = 5)
plot.roc(roc_combined$glm, 
         col = my_palette["glm"], 
         ylim = c(0, 1), 
         legacy.axes = TRUE,
         print.auc = TRUE, 
         print.auc.y = 0.7,
         print.auc.x = 0.3,
         auc.polygon.col = paste0(my_palette["glm"], "20"),
         lwd = 2)
plot.roc(roc_combined$xgb, 
         add = TRUE, 
         col = my_palette["xgb"],
         print.auc = TRUE, 
         print.auc.y = 0.6,
         print.auc.x = 0.3,
         auc.polygon.col = paste0(my_palette["xgb"], "20"),
         lwd = 2)
plot.roc(roc_combined$rf, 
         add = TRUE, 
         col = my_palette["rf"],
         print.auc = TRUE, 
         print.auc.y = 0.5,
         print.auc.x = 0.3,
         auc.polygon.col = paste0(my_palette["rf"], "20"),
         lwd = 2)
plot.roc(roc_combined$svm, 
         add = TRUE, 
         col = my_palette["svm"],
         print.auc = TRUE, 
         print.auc.y = 0.4,
         print.auc.x = 0.3,
         auc.polygon.col = paste0(my_palette["svm"], "20"),
         lwd = 2)
plot.roc(roc_combined$knn, 
         add = TRUE, 
         col = my_palette["knn"],
         print.auc = TRUE, 
         print.auc.y = 0.3,
         print.auc.x = 0.3,
         auc.polygon.col = paste0(my_palette["knn"], "20"),
         lwd = 2)
plot.roc(roc_combined$lasso, 
         add = TRUE, 
         col = my_palette["lasso"],
         print.auc = TRUE, 
         print.auc.y = 0.2,
         print.auc.x = 0.3,
         auc.polygon.col = paste0(my_palette["lasso"], "20"),
         lwd = 2)
plot.roc(roc_combined$gbm, 
         add = TRUE,
         col = my_palette["gbm"],
         print.auc = TRUE, 
         print.auc.y = 0.1,
         print.auc.x = 0.3,
         auc.polygon.col = paste0(my_palette["gbm"], "20"),
         lwd = 2)
dev.off()



# boxplot
final_long <- final_results %>%
  pivot_longer(cols = c(AUC, ACC, F1),
               names_to = "Metric",
               values_to = "Value")

model_colors <- c(
  "glm" = "#E58F8E",   # 红色
  "xgb" = "#ecb46c",   # 橙色
  "rf" = "#559d4f",    # 绿色
  "svm" = "#6A5ACD",   # 紫罗兰色 (SlateBlue)
  "knn" = "#20B2AA",   # 青色 (LightSeaGreen)
  "lasso" = "#9370DB",    # 紫色 (MediumPurple)
  "gbm" = "#FF6347"    # 番茄红 (Tomato)
)

ggplot(final_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8),outlier.shape = NA) +
  scale_fill_manual(values = model_colors) +
  labs(
    x = "Metric",
    y = "Score"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black"),
    legend.position = "top",
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )
ggsave("/work/DongZ/project/CosMx/CRC/figures_final/figure5/Metric_sc.pdf",width = 6, height = 4)













# 
# 
# seeds <- seq(420, 1000, by = 10)
# for (seed in seeds) {
#   set.seed(seed)
#   
#   # 数据划分
#   trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
#   X_train <- X[trainIndex, ]
#   y_train <- y[trainIndex]
#   X_test <- X[-trainIndex, ]
#   y_test <- y[-trainIndex]
#   
#   # 预处理
#   preProc <- preProcess(X_train, method = c("center", "scale"))
#   X_train <- predict(preProc, X_train)
#   X_test <- predict(preProc, X_test)
#   
#   # 训练控制
#   ctrl <- trainControl(
#     method = "cv",
#     number = 5,
#     verboseIter = TRUE,
#     classProbs = TRUE,
#     summaryFunction = twoClassSummary
#   )
#   
#   # 调参网格
#   rf_grid <- expand.grid(mtry = seq(1, ncol(X_train), by = 2))
#   xgb_grid <- expand.grid(
#     nrounds = c(50, 100, 150),
#     max_depth = c(1, 3, 5),
#     eta = c(0.01, 0.05),
#     gamma = c(0, 0.1, 1),
#     colsample_bytree = 0.8,
#     min_child_weight = 1,
#     subsample = 0.8
#   )
#   glmnet_grid <- expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-3, 2, length = 8))
#   
#   # 训练模型
#   models <- list(
#     rf = caret::train(X_train, y_train,
#                       method = "rf",
#                       trControl = ctrl,
#                       tuneGrid = rf_grid,
#                       metric = "ROC",
#                       importance = TRUE),
#     xgb = caret::train(X_train, y_train,
#                        method = "xgbTree",
#                        trControl = ctrl,
#                        tuneGrid = xgb_grid),
#     glm = caret::train(X_train, y_train,
#                        method = "glmnet",
#                        trControl = ctrl,
#                        tuneGrid = glmnet_grid)
#   )
#   
#   # 获取测试集预测概率
#   pred_prob_list <- list(
#     rf = predict(models$rf, X_test, type = "prob"),
#     xgb = predict(models$xgb, X_test, type = "prob"),
#     glm = predict(models$glm, X_test, type = "prob")
#   )
#   
#   # 计算测试集 ROC
#   roc_test <- lapply(pred_prob_list, function(pred_prob) {
#     roc(response = y_test, predictor = pred_prob[,2])
#   })
#   auc_test <- sapply(roc_test, auc)
#   
#   # 计算训练集预测概率
#   pred_train_prob <- list(
#     rf = predict(models$rf, X_train, type = "prob"),
#     xgb = predict(models$xgb, X_train, type = "prob"),
#     glm = predict(models$glm, X_train, type = "prob")
#   )
#   
#   # 计算训练集 ROC
#   roc_train <- lapply(pred_train_prob, function(pred_prob) {
#     roc(response = y_train, predictor = pred_prob[,2])
#   })
#   auc_train <- sapply(roc_train, auc)
#   
#   # 比较训练集和测试集 AUC
#   diff <- auc_train - auc_test
#   overfit <- any(diff > 0.2)  # 如果任一模型差异 > 0.2，认为过拟合
#   
#   # 输出结果
#   cat("Seed:", seed, "\n")
#   cat("Model | Train AUC | Test AUC | Diff\n")
#   cat("------|-----------|----------|------\n")
#   for (i in names(auc_train)) {
#     cat(i, "|", round(auc_train[i], 3), "|", round(auc_test[i], 3), "|", round(diff[i], 3), "\n")
#   }
#   cat("\n")
#   
#   # 如果没有过拟合，才绘制 ROC 图
#   if (!overfit) {
#     my_palette <- c("#E58F8E", "#ecb46c", "#559d4f", "#9cc4d4")
#     pdf(paste0("/work/DongZ/project/CosMx/CRC/figures_final/figure5/ROC_sc/ROC_", seed, ".pdf"), width = 4, height = 4)
#     plot.roc(roc_test$glm, col=my_palette[1], ylim=c(0,1), legacy.axes=TRUE, 
#              print.auc=TRUE, print.auc.y=0.4,
#              auc.polygon.col="#FF000020")
#     plot.roc(roc_test$xgb, add=TRUE, col=my_palette[2], 
#              print.auc=TRUE, print.auc.y=0.3)
#     plot.roc(roc_test$rf, add=TRUE, col=my_palette[3], 
#              print.auc=TRUE, print.auc.y=0.2)
#     dev.off()
#   } else {
#     cat("过拟合检测：AUC_train - AUC_test > 0.2，跳过 ROC 图绘制。\n")
#   }
# }

