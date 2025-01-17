#在用此代码作图前先调用 funky_heatmap_modify_func.R，修改其中一部分作图参数。

library(funkyheatmap)
library(dplyr)
library(tibble)
library(tidyr)

setwd('C:/Users/Li Shen/Desktop/scDrugPredict_benchmark/just4test/')
data <- read.csv('funkyheatmap.csv', header = T)


sciduc_data <- data %>%
  filter(Tool %in% c("sciduc_GDSC", "sciduc_CCLE")) %>%
  group_by(Dataset) %>%
  slice_max(AUC, with_ties = FALSE) %>%
  mutate(Tool = "sciduc")

dreep_data <- data %>%
  filter(Tool == "DREEP") %>%
  group_by(Dataset) %>%
  slice_max(AUC, with_ties = FALSE)

beyondcell_data <- data %>%
  filter(Tool == "beyondcell_PSc") %>%
  group_by(Dataset) %>%
  slice_max(AUC, with_ties = FALSE)

other_data <- data %>%
  filter(!Tool %in% c("sciduc_GDSC", "sciduc_CCLE", "DREEP", "beyondcell_PSc"))

final_data <- bind_rows(sciduc_data, dreep_data, beyondcell_data, other_data)

final_data <- final_data %>%
  mutate(Drug = ifelse(Tool == "DREEP", 
                       sub("^.*_(.*?):.*$", "\\1", Drug),  # 提取 _ 和 : 之间的药物名字
                       Drug))

final_data <- final_data %>%
  mutate(Drug = ifelse(Tool == "beyondcell_PSc", 
                       gsub("\\s*\\(.*\\)$", "", Drug),  # 去掉括号及括号中的内容
                       Drug))

final_data <- final_data %>%
  mutate(Drug = case_when(
    Drug == "NVP-TAE684" ~ "TAE684",
    TRUE ~ Drug  # 保留其他值不变
  )) #替换NVP-TAE684为TAE684


data_long <- final_data %>%
  pivot_longer(cols = c(AUC, ACC), names_to = "Metric", values_to = "Value") %>%
  unite("Tool_Metric", Tool, Metric, sep = "_") %>%
  pivot_wider(names_from = Tool_Metric, values_from = Value)

#将AUC和ACC分别聚合
data_long <- data_long %>%
  select(Dataset, Drug, 
         ends_with("_AUC"), 
         ends_with("_ACC"))


auc_cols <- grep("_AUC$", colnames(data_long), value = TRUE)
acc_cols <- grep("_ACC$", colnames(data_long), value = TRUE)

# 自定义排名函数：值为 NA 或 0 时不参与排名
custom_rank <- function(x) {
  x[is.na(x) | x == 0] <- NA  # 将 NA 或 0 保留为 NA
  ranks <- rank(-x, na.last = "keep", ties.method = "min")  # 按从大到小排序
  ranks[ranks > 3] <- NA  # 排名超过 3 的值设置为 NA
  return(ranks)
}


insert_rank_columns_by_row <- function(df, suffix = "", rank_suffix = "") {
  col_names <- grep(suffix, colnames(df), value = TRUE)  # 找到目标列
  
  # 对每一行计算排名
  rank_matrix <- t(apply(df[col_names], 1, custom_rank))  # 每行排名
  
  # 动态插入排名列
  for (col in rev(col_names)) {  # 从后往前插入，避免索引混乱
    rank_col <- paste0(rank_suffix, col)  # 构造排名列名
    col_index <- which(colnames(df) == col)  # 找到当前列索引
    
    # 插入排名列前后部分
    if (col_index < ncol(df)) {
      df <- cbind(
        df[, 1:col_index, drop = FALSE],  # 插入前部分
        setNames(as.data.frame(rank_matrix[, which(col_names == col), drop = FALSE]), rank_col),  # 插入排名列
        df[, (col_index + 1):ncol(df), drop = FALSE]  # 插入后部分
      )
    } else {
      df <- cbind(
        df[, 1:col_index, drop = FALSE],  # 插入前部分
        setNames(as.data.frame(rank_matrix[, which(col_names == col), drop = FALSE]), rank_col)  # 插入排名列
      )
    }
  }
  return(df)
}

data_long <- insert_rank_columns_by_row(data_long, suffix = "_AUC", rank_suffix = "ranked_")
data_long <- insert_rank_columns_by_row(data_long, suffix = "_ACC", rank_suffix = "ranked_")


#data_long[paste0("ranked_", auc_cols)] <- t(apply(data_long[auc_cols], 1, custom_rank))
#data_long[paste0("ranked_", acc_cols)] <- t(apply(data_long[acc_cols], 1, custom_rank))


# 将 AUC和ACC中的NA 替换为 0
#data_long[, 3:10][is.na(data_long[, 3:10])] <- 0
fill_na_conditionally <- function(df, start_col, end_col, exclude_prefix = "ranked") {
  # 筛选目标列名：列范围在 start_col 到 end_col 且不以 exclude_prefix 开头
  target_cols <- colnames(df)[start_col:end_col]
  target_cols <- target_cols[!grepl(paste0("^", exclude_prefix), target_cols)]
  
  # 填充 NA 为 0
  df[target_cols] <- lapply(df[target_cols], function(col) {
    ifelse(is.na(col), 0, col)
  })
  
  return(df)
}

data_long <- fill_na_conditionally(data_long, start_col = 3, end_col = 18, exclude_prefix = "ranked")

#一定要取消group
data_long <- data_long %>% ungroup()


# 创建 column_info
column_info <- data.frame(
  id = names(data_long),
  group = c("Dataset", "Drug",
            rep(c("AUC",""), length.out = 8),
            rep(c("ACC",""), length.out = 8)
            ),
  name = c('','',
           'ScIDUC','','DREEP','',
           'beyondcell','','scDr','',
           'ScIDUC','','DREEP','',
           'beyondcell','','scDr',''
          ),
  geom = c("text", "text",
           rep(c("bar","text"), length.out = 16)
           ), 
  palette = c(NA, NA,                            # 前两列不需要调色板
              rep(c("AUC_palette",NA), length.out = 8),           # AUC 用蓝色调色板
              rep(c("ACC_palette",NA), length.out = 8)
              ),
  options = I(c(
    list(list(width = 10, size = 2.5), list(width = 10, size = 2.5)),  
    rep(c(list(list(width = 2)), list(list(width = 2, overlay = TRUE, hjust = 0.15))), 8)
  ))
  
)


column_groups <- tribble(
  ~group,             ~palette,      
  "Dataset",             "overall",
  "Drug",             "overall",
  "AUC",              "AUC_palette",
  "ACC",              "ACC_palette",
)

palettes <- tribble(
  ~palette,             ~colours,
  "overall",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "AUC_palette",          grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(6, "Reds")))(101),
  "ACC_palette",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(101)
)


#################################### 手动画图例 ################################################
n_colors <- 5  # 需要的颜色数量
palette_names <- c("Blues", "Reds")  # 色卡名称
legend_labels <- c("AUC rank", "ACC rank")

# customized_legend <- function(n_colors = 5, palette_names = c('Blues', 'Reds'), 
#                               col_labels = c("", ""), label_size = 1.5, arrow = TRUE, legend_title = "Ranking") {
#   if (length(palette_names) != length(col_labels)) {
#     return('请确保颜色组的数量与标签数量一致！')
#   } else {
#     # 确定组的数量
#     n_groups <- length(palette_names)
#     
#     # 每组之间的间隔
#     row_spacing <- 1.1  # 调整组之间的垂直间距
#     
#     # 初始化 ggplot 对象
#     p <- ggplot() + theme_void() + 
#       theme(
#         legend.position = "none",  # 移除默认图例
#         axis.text = element_blank(),
#         axis.title = element_blank()
#       )
#     
#     # 通过 for 循环逐行生成热图
#     for (i in seq_along(palette_names)) {
#       # 获取当前组的颜色
#       legend_colors <- RColorBrewer::brewer.pal(n = 9, name = palette_names[i])
#       legend_colors <- rev(legend_colors[seq(1, 9, length.out = n_colors)])  # 选择颜色并反转顺序
#       
#       # 构建当前组的坐标和数据
#       df <- data.frame(
#         x = 1:n_colors,
#         y = (n_groups - i) * row_spacing,  # 确保行间有间隔
#         color = legend_colors
#       )
#       
#       # 添加当前组的热图方块
#       p <- p + geom_tile(data = df, aes(x = x, y = y, fill = color), color = "black", size = 0.8) + scale_fill_identity() 
#       
#       # 添加当前组的标签
#       p <- p + annotate(
#         "text", x = 0.5, y = df$y[1], label = col_labels[i], hjust = 1, size = label_size, fontface = "bold"
#       )
#     }
#     
#     # 如果需要箭头和标题
#     if (arrow) {
#       max_y <- n_groups - 1   # 箭头的 y 坐标起点
#       
#       # 添加箭头
#       p <- p + annotate(
#         "segment", x = 0.5, xend = n_colors + 0.5, y = max_y + 0.8, yend = max_y + 0.8,
#         arrow = arrow(length = unit(0.2, "inches"), type = 'closed'), size = 1.3, color = "black"
#       )
#       
#       # 添加箭头两端的标签和标题
#       p <- p + 
#         annotate("text", x = 0.3, y = max_y + 0.8, label = "4", size = 5, fontface = "bold") +
#         annotate("text", x = n_colors + 0.7, y = max_y + 0.8, label = "1", size = 5, fontface = "bold") +
#         annotate("text", x = n_colors / 2 + 0.5, y = max_y + 1, label = legend_title, size = 5, fontface = "bold")
#     }
#     
#     return(p)
#   }
# }

#另一种方法，就是需要手动在ppt中拼接
customized_legend <- function(n_colors = 5, palette_names = c('Blues', 'Reds'),
                              col_labels = c("", ""), label_size = 1.5, arrow = TRUE, legend_title = "") {
  if (length(palette_names) != length(col_labels)) {
    return('Please make sure the number of colors is consistant with that of labels!')
    break
  } else {
    # 确定总组数
    n_groups <- length(palette_names)

    # 设置绘图区域
    plot(NULL, xlim = c(0, n_colors + 1), ylim = c(0, n_groups * 2), type = "n", xlab = "", ylab = "", axes = FALSE)
    col_labels = col_labels
    # 遍历每组颜色并绘制
    for (i in seq_along(palette_names)) {
      # 获取当前组的颜色
      legend_colors <- RColorBrewer::brewer.pal(n = 9, name = palette_names[i])
      legend_colors <- rev(legend_colors[seq(1, 9, length.out = n_colors)])  # 按比例选取并反转颜色

      # 绘制矩形条
      rect(
        xleft = 1:n_colors - 0.5,
        ybottom = n_groups - i + 0.7,
        xright = 1:n_colors + 0.5,
        ytop = n_groups - i + 1.45,
        col = legend_colors,
        border = "black",
        lwd = 2
      )

      # 添加组名称标签
      text(x = 0.4, y = n_groups - i + 1, labels = col_labels[i], adj = 1, cex = label_size, font = 2)

      if (arrow) {


        arrows(
          x0 = 0.5,
          y0 = n_groups + 0.8,   # 箭头起点 y 坐标
          x1 = n_colors + 0.5,
          y1 = n_groups + 0.8,   # 箭头终点 y 坐标
          length = 0.2,
          angle = 20,
          col = "black",
          lwd = 3.5
        )

        text(
          x = 0.3,
          y = n_groups + 0.8,     # 调整 y 坐标到起点文字位置
          labels = "4",         # 起点文字
          adj = 0.5,            # 水平居中
          cex = 1.5,            # 字体大小
          font = 2              # 字体样式加粗
        )

        text(
          x = n_colors + 0.7,
          y = n_groups + 0.8,     # 调整 y 坐标到起点文字位置
          labels = "1",         # 起点文字
          adj = 0.5,            # 水平居中
          cex = 1.5,            # 字体大小
          font = 2              # 字体样式加粗
        )

        text(
          x = n_colors / 2 + 0.5,
          y = n_groups + 1.2,     # 调整 y 坐标到起点文字位置
          labels = legend_title,         # 起点文字
          adj = 0.5,            # 水平居中
          cex = 1.5,            # 字体大小
          font = 2              # 字体样式加粗
        )
      }
    }
  }
}

png('legend.png', units = 'in', width = 6, height = 3.5, res = 300)
customized_legend(n_colors = n_colors, palette_names = palette_names, label_size = 1.5, legend_title = "Ranking")
dev.off()
################################################################################################


png('funkyheatmap2.png', units = 'in', width = 6, height = 4, res = 300)
funky_heatmap(data = data_long, column_info = column_info, palettes = palettes,  
              add_abc = F, column_groups = column_groups, 
              position_args = position_arguments(
                col_annot_offset = 0.5, col_space = 0.5, expand_xmin = 2, expand_ymin = 6))

dev.off()

