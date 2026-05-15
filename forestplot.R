

library(readr)
DRD3data <- read_csv("DRD3data.csv")



library(readr)
dataDRD4 <- read_csv("dataDRD4.csv")



library(readr)
dataHTR1A <- read_csv("dataHTR1A.csv")


library(readr)
dataHTR2A <- read_csv("dataHTR2A.csv")




library(readr)
dataHTR2B <- read_csv("dataHTR2B.csv")




# 合并所有数据框
library(dplyr)

# 假设你的数据框已经被加载到环境中：DRD3data, dataDRD4, dataHTR1A, dataHTR2A, dataHTR2B

merged_data <- bind_rows(DRD3data, dataDRD4, dataHTR1A, dataHTR2A, dataHTR2B)

# 查看合并后的数据
head(merged_data)






data <- merged_data


# 加载 ggplot2 包
library(ggplot2)

# 假设已经有一个数据框 `data`，其中包含 OR 值及其置信区间（ci_lower 和 ci_upper）
# 确保数据框包含这些列：'odds_ratios', 'ci_lower', 'ci_upper', 'gene_name' 或其他标识符

# 创建一个包含 OR 和置信区间的数据框，假设有一列 'gene_name' 来标识每个基因/变异
forest_data <- data.frame(
  gene_name = data$`Protein Consequence`,  # 假设你有一个基因名列
  odds_ratio = data$odds_ratios,
  ci_lower = data$ci_lower,
  ci_upper = data$ci_upper,
  vep=data$`VEP Annotation`,
  p_value=data$p_values,
  GENE=merged_data$GENE
)


library(ggplot2)
library(dplyr)
library(stringr)

# 数据预处理：按VEP分组，组内按OR降序排列
forest_data <- forest_data %>%
  mutate(vep = factor(vep, 
                      levels = c("stop_gained", "frameshift_variant","splice_acceptor_variant", "missense_variant"),
                      ordered = TRUE)) %>%
  arrange(vep, desc(odds_ratio)) %>%
  mutate(
    rank = row_number(),
    gene_name = str_trunc(gene_name, 25, "right"),
    # 根据p_value创建星号标记（显著结果加*）
    star = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""  # 不显著则为空字符串
    ),
    # 更新显著性标识（基于p_value）
    is_significant = ifelse(p_value < 0.05, "Significant", "Not significant")
  )

# 创建自定义颜色映射
vep_colors <- c(
  "stop_gained" = "#D55E00",
  "frameshift_variant" = "#8C1B02",  # 使用更深的红色
  "missense_variant" = "#0072B2",
  "splice_acceptor_variant"="#9900CC"
)

# 创建森林图
ggplot(forest_data, aes(x = odds_ratio, y = rank)) +
  # 添加参考线
  geom_vline(
    xintercept = 1, 
    color = "#7f7f7f", 
    linetype = "dashed",
    linewidth = 0.8
  ) +
  # 置信区间水平线 - 按VEP注释分组着色
  geom_errorbarh(
    aes(xmin = ci_lower, xmax = ci_upper, color = vep),
    height = 0.15,
    linewidth = 0.8,
    alpha = 0.9
  ) +
  # 点估计值 - 显著结果高亮
  geom_point(
    aes(color = vep, fill = is_significant),
    size = 3.5,
    shape = 21,
    stroke = 1.0
  ) +
  # 添加显著性星号标记
  geom_text(
    aes(label = star),
    size = 6,  # 星号大小
    color = "black",
    nudge_y = 0.25,  # 垂直位置调整（在点上方）
    fontface = "bold"
  ) +
  # 反转Y轴并设置标签
  scale_y_reverse(
    breaks = 1:nrow(forest_data),
    labels = forest_data$gene_name,
    expand = expansion(add = 0.8)
  ) +
  # 对数坐标处理大范围OR值
  scale_x_log10(
    breaks = c(0.1, 0.5, 1, 2, 5, 10, 50, 100),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  # VEP注释颜色标度
  scale_color_manual(
    name = "Mutation type",
    values = vep_colors,
    guide = guide_legend(order = 1)
  ) +
  # 显著结果填充色标度
  scale_fill_manual(
    name = "Significance",
    values = c("Significant" = "#e74c3c", "Not significant" = "white"),
    guide = guide_legend(order = 2)
  ) +
  # 添加标签和标题
  labs(
    title = "Gene Variants Odds Ratios by VEP Category",
    x = "Odds Ratio (log scale)",
    y = "",
    caption = "***: p<0.001; **: p<0.01; *: p<0.05"
  ) +
  # 应用主题
  theme_minimal(base_size = 12) +
  theme(
    plot.margin = margin(1, 1.2, 0.8, 2, "cm"),
    panel.grid.major.x = element_line(color = "#f0f0f0"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(color = "#5d5d5d"),
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(
      hjust = 1, 
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(color = "#5d5d5d"),
    plot.caption = element_text(color = "#7f7f7f", size = 9, lineheight = 1.2),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  # 添加组间分隔线
  geom_hline(
    yintercept = c(
      max(which(forest_data$vep == "stop_gained")) + 0.5,
      max(which(forest_data$vep == "frameshift_variant")) + 0.5,
      max(which(forest_data$vep == "splice_acceptor_variant")) + 0.5
    ),
    color = "#b0b0b0",
    linetype = "dotted",
    alpha = 0.7
  ) +
  # 添加组标签注释
  annotate("text",
           x = 0.05, 
           y = mean(which(forest_data$vep == "stop_gained")),
           label = "Stop Gained",
           hjust = 0,
           size = 3.8,
           color = "#D55E00",
           fontface = "bold") +
  annotate("text",
           x = 0.05,
           y = mean(which(forest_data$vep == "frameshift_variant")),
           label = "Frameshift",
           hjust = 0,
           size = 3.8,
           color = "#8C1B02",
           fontface = "bold") +
  annotate("text",
           x = 0.05,
           y = mean(which(forest_data$vep == "splice_acceptor_variant")),
           label = "Splicing",
           hjust = 0,
           size = 3.8,
           color = "#9900CC",
           fontface = "bold") +
  annotate("text",
           x = 0.05,
           y = mean(which(forest_data$vep == "missense_variant")),
           label = "Missense",
           hjust = 0,
           size = 3.8,
           color = "#0072B2",
           fontface = "bold")










