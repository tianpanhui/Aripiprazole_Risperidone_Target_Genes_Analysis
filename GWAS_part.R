

library(qqman)
library(ggplot2)
library(dplyr)

diagnosis_age <-read.table("DRD3.biallele.qc.LD.filtered.diagnosis_age.diagnosis_age.glm.linear")
fsiq <-read.table("DRD3.biallele.qc.LD.filtered.fsiq.fsiq.glm.linear")
nviq <-read.table("DRD3.biallele.qc.LD.filtered.nviq.nviq.glm.linear")
viq <-read.table("DRD3.biallele.qc.LD.filtered.viq.viq.glm.linear")
scq <-read.table("DRD3.biallele.qc.LD.filtered.scq.scq_total_final_score.glm.linear")
rbsr <-read.table("DRD3.biallele.qc.LD.filtered.rbsr.rbsr_total_final_score.glm.linear")
vineland <-read.table("DRD3.biallele.qc.LD.filtered.vineland_abc_ss_latest.vineland_abc_ss_latest.glm.linear")
used_words_age_mos <-read.table("DRD3.biallele.qc.LD.filtered.used_words_age_mos.used_words_age_mos.glm.linear")
walked_age_mos <-read.table("DRD3.biallele.qc.LD.filtered.walked_age_mos.walked_age_mos.glm.linear")


diagnosis_age$source <- "diagnosis_age"
fsiq$source <- "fsiq"
nviq$source <- "nviq"
viq$source <- "viq"
scq$source <- "scq"
rbsr$source <- "rbsr"
vineland$source <- "vineland"
used_words_age_mos$source <- "used_words_age_mos"
walked_age_mos$source <- "walked_age_mos"


# Rename columns to match the format
names(diagnosis_age) <- c("CHR", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE", "source")
names(fsiq) <- names(diagnosis_age)
names(nviq) <- names(diagnosis_age)
names(viq) <- names(diagnosis_age)
names(scq) <- names(diagnosis_age)
names(rbsr) <- names(diagnosis_age)
names(vineland) <- names(diagnosis_age)
names(used_words_age_mos) <- names(diagnosis_age)
names(walked_age_mos) <- names(diagnosis_age)


# Apply Bonferroni and FDR corrections for each data frame
# Diagnosis Age
p_values <- diagnosis_age$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
diagnosis_age$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
diagnosis_age$p_adjusted_fdr <- p_adjusted_fdr

# FSIQ
p_values <- fsiq$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
fsiq$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
fsiq$p_adjusted_fdr <- p_adjusted_fdr

# NVIQ
p_values <- nviq$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
nviq$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
nviq$p_adjusted_fdr <- p_adjusted_fdr

# VIQ
p_values <- viq$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
viq$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
viq$p_adjusted_fdr <- p_adjusted_fdr

# SCQ
p_values <- scq$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
scq$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
scq$p_adjusted_fdr <- p_adjusted_fdr

# RBSR
p_values <- rbsr$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
rbsr$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
rbsr$p_adjusted_fdr <- p_adjusted_fdr

# Vineland
p_values <- vineland$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
vineland$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
vineland$p_adjusted_fdr <- p_adjusted_fdr

# Used Words Age Mos
p_values <- used_words_age_mos$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
used_words_age_mos$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
used_words_age_mos$p_adjusted_fdr <- p_adjusted_fdr

# Walked Age Mos
p_values <- walked_age_mos$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
walked_age_mos$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
walked_age_mos$p_adjusted_fdr <- p_adjusted_fdr






data <- rbind(diagnosis_age, fsiq,nviq,viq,scq,rbsr,vineland,used_words_age_mos,walked_age_mos)

data <- rbind( fsiq,nviq,viq,scq,rbsr)




















manhattan(scq,#数据
          col = c('#30A9DE','#EFDC05','#E53A40','#090707'),#交替使用颜色展示
          suggestiveline = -log10(1e-05),#－log10(1e－5)处添加"suggestive"横线
          genomewideline = -log10(5e-08),#－log10(5e－10)处添加"genome-wide sigificant"横线
          highlight = snpsOfInterest,#内置高亮的snp数据， 也可以对snpOfInterest进行设置
          annotatePval = 0.05,#标记p值小于0.05的点
          annotateTop = T,#如果为T，则仅批注低于注解阈值的每个染色体上的顶部点，为F则标记所有小于注解阈值的点。
          chr = "CHR",          # 染色体列
          bp = "POS",            # 物理位置列
          p = "P",              # p 值列
          snp = "ID",          # SNP ID 列
          main = "DRD3"#标题
)



qq(scq$P)





# Create the bubble plot
ggplot(scq, aes(x = POS, y = -log10(p_adjusted_fdr), size = abs(BETA), color = BETA)) +
  geom_point(alpha = 1, stroke = 1.5,shape=21) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Bubble Plot of POS vs Adjusted FDR Value with BETA as Size & Color",
    x = "POS Coordinate",
    y = "Adjusted FDR Value",
    color = "BETA Value"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
















# Assign different shapes to each source
shape_values <- c(
  "diagnosis_age" = 21,
  "fsiq" = 22,
  "nviq" = 23,
  "viq" = 24,
  "scq" = 25,
  "rbsr" = 7,
  "vineland" = 8,
  "used_words_age_mos" = 5,
  "walked_age_mos" = 6
)

# Create the volcano plot with different shapes and fill colors
# Set pvalue and logFC thresholds
cut_off_pvalue <- 0.01
cut_off_BETA <- 1

# Plot volcano plot with shape size representing BETA value
ggplot(data, aes(x = BETA, y = -log10(p_adjusted_fdr), shape = source, fill = source, size = abs(BETA))) +
  geom_point(alpha = 0.7) +
  scale_shape_manual(values = shape_values) +
  scale_fill_manual(values = c(
    "diagnosis_age" = "#ff9999",
    "fsiq" = "#99ccff",
    "nviq" = "#ffcc99",
    "viq" = "#ccff99",
    "scq" = "#ffccff",
    "rbsr" = "#99ffcc",
    "vineland" = "#ff9933",
    "used_words_age_mos" = "#66b3ff",
    "walked_age_mos" = "#c2c2f0"
  )) +
  geom_vline(xintercept = c(-cut_off_BETA, cut_off_BETA), lty = 4, col = "black", lwd = 1) +
  geom_hline(yintercept = -log10(cut_off_pvalue), lty = 4, col = "black", lwd = 1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand = expansion(mult = c(0, 1.5))) +
  labs(x = "BETA value", y = "-log10 (FDR)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right", 
        legend.title = element_blank(),
        text = element_text(size = 14))







# Create QQ plot with different shapes for each source
# Calculate theoretical quantiles
qq_data <- data %>% 
  arrange(P) %>% 
  mutate(theoretical = -log10(ppoints(n())), observed = -log10(P))

# Plot QQ plot
ggplot(qq_data, aes(x = theoretical, y = observed, shape = source, fill = source)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_shape_manual(values = shape_values) +
  scale_fill_manual(values = c(
    "diagnosis_age" = "#ff9999",
    "fsiq" = "#99ccff",
    "nviq" = "#ffcc99",
    "viq" = "#ccff99",
    "scq" = "#ffccff",
    "rbsr" = "#99ffcc",
    "vineland" = "#ff9933",
    "used_words_age_mos" = "#66b3ff",
    "walked_age_mos" = "#c2c2f0"
  )) +
  geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 1) +
  labs(x = "Theoretical Quantiles (-log10)", y = "Observed Quantiles (-log10)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right", 
        legend.title = element_blank(),
        text = element_text(size = 14))








############## 二分类变量 可视化  #############




cognitive_impairment_latest <-read.table("DRD3.biallele.qc.LD.filtered.cognitive_impairment_latest.cognitive_impairment_latest.glm.logistic.hybrid")
current_depend_adult <-read.table("DRD3.biallele.qc.LD.filtered.current_depend_adult.current_depend_adult.glm.logistic.hybrid")
dcdq_dcd <-read.table("DRD3.biallele.qc.LD.filtered.dcdq_dcd.dcdq_dcd.glm.logistic.hybrid")
regress_lang_y_n <-read.table("DRD3.biallele.qc.LD.filtered.regress_lang_y_n.regress_lang_y_n.glm.logistic.hybrid")
regress_other_y_n <-read.table("DRD3.biallele.qc.LD.filtered.regress_other_y_n.regress_other_y_n.glm.logistic.hybrid")



cognitive_impairment_latest$source <- "cognitive_impairment_latest"
current_depend_adult$source <- "current_depend_adult"
dcdq_dcd$source <- "dcdq_dcd"
regress_lang_y_n$source <- "regress_lang_y_n"
regress_other_y_n$source <- "regress_other_y_n"

cognitive_impairment_latest$colors <- "#1f77b4"
current_depend_adult$colors <- "#ff7f0e"
dcdq_dcd$colors <- "#2ca02c"
regress_lang_y_n$colors <- "#d62728"
regress_other_y_n$colors <- "#9467bd"














# Rename columns to match the format
names(cognitive_impairment_latest) <- c("CHR", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ", "FIRTH?","TEST", "OBS_CT", "OR", "LOG(OR)_SE", "Z_STAT", "P", "ERRCODE", "source","colors")
names(current_depend_adult) <- names(cognitive_impairment_latest)
names(dcdq_dcd) <- names(cognitive_impairment_latest)
names(regress_lang_y_n) <- names(cognitive_impairment_latest)
names(regress_other_y_n) <- names(cognitive_impairment_latest)




# Apply Bonferroni and FDR corrections for each data frame
# Diagnosis Age
p_values <- cognitive_impairment_latest$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
cognitive_impairment_latest$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
cognitive_impairment_latest$p_adjusted_fdr <- p_adjusted_fdr

# FSIQ
p_values <- current_depend_adult$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
current_depend_adult$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
current_depend_adult$p_adjusted_fdr <- p_adjusted_fdr

# NVIQ
p_values <- dcdq_dcd$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
dcdq_dcd$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
dcdq_dcd$p_adjusted_fdr <- p_adjusted_fdr

# VIQ
p_values <- regress_lang_y_n$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
regress_lang_y_n$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
regress_lang_y_n$p_adjusted_fdr <- p_adjusted_fdr

# SCQ
p_values <- regress_other_y_n$P
p_adjusted <- p.adjust(p_values, method = "bonferroni")
regress_other_y_n$p_adjusted <- p_adjusted
p_adjusted_fdr <- p.adjust(p_values, method = "fdr")
regress_other_y_n$p_adjusted_fdr <- p_adjusted_fdr




data_cat <- rbind(cognitive_impairment_latest, current_depend_adult,dcdq_dcd,regress_lang_y_n,regress_other_y_n)








# 安装并加载所需包
library(forestplot)

df <- data_cat
df <- data_cat[data_cat$P < 0.05,]
# 自定义颜色映射
source_colors <- c(
  "cognitive_impairment_latest" = "#1f77b4", # 蓝色
  "current_depend_adult" = "#ff7f0e",        # 橙色
  "dcdq_dcd" = "#2ca02c",                     # 绿色
  "regress_lang_y_n" = "#d62728",            # 红色
  "regress_other_y_n" = "#9467bd"           # 紫色
)

# 计算置信区间
df$lower_CI <- exp(log(df$OR) - 1.96 * df$`LOG(OR)_SE`)
df$upper_CI <- exp(log(df$OR) + 1.96 * df$`LOG(OR)_SE`)



# 创建森林图数据
tabletext <- cbind(
  c("Phenotype", df$source),
  c("OR", sprintf("%.2f", df$OR)),
  c("95% CI", paste(sprintf("%.2f", df$lower_CI), "-", sprintf("%.2f", df$upper_CI)))
)

# 绘制森林图
forestplot(
  labeltext = tabletext,
  mean = c(NA, df$OR),  # OR值（NA表示标题行）
  lower = c(NA, df$lower_CI),  # 下置信区间
  upper = c(NA, df$upper_CI),  # 上置信区间
  zero = 1,  # 指定零线（通常为1，用于OR值）
  col = forestplot::fpColors(lines = source_colors[df$source], box = source_colors[df$source]),
  xlab = "Odds Ratio (OR)",
  title = "DRD3"
)

library(forestplot)
library(grid)  # 需要引入 grid 包来使用 gpar()

forestplot(
  labeltext = tabletext,  # 标签文本
  mean = c(NA, df$OR),  # OR值（NA表示标题行）
  lower = c(NA, df$lower_CI),  # 下置信区间
  upper = c(NA, df$upper_CI),  # 上置信区间
  zero = 1,  # 指定零线（通常为1，用于OR值）
  col = forestplot::fpColors(lines = source_colors[df$source], box = source_colors[df$source]),
  xlab = "Odds Ratio (OR)",
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.6),  # 标签文字大小
    title = gpar(cex = 2),    # 标题文字大小
    xlab = gpar(cex = 1.2),   # x轴标签文字大小
    ticks = gpar(cex = 1.1),  # 刻度文字大小
    legend = gpar(cex = 1.3)  # 图例文字大小（如果有图例）
  )
)















library(dplyr)



# Assign different shapes to each source
shape_values <- c(
  "cognitive_impairment_latest" = 1,
  "current_depend_adult" = 2,
  "dcdq_dcd" = 3,
  "regress_lang_y_n" = 4,
  "regress_other_y_n" = 20
)


# Create QQ plot with different shapes for each source
# Calculate theoretical quantiles
qq_data <- data_cat %>% 
  arrange(P) %>% 
  mutate(theoretical = -log10(ppoints(n())), observed = -log10(P))

# Plot QQ plot
ggplot(qq_data, aes(x = theoretical, y = observed, shape = source, fill = source)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_shape_manual(values = shape_values) +
  scale_fill_manual(values = c(
    "cognitive_impairment_latest" = "#1f77b4",  # 蓝色
    "current_depend_adult" = "#ff7f0e",           # 橙色
    "dcdq_dcd" = "#2ca02c",           # 绿色
    "regress_lang_y_n" = "#d62728",            # 红色
    "regress_other_y_n" = "#9467bd"             # 紫色
  )) +
  geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 1) +
  labs(x = "Theoretical Quantiles (-log10)", y = "Observed Quantiles (-log10)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right", 
        legend.title = element_blank(),
        text = element_text(size = 14))

















# 使用 Ensembl Variation 数据集
mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "https://useast.ensembl.org")
options(timeout = 600)  # 设置超时时间为 600 秒 (10 分钟)


# 列出可用的属性，以确认是否有相关的 SNP 属性
attributes <- listAttributes(mart)
head(attributes)
print(attributes)
listFilters(mart)
# 创建包含染色体位置的列表
chromosome_list <- c("1", "2", "3")  # 替换为您的染色体编号
position_list <- c(123456, 234567, 345678)  # 替换为您的位置
# 获取对应 SNP 的注释信息
snp_annotation <- getBM(
  attributes = c("refsnp_id"),
  filters = c("chr_name", "start"),
  values = list(chromosome_list, position_list),
  mart = mart
)
# 查看注释结果
print(snp_annotation)







# 列出可用的属性，以确认是否有相关的 SNP 属性
attributes <- listAttributes(mart)
head(attributes)
print(attributes)
listFilters(mart)
snp_list <- c("rs123", "rs456", "rs789")  # 用您的 SNP rsID 替换这里
snp_annotation <- getBM(
     attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id", "associated_gene","ensembl_gene_name","variation_names","study_name"),
     filters = "snp_filter",
      values = snp_list,
      mart = mart
 )
print(snp_annotation)















