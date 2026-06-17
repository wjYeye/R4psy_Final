# =============================================================================
# Hurst & Kavanagh (2017). Life history strategies and psychopathology: the faster the life strategies, the more symptoms of psychopathology.Evolution and Human Behavior, 38, 1–8.
# 数据来源：作者公开补充材料（mmc1.csv），N = 138
# =============================================================================


# -- 0  依赖包 ---------------------------------------------------------------
library(tidyverse)
library(psych)
library(ppcor)
library(ggcorrplot)
library(lm.beta)
library(flextable)
library(officer)

theme_set(theme_classic(base_size = 12))


# -- 1  数据读取与变量构建（§ 2 Method）-------------------------------------
df_raw <- read.csv("D:\\桌面\\R_final\\0.raw data\\1-s2.0-S1090513816301118-mmc1.csv",
                   stringsAsFactors = FALSE)

df <- df_raw %>%
  dplyr::select(
    Age, Sex, Bio_Sib, Half_Sib, Step_Sib, Stepparent, Stepparent_N,  # 同时保留二值 Stepparent 与计数 Stepparent_N
    Attach_Total, AAQ_Avail, AAQ_Angry, AAQ_Goal,
    Aggresion_Total, AQ_Physical, AQ_Verbal, AQ_Anger, AQ_Hostility,
    DSM5_Total, MiniK_Total, HKSS_Total, SH_Total,
    SH_01:SH_17, starts_with("DSM5_")
  ) %>%
  dplyr::rename(
    age = Age, sex = Sex,
    bio_sib = Bio_Sib, half_sib = Half_Sib,
    step_sib = Step_Sib, step_par = Stepparent_N,  # 偏相关采用 Stepparent_N（继父母人数 0/1/2），与原文 Table 5 一致
    attach = Attach_Total, avail = AAQ_Avail,
    angry = AAQ_Angry, goal = AAQ_Goal,
    aggress = Aggresion_Total, phys = AQ_Physical,
    verbal = AQ_Verbal, anger = AQ_Anger, hostil = AQ_Hostility,
    psych = DSM5_Total, minik = MiniK_Total,
    hkss = HKSS_Total, selfharm = SH_Total
  ) %>%
  mutate(
    sex_f  = factor(sex, levels = c(1,2), labels = c("男性","女性")),
    sex_01 = ifelse(sex == 2, 1, 0)
  )

sh_cols <- paste0("SH_", sprintf("%02d", 1:17))
df[sh_cols] <- lapply(df[sh_cols], function(x) suppressWarnings(as.numeric(x)))


# -- 2  辅助函数 -------------------------------------------------------------

cor_fmt <- function(x, y) {
  res <- cor.test(x, y)
  r   <- round(res$estimate, 2)
  stars <- ifelse(res$p.value < .001, "***",
           ifelse(res$p.value < .01,  "**",
           ifelse(res$p.value < .05,  "*", "")))
  paste0(r, stars)
}

pcor_age <- function(x_var, y_var) {
  res <- ppcor::pcor.test(df[[x_var]], df[[y_var]], df[["age"]])
  list(pr = round(res$estimate, 2), p = round(res$p.value, 3),
       sig = res$p.value < .002)
}

pcor_block <- function(vars, labels) {
  map2_dfr(vars, labels, function(v, lbl) {
    pm <- pcor_age("minik", v)
    ph <- pcor_age("hkss",  v)
    tibble(变量 = lbl,
           pr_minik = pm$pr, p_minik = pm$p, sig_minik = pm$sig,
           pr_hkss  = ph$pr, p_hkss  = ph$p, sig_hkss  = ph$sig)
  })
}

run_moderation <- function(outcome, predictor) {
  fml <- as.formula(sprintf("%s ~ age + %s * sex_01", outcome, predictor))
  s   <- summary(lm(fml, data = df))$coefficients
  int <- paste0(predictor, ":sex_01")
  p   <- if (int %in% rownames(s)) round(s[int, "Pr(>|t|)"], 3) else NA
  label <- if (!is.na(p) && p < .05) "【显著】" else "【不显著，与论文一致】"
  cat(sprintf("  %-8s ~ %-8s × sex：交互项 p = %.3f  %s\n",
              outcome, predictor, p, label))
}

format_reg <- function(mod) {
  b    <- coef(mod)[-1]
  ci   <- confint(mod)[-1, ]
  t_v  <- summary(mod)$coefficients[-1, "t value"]
  beta <- lm.beta(mod)$standardized.coefficients[-1]
  p_v  <- summary(mod)$coefficients[-1, "Pr(>|t|)"]
  data.frame(
    预测变量 = c("年龄（协变量）","半亲兄弟姐妹（协变量）",
                 "继兄弟姐妹（协变量）","同父母兄弟姐妹（协变量）","依恋总分"),
    b = round(b,2), beta = round(beta,2), t_val = round(t_v,2),
    ci_low = round(ci[,1],2), ci_high = round(ci[,2],2),
    p_val = round(p_v,3), stringsAsFactors = FALSE)
}

add_sig <- function(p) {
  ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
}

ft_style <- function(ft, footer = NULL, caption = NULL) {
  ft <- ft %>%
    bold(part = "header") %>%
    align(j = -1, align = "center", part = "all") %>%
    align(j = 1,  align = "left",   part = "body") %>%
    border_outer(part = "all", border = fp_border(width = 1.5)) %>%
    autofit()
  if (!is.null(footer))  ft <- add_footer_lines(ft, footer)
  if (!is.null(caption)) ft <- set_caption(ft, caption)
  ft
}


# =============================================================================
# § 2.2  量表描述统计与信度 -> Table 1   
# =============================================================================

vars_t1 <- c("attach","avail","goal","angry",
             "aggress","phys","verbal","anger","hostil","minik","hkss")

desc_t1 <- psych::describe(df[, vars_t1]) %>%
  dplyr::select(mean, sd) %>% round(2)

# 信度系数说明：psych::alpha 同时给出 raw_alpha 与 std.alpha 两种 Cronbach's α。
#   经真实数据核验，std.alpha（标准化 α）在全部 11 个量表上与原文 Table 1 完全一致
#   （四舍五入后 11/11，δ=0%），表明原文报告的是标准化 α，故此处取 std.alpha。
#   check.keys=TRUE 用于自动核查反向题；本数据已规范计分、不触发反转，保留无害。
calc_alpha <- function(items) {
  suppressMessages(
    psych::alpha(df_raw[, items], check.keys = TRUE, warnings = FALSE)$total$std.alpha
  )
}

alpha_calc <- c(
  calc_alpha(paste0("AAQ_0", 1:9)),
  calc_alpha(paste0("AAQ_0", 4:6)),
  calc_alpha(paste0("AAQ_0", 7:9)),
  calc_alpha(paste0("AAQ_0", 1:3)),
  calc_alpha(paste0("BPAS_", sprintf("%02d", 1:29))),
  calc_alpha(paste0("BPAS_", sprintf("%02d", 1:9))),
  calc_alpha(paste0("BPAS_", sprintf("%02d", 10:14))),
  calc_alpha(paste0("BPAS_", sprintf("%02d", 15:21))),
  calc_alpha(paste0("BPAS_", sprintf("%02d", 22:29))),
  calc_alpha(paste0("MiniK_", sprintf("%02d", 1:20))),
  calc_alpha(paste0("HKSS_", sprintf("%02d", 1:22)))
)
alpha_calc <- round(alpha_calc, 2)

table1_df <- data.frame(
  量表 = c("依恋（AAQ）","  可及性","  目标校正伙伴关系","  愤怒苦恼",
           "攻击性（Buss-Perry）","  身体攻击","  言语攻击","  愤怒","  敌意",
           "生活史策略 Mini-K","生活史策略 HKSS"),
  M = desc_t1$mean, SD = desc_t1$sd, Alpha = alpha_calc,
  stringsAsFactors = FALSE)

ft_table1 <- flextable(table1_df) %>%
  set_header_labels(量表="量表", M="M", SD="SD", Alpha="α") %>%
  hline(i = c(4,9), border = fp_border(width=0.5, style="dashed")) %>%
  ft_style(footer = "注：M=均值；SD=标准差；α=Cronbach's α；分量表以缩进表示",
           caption = "Table 1. 各量表均值、标准差与 Cronbach's α")


# =============================================================================
# § 2.2.4  DSM-5 各条目达阈值百分比 -> Table 2
# =============================================================================
# 论文阈值：多数条目 ≥2；条目 11-13、21-23 为 ≥1（量表 0–4 计分）。
#
# 计分对齐说明：原始数据 DSM5_01–23 按 1–5 计分（1=无…5=严重），而论文阈值
#   定义在 0–4 量表上（多数条目 ≥2，条目 11-13、21-23 为 ≥1）。故先将数据
#   vals - 1 对齐到 0–4 量表，再套用论文阈值，所得百分比方与原文 Table 2 一致。
# -----------------------------------------------------------------------------

special_idx <- c(11,12,13,21,22,23)
threshold_pct <- sapply(1:23, function(i) {
  vals  <- suppressWarnings(as.numeric(df_raw[[paste0("DSM5_", sprintf("%02d",i))]]))
  vals0 <- vals - 1                       # 1–5 数据 → 对齐到论文 0–4 量表
  thr   <- ifelse(i %in% special_idx, 1, 2)   # 论文阈值：特殊≥1，其余≥2
  round(mean(vals0 >= thr, na.rm = TRUE) * 100, 1)
})

table2_df <- data.frame(
  条目 = paste0("DSM5_", sprintf("%02d",1:23)),
  症状描述 = c(
    "做事提不起兴趣或无法享受乐趣（抑郁）",
    "感到情绪低落、抑郁或绝望（抑郁）",
    "比平时更容易烦躁、暴躁或愤怒（愤怒）",
    "比平时睡眠少，但精力充沛（躁狂）",
    "比平时开始更多计划或做更多冒险的事（躁狂）",
    "感到紧张、焦虑、害怕、担忧或不安（焦虑）",
    "感到恐慌或极度恐惧（焦虑）",
    "回避令你感到焦虑的情境（焦虑）",
    "无法解释的疼痛（躯体症状）",
    "感觉疾病未受到足够重视（躯体症状）",
    "有伤害自己的想法（自杀意念）*",
    "听到别人听不到的声音（精神病）*",
    "感觉有人能听到你的想法（精神病）*",
    "睡眠问题影响整体睡眠质量（睡眠障碍）",
    "记忆或方向感问题（记忆）",
    "反复闯入性的不愉快想法（重复思维）",
    "感到被迫反复执行某些行为（重复行为）",
    "感到与自身或周围环境脱离（解离）",
    "不知道自己真正是谁或想要什么（人格功能）",
    "感觉与他人不亲近或无法享受人际关系（人格功能）",
    "两周内单次饮酒 ≥4 杯（物质使用）*",
    "吸烟或使用烟草制品（物质使用）*",
    "未经处方自行使用或过量使用药物（物质使用）*"),
  百分比 = threshold_pct,
  stringsAsFactors = FALSE)

ft_table2 <- flextable(table2_df) %>%
  set_header_labels(条目="条目", 症状描述="症状描述", 百分比="达阈值 (%)") %>%
  ft_style(footer = "* 该条目阈值为 ≥1（而非 ≥2），适用于自杀意念、精神病及物质使用条目",
           caption = "Table 2. DSM-5 各条目达阈值的被试百分比")


# =============================================================================
# § 2.2.5  自伤行为频率（按性别）-> Table 3
# =============================================================================
# 计数口径说明：SH_01–17 编码为 1=有、2=无、空=缺失；统计"发生人数"时只计
#   == 1（回答"有"的被试）。
#
# 列序对齐说明：codebook(mmc2.xlsx) 显示数据列与论文行序在第 5–7 项不一致——
#   数据 SH_05=咬伤、SH_06=刻图案、SH_07=用力抓伤；论文 Table 3 行序为
#   5=刻图案、6=用力抓伤、7=咬伤。故按行为名建立映射（见 sh_map），而非"行号=列号"。
# -----------------------------------------------------------------------------

women_df <- df[df$sex == 2, ]; n_women <- nrow(women_df)  # 92
men_df   <- df[df$sex == 1, ]; n_men   <- nrow(men_df)    # 46

# 论文 Table 3 行（行为名）-> 对应的数据列（依据 codebook）
sh_labels <- c(
  "割伤（Cutting）","用香烟烫伤","用打火机或火柴烫伤",
  "在皮肤上刻字","在皮肤上刻图案","用力抓伤皮肤","咬伤自己",
  "用砂纸摩擦皮肤","将酸液滴在皮肤上","用漂白剂或清洁剂接触皮肤",
  "将针、订书钉或大头针刺入皮肤","用玻璃摩擦皮肤","骨折（故意）",
  "撞头","打自己","阻止伤口愈合","其他自伤方式")
sh_map <- c("SH_01","SH_02","SH_03","SH_04",
            "SH_06",   # 在皮肤上刻图案  = 数据 SH_06
            "SH_07",   # 用力抓伤皮肤    = 数据 SH_07
            "SH_05",   # 咬伤自己        = 数据 SH_05
            "SH_08","SH_09","SH_10","SH_11","SH_12","SH_13",
            "SH_14","SH_15","SH_16","SH_17")

count_sh_yes <- function(sub) {
  sapply(sh_map, function(col) {
    sum(suppressWarnings(as.numeric(sub[[col]])) == 1, na.rm = TRUE)  # 只计回答"有"(==1)
  })
}
wn <- count_sh_yes(women_df); mn <- count_sh_yes(men_df)

table3_df <- data.frame(
  编号 = 1:17,
  自伤行为 = sh_labels,
  女性 = paste0(wn, " (", round(wn/n_women*100,1), ")"),
  男性 = paste0(mn, " (", round(mn/n_men  *100,1), ")"),
  stringsAsFactors = FALSE)

ft_table3 <- flextable(table3_df) %>%
  set_header_labels(
    编号="编号", 自伤行为="自伤行为",
    女性 = paste0("女性 (n=",n_women,")\nn (%)"),
    男性 = paste0("男性 (n=",n_men,  ")\nn (%)")
  ) %>%
  ft_style(caption = "Table 3. 各类自伤行为的发生频率（按性别分列）")


# =============================================================================
# § 3.1  年龄与性别的相关 -> Table 4   
# =============================================================================

vars_t4 <- c("attach","avail","angry","goal",
             "aggress","phys","verbal","anger","hostil","psych","selfharm")
labels_t4 <- c("依恋总分","  可及性","  愤怒苦恼","  目标校正伙伴关系",
               "攻击性总分","  身体攻击","  言语攻击","  愤怒","  敌意",
               "心理病理总分（DSM-5）","自伤总频率")

table4_df <- map2_dfr(vars_t4, labels_t4, function(v, lbl) {
  tibble(变量=lbl, 年龄=cor_fmt(df[[v]],df$age), 性别=cor_fmt(df[[v]],df$sex))
})

ft_table4 <- flextable(table4_df) %>%
  set_header_labels(变量="变量", 年龄="年龄 r", 性别="性别 r") %>%
  ft_style(footer = "* p<.05  ** p<.01  *** p<.001",
           caption = "Table 4. 各变量与年龄、性别的零阶相关系数")


# =============================================================================
# § 3.2  偏相关（控制年龄，Bonferroni p < .002）与回归 -> Table 5 / Table 6
# =============================================================================
# 继父母变量说明："是否有继父母"一行采用计数变量 Stepparent_N（继父母人数 0/1/2），
#   而非二值 Stepparent（0=无 / 1=有）。同时有继父+继母时 Stepparent=1 但
#   Stepparent_N=2，两者并非简单 +1 关系。论文 Table 5 的 "Presence of stepparent"
#   数值与 Stepparent_N 对齐：
#     Mini-K  原文 -.19 / 二值 -.17 / Stepparent_N -.19
#     HKSS    原文 -.06 / 二值 -.14 / Stepparent_N -.06
#   回归 Table 6 不含此变量，不受影响。
# -----------------------------------------------------------------------------

fam_results <- pcor_block(
  c("bio_sib","half_sib","step_sib","step_par"),
  c("同父同母兄弟姐妹数","半亲兄弟姐妹数","继兄弟姐妹数","是否有继父母"))

psych_results <- pcor_block(
  c("psych","selfharm"),
  c("心理病理总分（DSM-5）","自伤总频率"))

cat("\n§ 3.2.2  性别调节效应（自伤结果变量）：\n")
run_moderation("selfharm","minik")
run_moderation("selfharm","hkss")

attach_results <- pcor_block(
  c("attach","avail","angry","goal"),
  c("依恋总分","  可及性","  愤怒苦恼","  目标校正伙伴关系"))

mod_minik <- lm(minik ~ age + half_sib + step_sib + bio_sib + attach, data=df)
mod_hkss  <- lm(hkss  ~ age + half_sib + step_sib + bio_sib + attach, data=df)

reg_m <- format_reg(mod_minik); reg_h <- format_reg(mod_hkss)
r2_m  <- round(summary(mod_minik)$r.squared, 2)
r2_h  <- round(summary(mod_hkss)$r.squared,  2)
f_m   <- round(summary(mod_minik)$fstatistic[1], 2)
f_h   <- round(summary(mod_hkss)$fstatistic[1],  2)

table6_df <- data.frame(
  预测变量 = reg_m$预测变量,
  b_m=reg_m$b, beta_m=paste0(reg_m$beta, add_sig(reg_m$p_val)), t_m=reg_m$t_val,
  ci_lo_m=reg_m$ci_low, ci_hi_m=reg_m$ci_high,
  b_h=reg_h$b, beta_h=paste0(reg_h$beta, add_sig(reg_h$p_val)), t_h=reg_h$t_val,
  ci_lo_h=reg_h$ci_low, ci_hi_h=reg_h$ci_high,
  stringsAsFactors = FALSE)

ft_table6 <- flextable(table6_df) %>%
  set_header_labels(
    预测变量="预测变量",
    b_m="b", beta_m="β", t_m="t", ci_lo_m="CI下", ci_hi_m="CI上",
    b_h="b", beta_h="β", t_h="t", ci_lo_h="CI下", ci_hi_h="CI上") %>%
  add_header_row(values=c("","Mini-K（因变量）","HKSS（因变量）"), colwidths=c(1,5,5)) %>%
  bold(i=5, part="body") %>%
  hline(i=4, border=fp_border(width=0.5, style="dashed")) %>%
  ft_style(
    footer  = paste0("Mini-K: R²=",r2_m,", F(5,132)=",f_m,", p<.001  |  ",
                     "HKSS: R²=",r2_h,", F(5,132)=",f_h,", p<.001\n",
                     "* p<.05  ** p<.01  *** p<.001"),
    caption = "Table 6. 依恋（感知父母支持）预测生活史策略的回归分析")

aggr_results <- pcor_block(
  c("aggress","phys","verbal","anger","hostil"),
  c("攻击性总分","  身体攻击","  言语攻击","  愤怒","  敌意"))

cat("\n§ 3.2.4  性别调节效应（攻击性结果变量）：\n")
run_moderation("aggress","minik")
run_moderation("aggress","hkss")

all_pcor <- bind_rows(fam_results, psych_results, attach_results, aggr_results)
t5_data  <- dplyr::select(all_pcor, 变量, pr_minik, sig_minik, pr_hkss, sig_hkss)

ft_table5 <- flextable(dplyr::select(t5_data, 变量, pr_minik, pr_hkss)) %>%
  set_header_labels(变量="变量", pr_minik="Mini-K", pr_hkss="HKSS") %>%
  add_header_row(values=c("","偏相关系数 pr（控制年龄）"), colwidths=c(1,2)) %>%
  bold(i=which(t5_data$sig_minik), j=2, part="body") %>%
  bold(i=which(t5_data$sig_hkss),  j=3, part="body") %>%
  hline(i=c(4,6,10), border=fp_border(width=0.5, style="dashed")) %>%
  ft_style(footer  = paste0("注：加粗表示通过 Bonferroni 校正（p < .002，即 .05/30）；",
                            "“是否有继父母”采用继父母人数 Stepparent_N（与原文一致）"),
           caption = "Table 5. 生活史策略与各变量的偏相关系数（控制年龄）")


# =============================================================================
# § 自动对照：复现值 vs 作者值（并排打印 δ% = 100×|作者−复现|/|作者|）
# =============================================================================
# 评级规则（参考课程《可重复检验指南》）：
#   完全一致 (δ=0) / 偏差较小 (0<δ<10%) / 偏差较大 (δ≥10%)。
#   注意：当作者值本身极小（如 .01）时，很小的绝对差也会放大 δ%，此时看绝对差更合理。
# -----------------------------------------------------------------------------

delta_pct <- function(a, m) ifelse(is.na(a) | is.na(m) | a == 0, NA, 100 * abs(a - m) / abs(a))
rate_tag  <- function(a, m) {
  d <- delta_pct(a, m)
  if (is.na(d)) return("—")
  if (abs(a - m) < 5e-3) return("完全一致")
  if (d < 10) "偏差较小" else "偏差较大"
}
cmp_rows  <- list()
add_cmp   <- function(tbl, stat, author, mine, note = "") {
  cmp_rows[[length(cmp_rows) + 1]] <<- data.frame(
    表 = tbl, 统计量 = stat, 作者值 = author, 复现值 = round(mine, 3),
    `δ%` = round(delta_pct(author, mine), 2), 评级 = rate_tag(author, mine),
    备注 = note, check.names = FALSE, stringsAsFactors = FALSE)
}

## ---- Table 1：M / SD / α ----
t1_auth_M  <- c(15.54,5.63,4.77,5.14,59.90,16.81,11.51,14.30,17.27,14.62,87.34)
t1_auth_SD <- c(7.63,3.22,2.75,2.84,19.62,6.04,4.84,5.44,8.05,11.25,13.64)
t1_auth_al <- c(.92,.91,.93,.83,.93,.83,.84,.83,.90,.84,.91)
t1_name    <- c("依恋总分","可及性","目标校正","愤怒苦恼","攻击性总分","身体攻击",
                "言语攻击","愤怒","敌意","Mini-K","HKSS")
for (i in 1:11) {
  add_cmp("T1", paste0(t1_name[i],"-M"),  t1_auth_M[i],  desc_t1$mean[i])
  add_cmp("T1", paste0(t1_name[i],"-SD"), t1_auth_SD[i], desc_t1$sd[i])
  add_cmp("T1", paste0(t1_name[i],"-α"),  t1_auth_al[i], alpha_calc[i])
}

## ---- Table 2：DSM-5 达阈值百分比 ----
t2_auth <- c(34.1,29.7,36.2,28.3,18.8,35.5,21.0,34.1,29.0,15.2,18.8,7.2,6.5,
             31.2,15.2,20.3,10.1,13.8,30.4,26.8,34.1,12.3,9.4)
for (i in 1:23) add_cmp("T2", paste0("DSM5_",sprintf("%02d",i)), t2_auth[i], threshold_pct[i])

## ---- Table 3：自伤频率（按性别，复用脚本算好的 wn/mn）----
t3_auth_w <- c(13,5,3,6,1,10,6,0,0,0,3,1,1,1,4,3,4)
t3_auth_m <- c(6,1,3,0,0,2,2,2,0,0,5,0,0,5,2,1,8)
for (i in 1:17) {
  add_cmp("T3", paste0("#",i,"-女"), t3_auth_w[i], wn[i])
  add_cmp("T3", paste0("#",i,"-男"), t3_auth_m[i], mn[i])
}

## ---- Table 4：与年龄、性别的零阶相关 ----
t4_auth_age <- c(-.21,-.15,-.29,-.11,-.39,-.36,-.32,-.26,-.32,-.43,-.15)
t4_auth_sex <- c(.06,.02,.10,.03,.09,.01,-.06,.19,.13,.01,.01)
for (i in seq_along(vars_t4)) {
  add_cmp("T4", paste0(labels_t4[i],"-Age"), t4_auth_age[i], cor(df[[vars_t4[i]]], df$age))
  add_cmp("T4", paste0(labels_t4[i],"-Sex"), t4_auth_sex[i], cor(df[[vars_t4[i]]], df$sex))
}

## ---- Table 5：偏相关（复用脚本算好的 all_pcor）----
t5_auth <- data.frame(
  var = all_pcor$变量,
  a_minik = c(.01,-.10,-.05,-.19, -.51,-.04, -.42,-.39,-.44,-.39, -.46,-.44,-.03,-.31,-.52),
  a_hkss  = c(.03,-.20,-.08,-.06, -.41,-.17, -.41,-.37,-.44,-.24, -.43,-.34,-.10,-.25,-.52))
for (i in 1:nrow(all_pcor)) {
  nt <- if (all_pcor$变量[i] %in% c("  目标校正伙伴关系","  身体攻击"))
          "★原文该格疑似笔误：HKSS可复现而Mini-K不可，详见报告" else ""
  add_cmp("T5", paste0(trimws(all_pcor$变量[i]),"-MiniK"), t5_auth$a_minik[i], all_pcor$pr_minik[i], nt)
  add_cmp("T5", paste0(trimws(all_pcor$变量[i]),"-HKSS"),  t5_auth$a_hkss[i],  all_pcor$pr_hkss[i])
}

## ---- Table 6：回归（复用脚本算好的 reg_m / reg_h / R² / F）----
t6_auth_bm <- c(0.21,-0.01,0.21,-0.04,-0.60); t6_auth_bem <- c(.27,.01,.02,-.01,-.41)
t6_auth_bh <- c(0.20,-1.63,0.07,0.01,-0.67);  t6_auth_beh <- c(.21,-.10,.01,.01,-.37)
for (i in 1:5) {
  add_cmp("T6", paste0(reg_m$预测变量[i],"-Mini-K-b"),  t6_auth_bm[i],  reg_m$b[i])
  add_cmp("T6", paste0(reg_m$预测变量[i],"-Mini-K-β"),  t6_auth_bem[i], reg_m$beta[i])
  add_cmp("T6", paste0(reg_h$预测变量[i],"-HKSS-b"),    t6_auth_bh[i],  reg_h$b[i])
  add_cmp("T6", paste0(reg_h$预测变量[i],"-HKSS-β"),    t6_auth_beh[i], reg_h$beta[i])
}
add_cmp("T6", "Mini-K R²", 0.28, r2_m); add_cmp("T6", "Mini-K F", 10.49, f_m)
add_cmp("T6", "HKSS R²",   0.26, r2_h); add_cmp("T6", "HKSS F",    9.09, f_h)

## ---- 汇总与导出 ----
cmp_all <- do.call(rbind, cmp_rows)
cat("\n================== 复现 vs 作者 自动对照 ==================\n")
print(cmp_all, row.names = FALSE)
write.csv(cmp_all, "Replication_vs_Author_comparison.csv", row.names = FALSE, fileEncoding = "UTF-8")

valid <- cmp_all[!is.na(cmp_all$`δ%`), ]
cat(sprintf("\n【汇总】共 %d 项可比；δ<10%%（完全一致+偏差较小）占 %.1f%%\n",
            nrow(valid), 100 * mean(valid$评级 != "偏差较大")))
cat("✅ 已导出 Replication_vs_Author_comparison.csv\n")


# =============================================================================
# 导出所有表格到一个 Word 文件（按论文 Table 1–6 顺序）
# =============================================================================

doc <- read_docx() %>%
  body_add_par("Table 1. 各量表均值、标准差与 Cronbach's α", style="heading 2") %>%
  body_add_flextable(ft_table1) %>% body_add_break() %>%
  body_add_par("Table 2. DSM-5 各条目达阈值的被试百分比",   style="heading 2") %>%
  body_add_flextable(ft_table2) %>% body_add_break() %>%
  body_add_par("Table 3. 各类自伤行为的发生频率（按性别）", style="heading 2") %>%
  body_add_flextable(ft_table3) %>% body_add_break() %>%
  body_add_par("Table 4. 各变量与年龄、性别的零阶相关系数", style="heading 2") %>%
  body_add_flextable(ft_table4) %>% body_add_break() %>%
  body_add_par("Table 5. 生活史策略与各变量的偏相关系数",   style="heading 2") %>%
  body_add_flextable(ft_table5) %>% body_add_break() %>%
  body_add_par("Table 6. 依恋预测生活史策略的回归分析",     style="heading 2") %>%
  body_add_flextable(ft_table6)

print(doc, target = "Hurst_Kavanagh_2017_Tables_corrected.docx")
cat("\n✅ 已导出 Hurst_Kavanagh_2017_Tables_corrected.docx\n")
