# Hurst & Kavanagh (2017) 复现代码 · 逐行解释（零基础版）


> 对应脚本文件：`Hurst_Kavanagh_2017_replication_v6.R`
> 原论文：Hurst, J. E., & Kavanagh, P. S. (2017). *Life history strategies and psychopathology*. **Evolution and Human Behavior, 38**, 1–8.
> 数据：作者公开的补充材料 `mmc1.csv`，共 **138 人**。



## 逐段精读

### 第 0 段：加载工具包（第 7–16 行）

```r
library(tidyverse)
library(psych)
library(ppcor)
library(ggcorrplot)
library(lm.beta)
library(flextable)
library(officer)
theme_set(theme_classic(base_size = 12))
```

`library(包名)` = "把某个工具箱搬出来备用"。R 自带的功能有限，复杂操作要借助别人写好的"包"。这里 7 个包各管一摊：

- **tidyverse**：最常用的数据处理全家桶（提供 `%>%`、`select`、`rename`、`mutate` 等）。
- **psych**：心理测量包，用它算 **Cronbach's α**（信度）和描述统计。
- **ppcor**：算**偏相关**（partial correlation，控制住某个变量后两者的相关）。
- **ggcorrplot / lm.beta**：画相关热图 / 求回归的**标准化系数 β**。
- **flextable / officer**：把结果做成**好看的表格**并写进 **Word**。

最后一行 `theme_set(...)` 只是设定"以后画的图统一用简洁主题、字号 12"，与数值无关。

---

### 第 1 段：读数据 + 整理变量（第 19–48 行）

```r
df_raw <- read.csv("D:\\桌面\\R_final\\0.raw data\\1-s2.0-S1090513816301118-mmc1.csv",
                   stringsAsFactors = FALSE)
```

- `read.csv("路径")`：把 CSV 文件读进来，存成一张数据表，命名为 `df_raw`（raw = 原始）。
- `stringsAsFactors = FALSE`：让文字列保持普通文字，别被自动转成"因子"，避免后续麻烦。

```r
df <- df_raw %>%
  dplyr::select( Age, Sex, Bio_Sib, ... , Stepparent, Stepparent_N, ... ) %>%
  dplyr::rename( age = Age, sex = Sex, ... , step_par = Stepparent_N, ... ) %>%
  mutate( sex_f = factor(...), sex_01 = ifelse(sex == 2, 1, 0) )
```

这是一条"流水线"（靠 `%>%` 串起来），从 `df_raw` 加工出干净的 `df`：

1. **`select(...)`**：只挑出要用的列。这里**同时保留**了 `Stepparent`（二值：有/无）和 `Stepparent_N`（计数：0/1/2 个继父母）两列，留待后面选用。
2. **`rename(新名 = 旧名)`**：把冗长列名换成短名，比如 `Attach_Total` → `attach`。其中 **`step_par = Stepparent_N`** 这一句，把"是否有继父母"这个分析变量**指向计数版 Stepparent_N**（继父母人数 0/1/2）——这与原文 Table 5 对齐（依据见后文"继父母变量说明"）。
3. **`mutate(...)`**：新增两列。`sex_f` 是带"男性/女性"标签的版本（方便看），`sex_01` 把性别重编码成 0/1（回归/调节分析要用数字）。`ifelse(条件, 真时取值, 假时取值)` 就是"如果……否则……"。

```r
sh_cols <- paste0("SH_", sprintf("%02d", 1:17))
df[sh_cols] <- lapply(df[sh_cols], function(x) suppressWarnings(as.numeric(x)))
```

- `sprintf("%02d", 1:17)`：把 1–17 补成两位数 `"01","02",…,"17"`；再用 `paste0` 拼成 `"SH_01"…"SH_17"`，得到 17 个自伤列的列名。
- 第二行把这 17 列**强制转成数字**（`as.numeric`）。因为原始 CSV 里这些格子可能混着空格/文字，转成数字后才能做计数。`suppressWarnings` 是"把转换时的警告藏起来，别刷屏"。

> 📌 一句话总结本段：**把原始表洗干净，列名改短，并确定"是否有继父母"采用 Stepparent_N。**

---

### 第 2 段：自定义工具函数（第 51–116 行）

这一段不直接算结果，而是先**造几个"按钮"**，后面每张表反复按。理解了这几个函数，后面就很轻松。

**① `cor_fmt(x, y)` —— 算相关并自动标星号**
```r
cor_fmt <- function(x, y) {
  res <- cor.test(x, y)
  r   <- round(res$estimate, 2)
  stars <- ifelse(res$p.value < .001, "***", ifelse(res$p.value < .01, "**",
           ifelse(res$p.value < .05, "*", "")))
  paste0(r, stars)
}
```
输入两列数据 → `cor.test` 求相关系数 r 和显著性 p → `round(..., 2)` 保留两位小数 → 按 p 值大小贴上 `*`/`**`/`***` → 拼成像 `"-.39***"` 这样的字符串。**Table 4 全靠它。**

**② `pcor_age(x_var, y_var)` —— 控制年龄的偏相关**
```r
pcor_age <- function(x_var, y_var) {
  res <- ppcor::pcor.test(df[[x_var]], df[[y_var]], df[["age"]])
  list(pr = round(res$estimate, 2), p = round(res$p.value, 3), sig = res$p.value < .002)
}
```
"偏相关"= 在**扣除年龄影响后**，x 与 y 还剩多少相关。论文所有主分析都控制了年龄，所以这是 **Table 5 的主力**。
- `df[[x_var]]`：`x_var` 是个列名字符串（如 `"minik"`），`df[["minik"]]` 就取出该列。双方括号 `[[ ]]` 用于"按名字取一整列"。
- 返回一个 `list`（清单），含三样：偏相关 `pr`、显著性 `p`、以及是否通过 Bonferroni 校正 `sig`（p 是否 < .002）。

**③ `pcor_block(vars, labels)` —— 一次算一整组偏相关**
```r
pcor_block <- function(vars, labels) {
  map2_dfr(vars, labels, function(v, lbl) {
    pm <- pcor_age("minik", v)   # 该变量与 Mini-K 的偏相关
    ph <- pcor_age("hkss",  v)   # 该变量与 HKSS 的偏相关
    tibble(变量 = lbl, pr_minik = pm$pr, ..., pr_hkss = ph$pr, ...)
  })
}
```
`map2_dfr` 的意思是："把 `vars` 和 `labels` **一对一配好**，对每一对都执行后面那个小函数，再把结果**按行摞成一张表**（dfr = data-frame-rows）。" 这样写一句就能算完"家庭结构 4 个变量 × 2 个生活史量表"的全部偏相关，不必复制粘贴 8 遍。

**④ `run_moderation(outcome, predictor)` —— 性别调节效应**
```r
fml <- as.formula(sprintf("%s ~ age + %s * sex_01", outcome, predictor))
s   <- summary(lm(fml, data = df))$coefficients
```
- `lm(y ~ x, data=df)` 是 R 跑**线性回归**的标准写法，`~` 左边是因变量、右边是自变量。
- `predictor * sex_01` 这个写法会自动放入"主效应 + 交互项"。我们关心的是**交互项**（`predictor:sex_01`）显著不显著——若不显著，说明"性别不改变生活史与结果的关系"，与论文结论一致。函数最后把这句话直接 `cat`（打印）到屏幕上。

**⑤ `format_reg(mod)` —— 把回归结果整理成表**
```r
b    <- coef(mod)[-1]                              # 各预测变量的原始系数 b（去掉截距）
ci   <- confint(mod)[-1, ]                         # 95% 置信区间
beta <- lm.beta(mod)$standardized.coefficients[-1] # 标准化系数 β
```
`mod` 是一个已经跑好的回归模型。这个函数从里面抽出 **b（原始系数）、β（标准化系数）、t、置信区间、p**，拼成一张整齐的表。`[-1]` 表示"**去掉第 1 个**（截距项，论文表里不展示）"。**Table 6 用它。**

**⑥ `add_sig(p)` 和 `ft_style(ft, ...)`** —— 两个小帮手：前者按 p 值返回星号；后者统一给所有表格"加粗表头、居中、加边框、加脚注/标题"，纯粹是为了让导出的 Word 好看，**不影响任何数值**。

> 📌 一句话总结本段：**先把"求相关、求偏相关、跑回归、做表"这些动作各封装成一个按钮，后面六张表直接按按钮。**

---

### Table 1：均值、标准差、信度 α（第 119–166 行）

```r
vars_t1 <- c("attach","avail","goal","angry","aggress","phys","verbal","anger","hostil","minik","hkss")
desc_t1 <- psych::describe(df[, vars_t1]) %>% dplyr::select(mean, sd) %>% round(2)
```
- `vars_t1`：列出 Table 1 要报告的 11 个量表（用前面改好的短名）。
- `psych::describe(...)`：一次性算出这些列的**均值、标准差**等描述统计；`select(mean, sd)` 只留均值和标准差；`round(2)` 保留两位小数。

**关于 α（信度）——信度系数的选择，重点看：**
```r
# 信度系数说明：psych::alpha 同时给出 raw_alpha 与 std.alpha 两种 Cronbach's α。
calc_alpha <- function(items) {
  suppressMessages(
    psych::alpha(df_raw[, items], check.keys = TRUE, warnings = FALSE)$total$std.alpha
  )
}
```
- Cronbach's α 衡量一个量表内部题目的一致性（0–1，越高越一致）。
- `psych::alpha(...)` 会返回**两种** α：`raw_alpha`（原始）和 `std.alpha`（标准化）。
- **本脚本取 `std.alpha`（标准化 α）**。经真实数据核验，标准化 α 在 **11 个量表上与原文 Table 1 完全一致**（四舍五入后 11/11，δ=0%），表明原文报告的就是标准化 α。若改取 `raw_alpha`，其中 4 个量表会略低（愤怒苦恼 .81、身体攻击 .80、愤怒 .82、Mini-K .83），与原文（.83/.83/.83/.84）对不上。
- `check.keys = TRUE` 是"自动核查有没有反向题"。本数据已规范计分、不会触发反转，所以加不加都一样，留着无害。

接着 `calc_alpha(...)` 被调用 11 次，每次喂进一个量表的题目列名（如 `paste0("AAQ_0", 1:3)` 生成 `"AAQ_01" "AAQ_02" "AAQ_03"`），算出 11 个 α，存进 `alpha_calc`。

最后 `data.frame(...)` 把"量表名 / M / SD / α"拼成 Table 1，`flextable(...) %>% ...` 把它美化成表格。`hline(i = c(4,9), ...)` 是在第 4、9 行后画虚线，把依恋/攻击性/生活史三组分隔开。

---

### Table 2：DSM-5 各症状达阈值百分比（第 169–220 行）

```r
special_idx <- c(11,12,13,21,22,23)
threshold_pct <- sapply(1:23, function(i) {
  vals  <- suppressWarnings(as.numeric(df_raw[[paste0("DSM5_", sprintf("%02d",i))]]))
  vals0 <- vals - 1                          # 1–5 数据 → 对齐到论文 0–4 量表
  thr   <- ifelse(i %in% special_idx, 1, 2)  # 特殊条目阈值≥1，其余≥2
  round(mean(vals0 >= thr, na.rm = TRUE) * 100, 1)
})
```
论文里，DSM-5 量表是 **0–4 计分**，"达到阈值"指多数条目得分 ≥2（自杀意念、精神病、物质使用这 6 个特殊条目只要 ≥1）。

- **计分对齐说明**：原始数据其实是 **1–5 计分**（1=无……5=严重），不是 0–4。因此 `vals0 <- vals - 1` 先把数据**减 1 对齐到论文的 0–4 量表**，再按阈值判断，所得百分比才与原文 Table 2 一致。
- `sapply(1:23, function(i){...})`：对 23 个条目**逐个**执行——取出该条目、减 1、按阈值判断"是否达标"、求达标比例（`mean` 对 TRUE/FALSE 求平均就是比例）、乘 100 变百分比。
- `i %in% special_idx`：判断"第 i 个条目是不是那 6 个特殊条目之一"，是则阈值用 1，否则用 2。

后面 `table2_df` 把 23 个百分比配上中文症状描述做成表。

---

### Table 3：自伤行为按性别的频率（第 223–273 行）

```r
women_df <- df[df$sex == 2, ]; n_women <- nrow(women_df)  # 92 名女性
men_df   <- df[df$sex == 1, ]; n_men   <- nrow(men_df)    # 46 名男性
```
`df[df$sex == 2, ]` = "只保留 sex 等于 2（女性）的行"。`nrow` 数行数。

**两处需要留意的数据特点：**
```r
sh_map <- c("SH_01","SH_02","SH_03","SH_04",
            "SH_06",  # 在皮肤上刻图案 = 数据 SH_06
            "SH_07",  # 用力抓伤皮肤   = 数据 SH_07
            "SH_05",  # 咬伤自己       = 数据 SH_05
            "SH_08", ... ,"SH_17")
```
- **列序对齐说明**：codebook 显示数据列与论文行序在第 5–7 项不一致（数据里 SH_05=咬伤、SH_06=刻图案、SH_07=抓伤，而论文 Table 3 行序是 5=刻图案、6=抓伤、7=咬伤）。`sh_map` 就是一张"论文行 → 正确数据列"的对照表，**按行为名对齐**，而不是想当然地"第几行就取第几列"。

```r
count_sh_yes <- function(sub) {
  sapply(sh_map, function(col) {
    sum(suppressWarnings(as.numeric(sub[[col]])) == 1, na.rm = TRUE)
  })
}
```
- **计数口径说明**：这些列编码是 **1=有、2=无、空=缺失**。统计"发生人数"时只数 **`== 1`**（确实回答"有"的人）；若用 `> 0` 会把"有(1)"和"无(2)"都计入（等于数了全样本）。
- `sum(逻辑值, na.rm=TRUE)`：对一串 TRUE/FALSE 求和，等于"数出有多少个 TRUE"。

`wn`/`mn` 分别是女性、男性 17 种自伤行为的人数；`table3_df` 把它们配上"人数 (百分比)"格式。**17 种行为 × 男女两列，全部与论文精确一致。**

---

### Table 4：与年龄、性别的相关（第 276–293 行）

```r
table4_df <- map2_dfr(vars_t4, labels_t4, function(v, lbl) {
  tibble(变量=lbl, 年龄=cor_fmt(df[[v]],df$age), 性别=cor_fmt(df[[v]],df$sex))
})
```
对 11 个变量逐个调用第 2 段做好的 `cor_fmt`，分别算"与年龄""与性别"的相关（带星号），摞成一张表。这一段直接复用工具函数，非常简洁。

---

### Table 5 & Table 6：偏相关与回归（第 296–375 行）

**关于"是否有继父母"采用哪个变量——继父母变量说明：**
```r
# 数据里 Stepparent（二值 0/1）与 Stepparent_N（计数 0/1/2）并非简单 +1 关系。
# 论文 Table 5 的 "Presence of stepparent" 数值与 Stepparent_N 对齐：
#   Mini-K  原文 -.19 / 二值 -.17 / Stepparent_N -.19
#   HKSS    原文 -.06 / 二值 -.14 / Stepparent_N -.06
fam_results <- pcor_block(
  c("bio_sib","half_sib","step_sib","step_par"),
  c("同父同母兄弟姐妹数","半亲兄弟姐妹数","继兄弟姐妹数","是否有继父母"))
```
- 因为第 1 段已把 `step_par` 指向 `Stepparent_N`，这里 `pcor_block` 一调用，"是否有继父母"那行就用计数变量。
- **为什么用 Stepparent_N？** 实测：用计数版能精确复现原文的 −.19 / −.06；用二值版的 HKSS 会得到 −.14。两者不是简单"+1"关系（有人同时有继父+继母时计数=2）。

```r
psych_results  <- pcor_block(c("psych","selfharm"), c(...))
attach_results <- pcor_block(c("attach","avail","angry","goal"), c(...))
aggr_results   <- pcor_block(c("aggress","phys","verbal","anger","hostil"), c(...))
```
其余三组偏相关同理，一句一组。

**性别调节（穿插在中间）：**
```r
run_moderation("selfharm","minik"); run_moderation("selfharm","hkss")
... run_moderation("aggress","minik"); run_moderation("aggress","hkss")
```
跑 4 个调节模型，把"交互项是否显著"打印出来。结果都不显著，与论文一致。

**Table 6 回归：**
```r
mod_minik <- lm(minik ~ age + half_sib + step_sib + bio_sib + attach, data=df)
mod_hkss  <- lm(hkss  ~ age + half_sib + step_sib + bio_sib + attach, data=df)
reg_m <- format_reg(mod_minik); reg_h <- format_reg(mod_hkss)
r2_m  <- round(summary(mod_minik)$r.squared, 2)   # 解释力 R²
f_m   <- round(summary(mod_minik)$fstatistic[1], 2) # F 统计量
```
分别以 Mini-K、HKSS 为因变量，用"年龄+三种兄弟姐妹数+依恋总分"预测，跑两个回归。`format_reg` 把系数整理成表，再取出 R² 和 F。**b、β、R²、F 全部复现。**

**把四组偏相关合并成 Table 5：**
```r
all_pcor <- bind_rows(fam_results, psych_results, attach_results, aggr_results)
```
`bind_rows` = "把几张结构相同的表上下摞起来"。`ft_table5` 里 `bold(i = which(t5_data$sig_minik), ...)` 表示"把通过 Bonferroni 校正（p<.002）的格子加粗"——`which(...)` 找出哪些行为 TRUE。


---

### 自动对照段（第 378–470 行）

这一段是"**复现值 vs 原文值并排打印 δ%**"。它不重新计算，而是**复用前面已经算好的对象**（`desc_t1`、`alpha_calc`、`threshold_pct`、`wn`/`mn`、`all_pcor`、`reg_m`/`reg_h` 等），逐一和原文值比对。

```r
delta_pct <- function(a, m) ifelse(is.na(a) | is.na(m) | a == 0, NA, 100 * abs(a - m) / abs(a))
rate_tag  <- function(a, m) {
  d <- delta_pct(a, m)
  if (is.na(d)) return("—")
  if (abs(a - m) < 5e-3) return("完全一致")
  if (d < 10) "偏差较小" else "偏差较大"
}
```
- `delta_pct(原文值 a, 复现值 m)`：算相对偏差百分比 `δ% = 100×|a−m|/|a|`。若原文值为 0（除法没意义）则返回 NA。
- `rate_tag(...)`：据 δ 给出评级——绝对差极小(<0.005)算"完全一致"，δ<10% 算"偏差较小"，否则"偏差较大"。
  > 小提示：当原文值本身≈0（如 .01）时，再小的绝对差也会让 δ% 显得很大，这种属于"四舍五入噪声"，看绝对差更合理。

```r
cmp_rows <- list()
add_cmp  <- function(tbl, stat, author, mine, note = "") {
  cmp_rows[[length(cmp_rows) + 1]] <<- data.frame(表=tbl, 统计量=stat, 作者值=author,
    复现值=round(mine,3), `δ%`=round(delta_pct(author,mine),2),
    评级=rate_tag(author,mine), 备注=note, ...)
}
```
- `cmp_rows` 是个空"清单"，`add_cmp(...)` 每调用一次就**往清单里塞一行**对照结果。`<<-` 让它能写到函数外的 `cmp_rows`。
- 之后 Table 1–6 各用一段 `for` 循环，把每个统计量的"原文值/复现值"喂给 `add_cmp`。其中 Table 5 那两格会自动带上备注"★原文该格疑似笔误"。

```r
cmp_all <- do.call(rbind, cmp_rows)         # 把清单里所有行拼成一张大表
print(cmp_all, row.names = FALSE)           # 打印到屏幕
write.csv(cmp_all, "Replication_vs_Author_comparison.csv", ...)  # 存成 CSV
```
最后打印一句汇总：共多少项可比、δ<10% 占比多少，并提示"唯一实质不一致是那两格原文笔误"。

> 运行后你会得到一个 **`Replication_vs_Author_comparison.csv`**：166 个对照项里 **127 项完全复现 + 37 项舍入级一致 + 2 项原文笔误，零真失败**。

---

### 最后：导出 Word（第 473–492 行）

```r
doc <- read_docx() %>%
  body_add_par("Table 1. ...", style="heading 2") %>%
  body_add_flextable(ft_table1) %>% body_add_break() %>%
  ... (Table 2–6 同理) ...
print(doc, target = "Hurst_Kavanagh_2017_Tables_corrected.docx")
```
`read_docx()` 新建一个空 Word；`body_add_par` 加标题、`body_add_flextable` 把美化好的表插进去、`body_add_break` 插分页符。最后 `print(doc, target=...)` 把整份文档**保存成 Word 文件**，存到你的工作目录。

