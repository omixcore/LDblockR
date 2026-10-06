# LDblockR 0.0.1 使用说明

本文件对应源码包 `LDblockR_0.0.1_source.zip` 和安装包
`LDblockR_0.0.1.tar.gz`。项目主页为
<https://gitee.com/anhuikylin/LDblockR>；GitHub 镜像为
<https://github.com/anhuikylin/LDblockR>。包为纯 R 实现，要求 R ≥ 4.1.0，不依赖 Java、
TASSEL、PLINK 或额外 CRAN 包。

## 安装

```r
install.packages("LDblockR_0.0.1.tar.gz", repos = NULL, type = "source")
library(LDblockR)
packageVersion("LDblockR")
```

## 内置示例与来源

所有示例文件随包安装。官方 TASSEL/rTASSEL 五文件可这样取得：

```r
paths <- example_data("tassel")
paths
stopifnot(all(file.exists(paths)))
file.show(paths[["SOURCE"]])
file.show(paths[["PROVENANCE"]])
```

`mdp_genotype.hmp.txt`、`mdp_phenotype.txt`、`mdp_kinship.txt`、
`mdp_population_structure.txt` 和 `mdp_traits.txt` 是
maize-genetics/rTASSEL 固定提交
`e8cbfcce47422319b5bfe8376cf5550ebd3b365e` 的原始文件副本。原始下载地址、
字节数和 SHA-256 在 `PROVENANCE.json` 中；上游许可通知在
`LICENSE.rTASSEL.txt` 中。其它 VCF、HapMap、矩阵、GFF3、PLINK 和
LDBlockShow 参考输出均位于 `system.file("extdata", package="LDblockR")`。

内置区间 GWAS 示例可以直接这样读取：

```r
regional <- example_data("regional")
stopifnot(all(file.exists(regional[names(regional) != "plink"])))
gwas <- read_gwas(regional[["regional_gwas"]])
regional[["regional_gwas"]]
```

`regional_gwas.tsv` 随包安装；也可以使用 `regional[["regional_gwas"]]` 这个绝对路径。
如果当前安装包报告示例数据不完整，不能用任意占位文件名代替；应重新安装完整源码包，或传入
自己电脑上已经存在的真实输入文件路径。

## 混合线性模型 GWAS

下面的流程使用 `mdp_phenotype.txt` 的 `EarHT`，保留两个地点的重复观测，并将
`location`、`Q1`、`Q2`、`Q3` 作为固定效应：

```r
library(LDblockR)
paths <- example_data("tassel")
pheno <- read_tassel_phenotype(paths[["mdp_phenotype"]])

gwas <- gwas_mlm(
  genotype = paths[["mdp_genotype"]],
  phenotype = pheno,
  trait = "EarHT",
  covariates = c("location", "Q1", "Q2", "Q3"),
  replicate = "expand",     # 保留 A/B 地点重复观测
  min_maf = 0.05,
  max_missing = 0.20,
  chunk_size = 256L,
  verbose = TRUE
)

head(gwas[order(gwas$p), c("id", "chr", "pos", "p", "PVE", "model_PVE")])
attr(gwas, "model")$model_PVE
```

## Manhattan 图与 Q-Q 图

两类诊断图都直接使用同一份 `gwas_mlm()` 结果，不重新计算关联统计量，因此图形
与区域 LD 主图完全一致。输入也可以是 `read_gwas()` 能识别的 TSV 文件或普通
数据框。函数不依赖 `ggplot2` 或其他外部绘图库。

```r
cutline <- -log10(0.05 / nrow(gwas))
lead <- gwas$id[which.max(gwas$logp)]

manhattan <- plot_manhattan(
  gwas,
  cutline = cutline,                 # -log10(P)；NULL 默认 Bonferroni
  color_by = "significance",        # 也可为 chromosome 或 pve
  point_size_by = "PVE",             # 可改为 none
  highlight = lead, label_top = 5,
  title = "Ear height mixed-model GWAS",
  draw = FALSE
)

qq <- plot_qq(
  gwas, confidence = 0.95, highlight = lead,
  title = "Ear height mixed-model GWAS Q-Q plot",
  draw = FALSE
)

save_gwas_plot(manhattan, "Figure_EarHT_Manhattan.pdf", width = 8.5, height = 5.4)
save_gwas_plot(manhattan, "Figure_EarHT_Manhattan.svg", width = 8.5, height = 5.4)
save_gwas_plot(manhattan, "Figure_EarHT_Manhattan.png", width = 8.5, height = 5.4, dpi = 400)
save_gwas_plot(qq, "Figure_EarHT_QQ.pdf", width = 6.8, height = 6.2)
save_gwas_plot(qq, "Figure_EarHT_QQ.svg", width = 6.8, height = 6.2)
save_gwas_plot(qq, "Figure_EarHT_QQ.png", width = 6.8, height = 6.2, dpi = 400)
```

### Manhattan 图参数

- `color_by = "chromosome"`：按染色体交替着色，适合标准全基因组图。
- `color_by = "significance"`：超过 `cutline` 的点为红色，其余为蓝色，和区域
  LD 主图一致。
- `color_by = "pve"`：连续颜色编码 marker-level `PVE`；`point_size_by = "PVE"`
  还会用点大小编码 PVE，并自动添加 PVE 图例。
- `highlight` 可传 SNP ID、`chr:position` 字符串，或包含 ID/chr/pos 的数据框；
  `label_top` 可标注最显著的前若干 SNP。

`plot_manhattan()` 的参数设计兼容原有 `cutline`/`palette`，并提供 CMplot 风格的
逐染色体配色和版式控制：

```r
# 1) 为 gwas 中实际出现的每条染色体生成颜色和标签。
#    含有 10 条染色体的 gwas 不能只传 3 个颜色或 3 个标签。
chromosomes <- unique(as.character(gwas$chr))
chr_col <- setNames(grDevices::hcl.colors(length(chromosomes), palette = "Dark 3"),
                     chromosomes)
chr_labels <- setNames(
  ifelse(grepl("^chr", chromosomes, ignore.case = TRUE), chromosomes,
         paste0("Chr ", chromosomes)),
  chromosomes
)
m1 <- plot_manhattan(
  gwas,
  threshold = c(-log10(0.05 / nrow(gwas)), 4),
  suggestive = 2,
  chr_colors = chr_col,
  color_cycle = FALSE,
  chr_labels = chr_labels,
  color_significant = TRUE,
  point_size_by = "PVE", signal_cex = 1.35,
  threshold_col = c("#D73027", "#762A83"),
  threshold_lty = c(2, 3), threshold_lwd = c(1.2, 0.9),
  draw = FALSE
)

# 2) 两色循环；染色体多于颜色数时自动循环
m2 <- plot_manhattan(gwas, chr_colors = c("#2C7FB8", "#7B3294"),
                     color_cycle = TRUE, draw = FALSE)
```

常用参数含义：

- `threshold`/`suggestive`：一个或多个 `-log10(P)` 阈值；`threshold` 会覆盖
  旧参数 `cutline`。对应的颜色、线型和线宽可用 `threshold_col`、`threshold_lty`、
  `threshold_lwd` 与 `suggestive_col` 等设置。
- `chr_colors`：命名向量按染色体标签匹配；不命名时按染色体排序使用。设
  `color_cycle = FALSE` 可强制检查颜色是否完整，避免颜色错配。
- `chr_labels`、`gap_ratio`/`gap_bp`/`gap_mode`：分别控制横轴标签和染色体间空隙。
  `gap_mode = "fixed"`（默认）在所有染色体边界使用同一个间隔；如果要明确指定
  固定距离，可设置 `gap_bp`，其单位与 `pos` 相同。例如，物理坐标为 bp 时，
  `gap_bp = 5000000` 表示每条染色体之间固定间隔 5 Mb。`gap_mode = "relative"`
  仅在需要按前一条染色体宽度缩放间隔时使用。
- `chr_lengths`：传入命名或按顺序排列的物理染色体长度（bp），使横坐标中每条
  染色体的宽度按真实长度比例显示。默认 `width_mode = "length"` 使用每条染色体
  的最大观测坐标；`width_mode = "observed"` 按标记跨度，`"equal"` 等宽。
- `pch`、`point_size`、`point_cex`、`signal_cex`、`alpha`：控制点形状、大小、
  显著点放大和透明度；`point_size_by = "PVE"` 仍可用点大小表达 marker-level PVE。
- `ylim`/`xlim`、`show_grid`、`axis_cex`、`label_cex`、`title_cex`：控制投稿版
  坐标范围、网格和字体。`show_threshold_legend = FALSE` 可隐藏阈值线图例。

如果不确定颜色与染色体的对应关系，优先使用带名字的 `chr_colors`；如果使用不命名
颜色并设置 `color_cycle = FALSE`，颜色数量少于染色体数会直接报错。

例如使用参考基因组长度和独立染色体配色：

```r
# 可运行示例：用当前 GWAS 中每条染色体的最大观测坐标生成完整向量。
# 正式分析应替换为真实参考基因组长度。
chromosomes <- unique(as.character(gwas$chr))
chr_lengths <- tapply(gwas$pos, as.character(gwas$chr), max, na.rm = TRUE)
chr_lengths <- chr_lengths[chromosomes]
m_length <- plot_manhattan(
  gwas, chr_lengths = chr_lengths, width_mode = "length",
  chr_colors = chr_col,
  color_cycle = FALSE, gap_bp = 5000000, draw = FALSE
)
```

如果使用参考基因组长度，必须为 `gwas$chr` 中的每条染色体提供一个值。如果只绘制 1、2、3
号染色体，应先对子集后的 `gwas` 绘图，不能把三组参数直接传给十染色体数据。

### Q-Q 图解释

Q-Q 图横轴是理论 `-log10(P)`，纵轴是观测 `-log10(P)`；紫色透明带是指定置信度
下的 order-statistic 置信包络，灰/黑色对角线是完全符合均匀零分布的参考线。图例
中的 `lambda` 使用 1 df 卡方变换后的中位数比值计算；它只用于诊断整体膨胀，不替代
混合模型的 kinship 校正。

### 与区域 LD 主图联用

三种图使用同一 `gwas` 对象即可：

```r
region_info <- gwas_ld_region(gwas, cutline = cutline, flank = 250000,
                              min_width = 100000)
gwas_region <- region_info$selected
```

将 `gwas_region` 传给 `plot_ld(..., gwas = gwas_region)`；这样 Manhattan/Q-Q 是
全基因组诊断，LD 图是阈值点所在局部区间，三者不会混用 MAF 代替 GWAS 统计量。

模型为 `y = X beta + Zu + e`，`u ~ N(0, K sigma_g^2)`。实现采用与 GAPIT
P3D/MLM 相同的统计原理，但代码独立：只对零模型估计一次
`delta = sigma_e^2 / sigma_g^2`，只做一次 kinship 特征分解，然后以 BLAS 友好的
marker 块进行白化、固定效应残差化和一自由度 GLS 检验。

- `PVE`：单个 SNP 的部分表型方差解释率（加 SNP 后加权残差平方和的相对下降）。
- `model_PVE`：零模型的方差组分估计 `1/(1+delta)`，在结果列中重复提供，也保存在
  `attr(gwas, "model")$model_PVE`。
- `replicate = "mean"`：按 Taxa 平均重复表型；默认 `"expand"` 保留重复观测。
- `kinship = paths[["mdp_kinship"]]`：改用官方 TASSEL kinship；样本按 Taxa 自动取交集。

若输入是矩阵，要求样本在行、marker 在列、编码为 `0/1/2/NA`，并用 `map=`
提供 `chr`、`pos`、`id`：

```r
gwas2 <- gwas_mlm(G, phenotype, trait = "height", map = snp_map,
                  kinship = K, covariates = PCs, n_pc = 0)
```

## 投稿版 LDheatmap

先用超过阈值线的 GWAS 点定义区域，再将同一局部 GWAS 表传给热图。这样不会把
全基因组中不属于 LD 矩阵的点画到顶部：

```r
cutline <- -log10(0.05 / nrow(gwas))
region_info <- gwas_ld_region(gwas, cutline = cutline, flank = 250000,
                              min_width = 100000)
region <- region_info$region
lead <- region_info$lead_id
z <- gwas[match(lead, gwas$id), , drop = FALSE]
gwas_region <- region_info$selected

x <- read_hapmap(paths[["mdp_genotype"]], region = region,
                 min_maf = 0.05, max_missing = 0.20, quiet = TRUE)
ld <- ld_compute(x, measure = "r2", r2_method = "dosage")
blocks <- detect_ld_blocks(ld, method = "solid_spine", metric = "r2", spine_cut = 0.8)
tags <- select_tag_snps(ld, threshold = 0.8, blocks = blocks)

p <- plot_ld(ld, metric = "r2", gwas = gwas_region, lead = lead,
              blocks = blocks, tags = tags,
              special = data.frame(id = lead, pos = z$pos, label = "Lead SNP"),
              cutline = cutline, show_connectors = TRUE,
              show_snp_connectors = TRUE,
              heatmap_colors = c("#FFFDF2", "#FDBB55", "#B40426"),
              palette = "publication", position_scale = "index",
              show_values = TRUE,
              show_maf = FALSE, gwas_color_by = "significance",
              show_snp_labels = TRUE, snp_label_type = "id_position",
              max_snp_labels = Inf,
              title = "Ear height mixed-model GWAS and regional LD", draw = FALSE)
save_ld_plot(p, "Figure_EarHT_LD.pdf", width = 8.5, height = 7.2)
save_ld_plot(p, "Figure_EarHT_LD.svg", width = 8.5, height = 7.2)
```

`gabriel`、`solid_spine`、`strong` 等数据驱动方法在阈值严格或区域较小时可能返回
0 行区块，这是正常结果。计算单倍型前不要直接访问 `blocks$snps[1]`，应先判断
`nrow(blocks) > 0`；如果需要复现内置区域示例，可读取 `regional[["fixed_blocks"]]`
并使用 `detect_ld_blocks(ld, method = "fixed", fixed = fixed)`，且固定区间的染色体和坐标必须
与 `ld$data$variants` 一致。

这里的 `gwas` 与 `paths[["mdp_genotype"]]` 必须来自同一基因型面板；如果使用内置区域
示例的 `regional_gwas.tsv`，应同时使用 `example_data("regional")` 返回的
`regional[["hapmap"]]`，不要与 TASSEL 的 `mdp_genotype.hmp.txt` 混用。区域读取会兼容
`1` 与 `chr1` 这类染色体标签差异，但不会生成基因型文件中不存在的 marker。

布局与参考图对应：

1. 顶部区间 GWAS 的 `-log10(P)` 轨道、红色虚线阈值、红色显著点、蓝色非显著点；有 `PVE` 时点
   大小编码 PVE，并给出 PVE 图例。
2. 中间可选彩色区域/基因条、基因模型、基因名、MAF 轨道以及 `Start / chr Region / End` 标注。
3. SNP 标注轨道与热图共用同一组 SNP 索引顺序，默认逐个显示 `SNP ID` 和物理位置；每个 marker
   在热图顶部基线都有对应短刻度，并用绿色浅色虚线直接连接到热图顶部，因此可以直接追踪 SNP
   在哪个热块中。
4. Block 轨道位于 LD 热图正上方，Block 名称和 V 形边界按同一组 SNP 索引顺序绘制。
5. 底部倒三角 LD 热图（SNP 共线基准位于顶部），默认采用等间距 SNP 索引并恢复紧凑的正方形单元格布局；每个 Block 对应的 LD 子三角会用同色边界勾出，并在热块内部标注 `Block1`、`Block2` 等名称。保留黑/白单元格边界、两位小数格内数值、tagSNP
   指示线、特殊 SNP 标签和连续 `R² color key`，但不绘制热图外部矩形框或总三角外围线。
   GWAS 点通过加粗连线连接到同一热图基准 marker，基因组位置坐标移到上方 GWAS
   轨道，热图下方不再绘制横坐标轴。

`plot_ld()` 会优先按 marker ID、再按 `chr + position` 精确匹配 GWAS 与 LD 变异；
未进入 LD 矩阵的行不会被硬塞到图中，而是给出警告并丢弃。若要关闭连接线，可
设置 `show_connectors = FALSE`。

MAF 轨道默认关闭；旧式频率轨道可显式使用 `show_maf = TRUE`。有真实同组装版本
的 GFF3 时，传入 `genes = "annotation.gff3"` 即可绘制基因模型；没有注释时仍可
保留顶部 GWAS 和底部 LD 两个核心轨道。

若希望直接使用与参考图一致的综合显示功能，可调用 `plot_ld_region()`。该函数默认
打开 MAF、关键 SNP 图例和 SNP 到热图的虚线引导，并按 lead SNP 的 `r²` 给区域 GWAS
点着色；它仍然返回标准的 `ld_plot` 对象，因此可以继续用 `save_ld_plot()` 输出
PDF、SVG、PNG 或 TIFF。

```r
p_region <- plot_ld_region(
  ld, gwas = gwas_region, genes = regional[["gff3"]],
  blocks = blocks, tags = tags, lead = lead,
  cutline = cutline, draw = FALSE
)
save_ld_plot(p_region, "Figure_integrated_regional_LD.pdf", width = 8.5, height = 8)
save_ld_plot(p_region, "Figure_integrated_regional_LD.png", width = 8.5, height = 8, dpi = 400)
```

关键 SNP 轨道中，紫色菱形为 lead SNP，橙色圆点为 tag SNP，蓝色三角为 block
boundary，红色星号为 `special` 变异；可用 `key_snp_colors` 修改这四类颜色。

SNP 标注轨道默认最多显示 120 个标签；大区域可使用
`max_snp_labels = Inf` 显示全部标签，或使用 `snp_label_type = "id"`、
`snp_label_type = "position"` 简化标签；`show_snp_labels = FALSE` 可关闭整个轨道，
`show_snp_connectors = FALSE` 可关闭 SNP 到热块的虚线。

热图颜色可用 `heatmap_colors` 自由指定 2 个或更多颜色，例如
`heatmap_colors = c("white", "skyblue", "navy")`；它会覆盖 `palette`，并同步更新热图与
下方的 `R² color key`。不提供 `heatmap_colors` 时，原有的 `palette` 参数保持有效。

如需按照真实碱基距离拉伸 SNP，可显式设置 `position_scale = "physical"`；此时
热块不保证为正方形。

## 一键示例脚本

源码包中的 `scripts/reproduce_gwas_example.R` 会输出 GWAS TSV、`sessionInfo()`
以及 PDF/SVG/PNG 主图：

```bash
Rscript scripts/reproduce_gwas_example.R EarHT_results
```

如有同一组装版本 GFF3，可作为第二个参数传入：

```bash
Rscript scripts/reproduce_gwas_example.R EarHT_results annotation.gff3
```

参考图风格的综合区域图可用内置示例一键复现：

```bash
Rscript scripts/reproduce_integrated_region.R integrated_region_results
```

论文 Figure 2 的复现区分两类数据来源：区域示范图使用仓库内 `inst/extdata/` 的合成基因型、GWAS
和注释文件；数值验证使用 861 对合成 SNP 与归档的 LDBlockShow 1.41 结果进行比较。官方 TASSEL
玉米示例用于独立的玉米 QC 和区域分析示例。因此不需要寻找未随包提供的
`population.vcf.gz` 或 `maize507.hmp.txt.gz` 等外部文件。从源码目录运行：

```bash
Rscript scripts/reproduce_manuscript.R LDblockR_results
```

Figure 2 输出在 `LDblockR_results/figures/Figure_2_reference_results.pdf`、`.svg` 和 `.png`；
861 对 SNP 的数值验证表输出在 `LDblockR_results/results/` 下的 `R_numerical_validation.tsv`
和 `R_synthetic_pair_validation.tsv`。

## 复现与投稿核验

正式投稿前保存：`gwas_mlm_result` 原始 TSV、`attr(gwas,"model")`、完整命令、
`sessionInfo()`、安装日志、`R CMD build`/`R CMD check` 日志，以及 LDBlockShow
861 对 SNP 比较的原始输出、软件版本与参数。固定源码版本应同时记录 Git commit、
安装包 SHA-256 和 `PROVENANCE.json`。
