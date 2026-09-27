
import importlib
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

import sccoda.datasets as scd

celltype_order = [
    "CD4NC",
    "CD4ET",
    "CD8NC",
    "CD8ET",
    "Prolif",
    "NK",
    "NKR",
    "BIN",
    "BMem",
    "ABC",
    "Plasma",
    "MonoC",
    "MonoNCI",
    "MonoNC",
    "Neu",
    "LDG",
    "cDC",
    "pDC",
    "Mega"
]

celltype_order = [
    "CD4T",
    "CD8T",
    "NK",
    "B",
    "Mono",
    "Neu",
    "Mega"
]

composition_file = (
"~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/"
"immune_cell_composition_matrix_majorCelltype.csv"
)
composition_file = (
"~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/"
"immune_cell_composition_matrix_Celltype_ABC.csv"
)

clinical_file = (
"~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/"
"clinical_data_Covariates.txt"
)

out_dir = (
"~/scRNA/B/10HC11LN/DEG/DEG_pseudoBulk_Deseq2/scCODA/"
)

composition = pd.read_csv(
    composition_file,
    index_col=0
)
print(composition.head())
print(composition.shape)
# clinical data
clinical = pd.read_csv(
    clinical_file,
    sep="\t"
)
#clinical.index = clinical["sample"]
print(clinical.head())

###6. 合并clinical和cell count
composition.index.name="sample"
data_df = composition.merge(
    clinical,
    left_index=True,
    right_on="sample",
    how="left"
)
data_df.index = data_df["sample"]
#下面是之前只保留一个Group的整合代码##
#data_df = pd.concat(
#    [
#        clinical[["Group"]],
#        composition
#    ],
#    axis=1
#)
print(data_df.head())
print(data_df.columns)
#7. 设置Group reference
data_df["Group"] = pd.Categorical(
    data_df["Group"],
    categories=[
        "HC",
        "LN"
    ]
)
#8. 转换为scCODA AnnData对象
sccoda_data = dat.from_pandas(
    data_df,
    covariate_columns=[
        "sample",
        "Group",
        "Age",
        "Gender",
        "Batch"
    ]
)
sccoda_data = sccoda_data[
    :,
    celltype_order
]
print(sccoda_data)
####画细胞比例box_Plot
viz.boxplots(sccoda_data,
             figsize=(8, 3.5),
             feature_name="Group")
plt.savefig(out_dir+"scCODA_boxplot_Major_celltype.pdf")
########画堆积柱状图
# Stacked barplot for each sample
viz.stacked_barplot(sccoda_data,  figsize=(8, 6), feature_name="sample")
plt.savefig(out_dir+"scCODA_Stacked_boxplot_sample.pdf")

# Stacked barplot for the levels of "Condition"
viz.stacked_barplot(sccoda_data, feature_name="Group")
plt.savefig(out_dir+"scCODA_Stacked_boxplot_Group.pdf")
#############Grouped boxplots
viz.boxplots(
    sccoda_data,
    feature_name="Group",
    plot_facets=False,
    y_scale="relative",##显示的是proportion
    add_dots=False,
)
plt.savefig(out_dir+"scCODA_boxplot.pdf")##与上面单独画的boxplot一致

# Grouped boxplots. Facets, log scale, added dots and custom color palette.
viz.boxplots(
    sccoda_data,
    feature_name="Group",
    plot_facets=True,
    y_scale="log",
    add_dots=True,
    cmap="Reds",
)

fig = plt.gcf()

fig.set_size_inches(8,12)

fig.savefig(
    out_dir+"scCODA_Major_celltype_order_Facet_boxplot.pdf",
    bbox_inches="tight"
)

plt.close(fig)
###################设定细胞排列顺序和颜色#########

viz.boxplots(
    sccoda_data,
    feature_name="Group",
    figsize=(12,4),
    y_scale="relative"
)
plt.savefig(
    out_dir+"scCODA_boxplot_allcelltype_order.pdf",
    bbox_inches="tight"
)

#9.运行CompositionalAnalysis（官方方法）
analysis = mod.CompositionalAnalysis(
    sccoda_data,
    formula="Group"
)
#10.运行MCMC
results = analysis.sample_hmc(
    num_results=20000
)
#11.查看结果
results.summary()
#12.保存结果
effect_df = results.effect_df

effect_df.to_csv(
    out_dir+"scCODA_effect_df.csv"
)
credible = results.credible_effects()

credible.to_csv(
    out_dir+"scCODA_credible_effects.csv"
)
#保存完整summary文本
with open(
    out_dir+"scCODA_summary_extended.txt",
    "w"
) as f:
    f.write(
        str(
            results.summary_extended()
        )
    )
#####13.画图
library_data = results.effect_df
effect = results.effect_df.reset_index()
effect_plot = effect[
    effect["Covariate"]=="Group[T.LN]"
]


plt.figure(figsize=(6,4))

plt.barh(
    effect_plot["Cell Type"],
    effect_plot["Final Parameter"]
)

plt.axvline(
    0,
    linestyle="--"
)

plt.xlabel(
    "scCODA effect size"
)

plt.tight_layout()

plt.savefig(
    out_dir+"scCODA_effect_plot.pdf"
)