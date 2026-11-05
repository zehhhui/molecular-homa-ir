# 分子 HOMA-IR 建模与贡献分解代码

本代码包只包含两部分：

1. **构建分子 HOMA-IR 模型**：代谢组和蛋白组分别拟合 elastic net，在外层训练集内用交叉拟合预测进行 OLS 整合；将冻结的基线模型用于同一受试者的 week 24 数据。
2. **计算 mol-HOMA-IR 降低的分子贡献率**：计算每个分子的个体贡献，再在 CMAP 和 CMAS 内分别估计两个剂量相对对照的调整后贡献，进行 bootstrap，并筛选和排序。

对应当前论文的 **705 人建模、593 人纵向分析、不调整遗传 PC** 版本。完整分析得到 CMAP 85 个、CMAS 80 个、两试验共有 52 个分子；这些数量是分析结果，程序不会强制得到这些数量。

## 运行

安装依赖后，从代码包根目录运行：

```bash
conda env create -f environment.yml
conda activate molecular-homa-ir
bash tests/run_smoke.sh
```

示例使用完全模拟的数据，并缩小交叉验证和 bootstrap 次数，用来检查流程是否可运行。示例的结果没有临床含义，也不应与论文数值比较。

正式数据的运行顺序及参数见 [README.md](README.md)。输入格式见 [docs/data_requirements.md](docs/data_requirements.md)。真实数据、个体标识、原始 fold 分配、拟合模型和个体结果均不在发布包内。

## 贡献率的含义

先在全部分子上分解每个人的分数变化，再拟合：

```text
分子贡献 ~ 治疗组 + 基线 mol-HOMA-IR + 性别 + 年龄
```

分别计算 CG32 对照差异和 CG48 对照差异。贡献率的分母是该对比中**全部模型分子的降低方向贡献之和**，不只是 148 个候选、85/80 个入选分子或图中前 10 个分子。

随后要求每次重复都至少有一折入选，并在对应试验的两个剂量对照差异中均呈降低方向。这个步骤没有显著性筛选。最终按两个剂量的平均贡献率排序。

这是治疗与对照的差异分解，不能自动解释为每个分子在治疗组内显著改变、介导疗效或独立于脂肪分布变化。

## 许可与引用

代码采用 [MIT 许可](LICENSE)，版权署名为 **2026 Zehui Cao and contributors**。版本为 **1.0.0**，软件创建者为 **Zehui Cao（Fudan University）**，GitHub 首次发布日期为 **2026-11-05**。`CITATION.cff` 和 `.zenodo.json` 使用一致的引用信息。

- [GitHub 仓库](https://github.com/zehhhui/molecular-homa-ir)
- [Zenodo 归档](https://doi.org/10.5281/zenodo.23991482)：DOI 已预留，归档记录待正式发布。
- [Code Ocean capsule](https://codeocean.com/capsule/5336449/tree)：已创建，公开访问及平台上的可复现运行尚待核实。

Code Ocean 的运行和环境说明见 [docs/code_ocean.md](docs/code_ocean.md)。其公开演示使用合成数据。请勿同时上传项目的真实输入数据、`results/`、个体预测或模型 RDS。更完整的实现来源和版本界限见 [docs/provenance.md](docs/provenance.md)。
