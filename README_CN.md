# genodiff — 自研基因组比对变异检测工具

[![Rust](https://img.shields.io/badge/language-Rust-DEA584.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-45%20passed-brightgreen.svg)](#测试)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20Windows-lightgrey.svg)](#安装)

**零依赖、单二进制的单倍体基因组比对与变异检测工具。**

genodiff 输入两个（或多个）组装基因组 FASTA，输出包含 SNP、小 InDel、
结构变异（DEL / INS / INV / DUP / BND）的标准 VCF 4.2 文件——
一个纯 Rust 静态二进制，**无任何外部依赖**。

> 📖 English documentation: [README.md](README.md)
> 📖 完整技术文档（底层逻辑 / 程序图 / 创新点）: [docs/PIPELINE.md](docs/PIPELINE.md)

---

## 为什么选择 genodiff？

| | genodiff | minimap2 + Syri | MUM&Co |
|---|---|---|---|
| 外部依赖 | **0** | ≥4 | MUMmer + Python |
| 一条命令 FASTA→VCF | **是** | 多步流程 | 部分 |
| SV 类型 | DEL/INS/INV/DUP/BND | DEL/INS/INV/DUP/TRA | 部分 |
| 直接 VCF（已左规范化） | **是** | 自有格式 + VCF | VCF |
| 内置模拟器 + 真值 VCF | **是** | 无 | 无 |
| 内置多样本合并 | **是** | 需 bcftools | 无 |
| 安装 | `cargo build` | 多工具安装 | MUMmer 编译 |

### 核心创新

1. **单二进制零依赖**：纯 Rust 静态二进制 FASTA→VCF 全流程
2. **组装比较专用**：绕开读段比对参数体系，共线性锚点链为一等公民
3. **patience 断链容忍**：800 Mb 真实数据无效断链 1.03×10⁷ → 5.2×10⁵
4. **对称间隙倒位判据**：桥接任意长度共享重复侧翼
5. **种子重叠滑移 indel**：恢复常规间隙比对遗漏的同聚物/短重复 indel
6. **重复吸引子过滤**：自动识别并剔除旁系同源比对
7. **有界内存流式设计**：800 Mb × 800 Mb 全流程 117 s / ~4 GB（48 线程）
8. **内建验证闭环**：模拟器 + 真值 VCF + 45 项测试

详细说明见 [docs/PIPELINE.md](docs/PIPELINE.md)。

## 安装

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"
git clone https://github.com/<your-username>/genodiff.git
cd genodiff && cargo build --release && cargo test --release
```

## 快速开始

```bash
G=target/release/genodiff

# 一条命令 FASTA → VCF
$G call ref.fa query.fa -o out.vcf --threads 48 --sample query

# 或分步：索引 + PAF
$G index ref.fa ref.gdx
$G map ref.gdx query.fa --threads 48 > out.paf

# 批量（18 基因组 → 多样本 VCF）
bash scripts/run_all.sh

# 可视化
python scripts/plot_syntenic_dotplot.py out.paf -o dotplot.png --top 8
python scripts/plot_syntenic_linear.py  out.paf -o synteny.svg --top 8
python scripts/plot_syntenic_circos.py  out.paf -o circos.png  --top 8
```

## 精度验证

模拟基因组（5 Mb，5530 变异）：**SNP 100%/100%，InDel 100%/100%，
四类 SV 100% 召回**。真实数据（800 Mb × 800 Mb）：call 117 s，
SNP=597 万，SV 合计 ~55,000。

## 已知限制

1. SV 断点精度链端级（~数百 bp）
2. 贪心链式在多拷贝区域可能选次优拷贝
3. 复杂重组（>2kb 间隙且无共享侧翼）仅报链断点
4. 着丝粒/高重复区种子层不可比对（所有方法的共同边界）

## 许可

[MIT](LICENSE)
