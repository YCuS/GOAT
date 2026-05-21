# GOAT

GOAT 是一个基于 PWM 的 motif 搜索小工具链。本版本回到简洁的原始流程：

1. 将 MEME 或简单 motif 格式转换为标准 PWM；
2. 仅根据 PWM 模拟得到固定阈值；
3. 用固定阈值在 FASTA 序列中滑窗搜索 motif。

本版本不再使用 GC 含量、motif clustering、relaxed/strict cluster 规则或 gap 对齐。

## 环境要求

- 支持 C++17 的编译器，例如 `g++` 或 `clang++`
- `jq`，用于读取 `config.json`
- Bash

## 编译

```bash
bash compile.sh
```

会生成三个可执行文件：

- `meme2pwm`：将 MEME 格式或简单 `Motif:` 格式转换为归一化 PWM；
- `get_thr`：根据 PWM 模拟计算 full motif 和 core region 的固定阈值；
- `search_motif`：在 FASTA 序列中搜索 motif，并输出清晰的 TSV 结果。

## 运行

先根据自己的数据修改 `config.json`，然后运行：

```bash
bash run.sh
```

默认配置使用 `examples/` 目录中的示例文件，因此编译后可以直接测试。

## 配置说明

`config.json` 中常用字段如下：

- `paths.original_pwm_file`：输入 motif 文件，可为 MEME 格式或简单 `Motif:` 格式；
- `paths.sequence_file`：输入 FASTA 文件，或包含 `.fa`、`.fasta`、`.fna` 的目录；
- `paths.results_file`：结果 TSV 文件；批量模式下也可作为输出目录；
- `paths.threshold_file`：保存阈值的 TSV 文件；
- `parameters.thr1`：full motif 阈值百分位，默认 `95.0`；
- `parameters.thr2`：core region 阈值百分位，默认 `75.0`；
- `parameters.simulation_iterations`：每个 motif 的模拟次数；
- `settings.core_length`：core window 长度，程序会选择信息量最高的一段 core。

## Motif 输入格式

简单格式如下，每一行依次表示 `A C G T`：

```text
Motif:motif_name
0.25 0.25 0.25 0.25
0.90 0.03 0.04 0.03
```

每一行会自动归一化。

## 输出格式

阈值文件列：

```text
motif_id  full_threshold  core_threshold  core_start  motif_length
```

搜索结果列：

```text
sequence_id  motif_id  strand  start  end  full_score  core_score  matched_sequence
```

坐标为输入 FASTA 上的 1-based 闭区间。`-` 链命中的 `matched_sequence` 按 motif 方向输出。

分数是 PWM 的负对数似然，越小表示匹配越好。只有 full motif 分数和 core region 分数都通过阈值时，才会输出为命中结果。
