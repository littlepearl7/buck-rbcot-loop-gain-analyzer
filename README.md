# BUCK 变换器 RBCOT 环路增益分析工具

## 概述

本工具用于分析 BUCK 变换器在 **RBCOT（Ripple-Based Constant On-Time）** 控制下的环路增益特性，实现了论文 *"Accurate Loop Gain Model of Ripple-Based Constant on-time Controlled Buck Converters"* 中的三种补偿模型：

- **ESR 纹波补偿模型**
- **固定斜坡补偿模型**
- **SW 滤波器补偿模型**

## 文件结构

```
BUCK/
├── model_val_acc2approx.py          # ESR近似模型验证（主入口）
├── model_val_acc2imp.py             # 完整三模型GUI分析器
├── rbcot_formulas.py                # 核心公式库（传递函数/Bode分析）
```

## 核心文件说明

| 文件 | 功能 |
|------|------|
| `rbcot_formulas.py` | 传递函数实现、Padé近似、Bode分析工具（计算引擎） |
| `model_val_acc2imp.py` | 完整GUI，支持三种补偿模型切换与对比 |
| `model_val_acc2approx.py` | ESR近似模型验证入口（继承主GUI，重载ESR计算链路） |

## `model_val_acc2approx.py` 详细说明

`model_val_acc2approx.py` 是 **ESR近似公式验证** 的专用入口脚本，重点不是三模型切换，而是将ESR模型改为论文近似解并与Bode数值结果对照显示。

- 基于 `model_val_acc2imp.py` 的 `BUCKBodeAnalyzer` 继承实现（复用UI框架）
- 重载 `calculate_esr_model()`，调用 `rbcot_formulas.py` 的 `esr_approx_bode()`
- 在结果区同时显示：
  - `Fc`（Bode 0dB插值得到的实际穿越频率）
  - `Eq.(13)`（论文公式(13)近似交叉频率）
  - `PM`（Bode相位计算值）
  - `Eq.(14)`（论文公式(14)近似相位裕度）
  - `fz`（ESR零点频率）
- 若无0dB穿越，会退化显示峰值增益频点和近似公式结果，便于判断设计偏差

**适用场景**：快速验证“近似设计公式”和“频域扫频结果”是否一致时，优先运行该脚本。

## 快速开始

### 安装依赖

```bash
pip install numpy matplotlib control pillow
```

### 运行ESR近似模型验证

```bash
python model_val_acc2approx.py
```

### 运行完整三模型分析

```bash
python model_val_acc2imp.py
```

## 模型说明

### 1. ESR纹波补偿模型

利用输出电容ESR产生的纹波作为PWM斜坡，是最自然的补偿方式。适用于ESR较大的电解电容场景。

**核心传递函数：**

```
F_loop(s) = A₀·(1+s/ωz) / [(1+s/(Qp1·ωp1)+s²/ωp1²)·(1+s/(Qp2·ωp2)+s²/ωp2²)]
```

### 2. 固定斜坡补偿模型

通过外部电路产生固定斜率斜坡，适用于陶瓷电容等小ESR场景。无需电感电流检测，但高占空比下性能下降。

**核心传递函数：**

```
F_loop(s) = A₀_s·(1+s/(Qxs·ωxs)+s²/ωxs²) / [(1+s/(Qps1·ωps1)+s²/ωps1²)·(1+s/(Qps2·ωps2)+s²/ωps2²)]
```

### 3. SW滤波器补偿模型

利用SW节点RC滤波器产生斜坡，兼具采样和补偿功能。复杂度最高，设计需谨慎。

**核心传递函数：**

```
F_loop(s) = A₀_r·(1+s/ωzr) / [(1+s/(Qpr1·ωpr1)+s²/ωpr1²)·(1+s/(Qpr2·ωpr2)+s²/ωpr2²)]
```

## 设计指导

| 模型 | 优化参数 | 推荐范围 |
|------|----------|----------|
| ESR补偿 | Rc | D<0.2: Rc = Tsw/(π·Co)<br>D>0.5: Rc = 4Tsw/Co |
| 固定斜坡 | Vslope | Vslope_opt = 2.3·Tsw²·Vout/(π²·L·C) |
| SW滤波器 | Rx·Cx | π·L·C/(5Tsw) < Rx·Cx < L·C/√((Tsw/π)²+(Ton/π)²) |

## 参数默认值（论文Fig.3基准）

| 参数 | 默认值 | 说明 |
|------|--------|------|
| Vin | 12 V | 输入电压 |
| Vout | 1.2 V | 输出电压 |
| Fsw | 400 kHz | 开关频率 |
| L | 660 nH | 功率电感 |
| Cout | 250 μF | 输出电容 |
| Rc | 5 mΩ | ESR |
| Rload | 0.1 Ω | 负载电阻 |
| Vslope | 10 mV/周期 | 固定斜坡幅度 |
| Rx | 1 MΩ | SW滤波器电阻 |
| Cx | 200 pF | SW滤波器电容 |

## 输出说明

- **增益曲线**：显示环路增益幅频特性，标注0dB穿越频率（Fc）和开关频率（fs/fs2）
- **相位曲线**：显示环路相位响应，标注相位裕度（PM）和-180°参考线
- **计算结果面板**：显示各模型的Fc和PM数值
- **设计指导面板**：根据当前参数给出优化建议

## 依赖版本

- Python >= 3.7
- numpy >= 1.19
- matplotlib >= 3.3
- control >= 0.9.0
- Pillow >= 8.0
- tkinter（通常随Python提供）

## 参考文献

Lu, D., Zeng, X., & Hong, Z. (2023). *Accurate Loop Gain Model of Ripple-Based Constant On-time Controlled Buck Converters*. IEEE Transactions on Power Electronics, 38(6), 7034-7047.
