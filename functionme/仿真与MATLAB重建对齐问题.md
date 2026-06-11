# Python 仿真全息图与 MATLAB 重建管线对齐问题

> 本文档记录 Python 端 `generate_hologram_matlab_validation_new.py` 生成的合成全息图与 MATLAB 端 `Fresnel_reconstructFastStats.m` 重建管线之间**除 UI 可调参数外**的代码级不对齐问题。
>
> 生成时间：2026-06-11　|　最后更新：2026-06-11（P1/P2/P3 已通过「仿真模式」开关解决）

---

## 0. 修复状态总览

| 问题 | 状态 | 解决方案 |
|---|---|---|
| P1 硬编码粒径筛选 | ✅ 已修复 | `simulationMode=true` 时自动放宽 `minDiamPx/maxDiamPx/minCircularity` |
| P2 8-bit 编码断裂 | ✅ 已修复 | `simulationMode=true` 时 `AS_preprocessHologram` 做 `/255 - DC` |
| P3 对比度反转 | ✅ 已修复 | `simulationMode=true` 时默认 `invertContrast=false` |
| P4 传播模型差异 | ⚠️ 未修复 | 仅文档记录，实际影响较小 |
| P5 Born 近似 vs 实验特征 | ⚠️ 未修复 | 仿真模式下分水岭 H-minima 已降至 0.3 |
| P6 坐标系统 | ⚠️ 无需修复 | 对比真值时注意 +1 即可 |

**使用方式**：在 MATLAB GUI 左栏勾选「Simulation mode」复选框，即可自动应用所有仿真适配。不勾选时完全保持原有实验流程不变。

---

## 目录

- [1. 问题总览](#1-问题总览)
- [2. 致命问题：硬编码粒径筛选器](#2-致命问题硬编码粒径筛选器)
- [3. 严重问题：8-bit 编码/解码链路断裂](#3-严重问题8-bit-编码解码链路断裂)
- [4. 严重问题：MIP 对比度反转假设](#4-严重问题mip-对比度反转假设)
- [5. 中等问题：传播模型差异](#5-中等问题传播模型差异)
- [6. 中等问题：全息图形成模型的物理假设差异](#6-中等问题全息图形成模型的物理假设差异)
- [7. 轻微问题：坐标系统与索引约定](#7-轻微问题坐标系统与索引约定)
- [8. 影响链路总结](#8-影响链路总结)
- [9. 修复优先级建议](#9-修复优先级建议)

---

## 1. 问题总览

| 编号 | 问题 | 严重程度 | 影响 | 是否需要改代码 |
|---|---|---|---|---|
| P1 | 硬编码 20-35 px 粒径筛选 | 🔴 致命 | 仿真粒子全部被丢弃 | ✅ MATLAB |
| P2 | 8-bit BMP 编码/解码链路断裂 | 🔴 严重 | 全息图物理强度信息丢失，DC 分量异常 | ✅ 双端 |
| P3 | `I = max(Amp) - Amp` 对比度反转假设 | 🔴 严重 | MIP 中粒子信号弱/不可见 | ✅ MATLAB |
| P4 | 角谱 vs 菲涅尔传播模型 | 🟡 中等 | 大角度/高 NA 时重建精度下降 | ⚠️ 视情况 |
| P5 | Born 近似 vs 实验全息图特征 | 🟡 中等 | 启发式阈值行为不一致 | ⚠️ MATLAB |
| P6 | 坐标系统（0-indexed vs 1-indexed） | 🟢 轻微 | 真值对比时偏移 | ⚠️ 后处理 |

---

## 2. 致命问题：硬编码粒径筛选器

### 位置

`MATLAB版/Fresnel_reconstructFastStats.m` — 内部函数 `iDetectCandidates`，第 384 行：

```matlab
keep_idx = (all_diams >= 20) & (all_diams <= 35) & (circularity > 0.6);
```

### 问题描述

该行在 MIP 投影 → 归一化 → 二值化 → 分水岭分割 → `regionprops` 之后，用 **像素直径 20~35 px** 和 **圆度 > 0.6** 作为硬过滤条件。

`iDetectCandidates` 的完整处理链：

```matlab
% 第 353-362 行：二值化 + 分水岭
img2D = AS_normalizeImage(mip2D);
bw = imfill(bwareaopen(imbinarize(img2D, level), 3), 'holes');
D = -bwdist(~bw);
mask = imextendedmin(D, 0.8);
D2 = imimposemin(D, mask);
L_watershed = watershed(D2);
bw(L_watershed == 0) = 0;

% 第 364 行：提取候选
stats_raw = regionprops(bw, img2D, 'Centroid', 'WeightedCentroid', ...
    'EquivDiameter', 'BoundingBox', 'Area', 'Perimeter');

% 第 380-384 行：硬过滤 ← 问题所在
all_diams = [stats_raw.EquivDiameter];
keep_idx = (all_diams >= 20) & (all_diams <= 35) & (circularity > 0.6);
stats = stats_raw(keep_idx);
```

### 数据验证

以已生成的仿真数据集为例：

**数据集 `0001`（Python 默认参数：d=10-20 µm, pixel=2.2 µm）**：

| 粒子 | 物理直径 (µm) | 像素直径 (px) | 20-35 px 筛选 |
|---|---|---|---|
| 1 | 16.26 | **7.4** | ❌ |
| 2 | 17.08 | **7.8** | ❌ |
| 5 | 17.19 | **7.8** | ❌ |
| 12 | 19.95 | **9.1** | ❌ |

**数据集 `0011`（d=100 µm, pixel=2.2 µm）**：

| 粒子 | 物理直径 (µm) | 像素直径 (px) | 20-35 px 筛选 |
|---|---|---|---|
| 全部 15 个 | 100.0 | **≈45.5** | ❌ |

> ⚠️ 这两个典型参数组合下，**100% 的仿真粒子都会被过滤掉**，导致重建结果为空。

### 根本原因

`[20, 35]` 这个区间是为特定实验条件硬编码的：
- 预期像元尺寸 3.45 µm，对应物理直径 **69–121 µm**
- 这个条件写在 `iDetectCandidates` 内部，**无法通过 UI 参数调整**
- 与 `params.roiMinRadius`、`params.roiScale` 等可配置参数不同，此筛选器完全不受外部控制

### 为什么实验全息图需要这个筛选器

真实实验全息图的 MIP 中，噪声、离焦粒子拖影、光学伪影会产生大量小尺寸假阳性候选。20-35 px 筛选器配合圆度 > 0.6 是为了排除这些。但仿真全息图的噪声特征完全不同，不能直接复用同一套阈值。

---

## 3. 严重问题：8-bit 编码/解码链路断裂

### 3.1 Python 端：`_to_uint8` 的非线性编码

**位置**：`仿真图像生成/generate_hologram_matlab_validation_new.py`，第 9-22 行：

```python
def _to_uint8(image):
    clipped = image.copy()
    low = np.percentile(clipped, 0.1)      # ← 裁剪下界
    high = np.percentile(clipped, 99.9)    # ← 裁剪上界
    clipped[clipped < low] = low
    clipped[clipped > high] = high
    image_min = clipped.min()
    image_max = clipped.max()
    normalized = (clipped - image_min) / (image_max - image_min)
    return np.uint8(normalized * 255)
```

**问题**：

1. **百分位截断丢失极值**：0.1% 和 99.9% 的百分位裁剪会抹掉全息图干涉条纹中最亮和最暗的像素，这些恰好是携带粒子高频信息的关键区域
2. **线性映射到 [0,255] 后丢失绝对物理标度**：原始全息图的物理强度值域（如 `|1 - U_total|²` 可能分布在 [0.5, 1.5] 区间）被压缩到 [0, 255]，无法恢复
3. **DC 分量被固定到 128 附近**：因为映射是全局线性的，全息图的平均强度被映射到 128 左右，而这个值对 MATLAB 端有显著影响

### 3.2 MATLAB 端：uint8 像素值直接当复数场

**位置**：`MATLAB版/VolumeGUI_AngularSpectrum.m`，第 630 行：

```matlab
holoPre = dataImg;   % dataImg 来自 imread()，类型 uint8，值域 [0, 255]
```

然后传入 `Fresnel_reconstructFastStats`：

```matlab
% Fresnel_reconstructFastStats.m，第 25 行
holo = single(holoPre);   % → 值域仍是 [0, 255]，而非物理强度
```

随后直接做 FFT：

```matlab
% 第 78 行
U0f = fft2(holo);   % 对一个均值 ≈128 的矩阵做 FFT
```

**物理后果**：

- 全息图 `I(x,y) = |R + O|²` 的真实 DC 分量是 `|R|² + ⟨|O|²⟩` ≈ 1 + ε，很小的偏移
- 但 `_to_uint8` 映射后 DC ≈ 128（uint8 的中值），经过 FFT 后在零频产生一个 **巨大的尖峰**
- 这个 DC 尖峰传播到所有 z 层，在重建体中产生强烈的背景偏移
- 粒子信号（由干涉条纹的边带携带）相对于 DC 的幅度比被严重压缩

### 3.3 Python 开发者的原始意图（未被执行）

**位置**：`generate_hologram_matlab_validation_new.py`，第 150-153 行：

```python
print(" 1. I = imread('hologram.bmp');")
print(" 2. I = double(I) / 255.0;")
print(" 3. I = I - mean(I(:));  % 在 Matlab 内减去直流")
```

这些建议是**正确的**：除以 255 恢复到 [0,1]，减均值去除 DC。但是：

1. 这只是 `print` 输出，没有被 MATLAB 代码执行
2. `_to_uint8` 已经不可逆地改变了数据的分布（百分位截断 + 线性映射），即使执行这三步也无法恢复原始物理值
3. MATLAB GUI 的 `onReconstruct` 完全没有这些步骤

### 3.4 正确的数据流应该是什么

```
方案 A（推荐：绕过 BMP，保存浮点数据）：
  Python: np.save('hologram.npy', clean_hologram.astype(np.float32))
  MATLAB: holo = single(readNPY('hologram.npy'));  % 直接获得物理强度值

方案 B（修复 BMP 管线，精度有限）：
  Python: 不做 _to_uint8，将物理值线性映射到 uint16（65536 级）
          或先做 log 压缩再映射以保留动态范围
  MATLAB: imread → double → 逆向映射恢复近似物理值 → 减 DC
```

---

## 4. 严重问题：MIP 对比度反转假设

### 位置

`MATLAB版/Fresnel_reconstructFastStats.m`，GPU 路径第 89-96 行（CPU 路径同理）：

```matlab
H = exp(1i * k_offset * zBatch + phaseBase .* zBatch);
Amp = abs(ifft2(U0f .* single(H)));
sliceMax = max(max(Amp, [], 1), [], 2);
I = bsxfun(@minus, sliceMax, Amp);    % ← 对比度反转
```

### 物理含义

`I = max(Amp) - Amp` 将振幅图反转：原本振幅**低**的像素（粒子阴影区）变成**高**值，在 MIP 累积 (`max over z`) 中成为亮区。

### 对实验全息图成立的假设

实验全息图中，粒子在重建振幅中表现为**暗区**（因为粒子遮挡/散射了照明光），反转后粒子变亮。同时背景（无粒子区域）的振幅波动是均匀的，反转后变成接近 0 的值。

### 对仿真全息图的失效模式

仿真全息图的粒子在重建振幅中的表现取决于：

1. **全息图的整体 DC 偏移**：经过 `_to_uint8` 后 DC ≈ 128，传播后会产生强的背景振幅
2. **粒子是"完美"不透明圆盘**：Born 近似下的散射场没有实验粒子那样的复杂衍射结构
3. **没有实验噪声基底**：反转后 `I = max - Amp`，如果整个切片的 Amp 都接近均匀（仿真全息图的重建），则 `I` 中几乎没有粒子信号

**可能的后果**：

- 粒子在 MIP 中信号很弱，被背景波动淹没
- `iDetectCandidates` 中 `imbinarize(img2D, level)` 可能完全找不到粒子区域（因为对比度太低）
- 即使粒子被检测到，`AS_refineMeasurementROI` 的"能量质心"算法会因为弱信号而定位不准

### 需要注意的关联影响

`iDetectCandidates` 第 353 行的归一化：

```matlab
img2D = AS_normalizeImage(mip2D);   % → [0, 1]
level = min(max(params.xyThreshold, 0), 1);
bw = imfill(bwareaopen(imbinarize(img2D, level), 3), 'holes');
```

如果 MIP 中粒子信号本来就弱（因为上面的反转假设不成立），那么即使降低 `xyThreshold`，二值化结果也会很不可靠——要么阈值太高漏掉粒子，要么阈值太低引入大量噪声。

---

## 5. 中等问题：传播模型差异

### 5.1 两端的传递函数

**Python 生成（角谱法，严格解）**：

```python
# generate_hologram_matlab_validation_new.py，第 109-110 行
sqrt_arg = 1 - (wavelength * fx_grid) ** 2 - (wavelength * fy_grid) ** 2
valid_mask = sqrt_arg >= 0
transfer[valid_mask] = np.exp(1j * k * zi * np.sqrt(sqrt_arg[valid_mask]))
```

传递函数：$$H_{AS}(z) = \exp\left(j \cdot k \cdot z \cdot \sqrt{1 - \lambda^2(f_x^2 + f_y^2)}\right)$$

自动滤除倏逝波（$$1 - \lambda^2(f_x^2 + f_y^2) < 0$$ 的区域直接置零）。

**MATLAB 重建（菲涅尔傍轴近似）**：

```matlab
% Fresnel_reconstructFastStats.m，第 73-74 行
phaseBase = (-1i * pi * lambda) .* (FX.^2 + FY.^2);
k_offset = (2 * pi / lambda);
H = exp(1i * k_offset * zBatch + phaseBase .* zBatch);
```

传递函数：$$H_F(z) = \exp(jkz) \cdot \exp(-j\pi\lambda z(f_x^2 + f_y^2))$$

### 5.2 差异分析

菲涅尔近似来自角谱的一阶泰勒展开：

$$\sqrt{1 - \lambda^2(f_x^2+f_y^2)} \approx 1 - \frac{1}{2}\lambda^2(f_x^2+f_y^2)$$

代入得：$$\exp(jkz \cdot 1 - jkz \cdot \frac{1}{2}\lambda^2(f_x^2+f_y^2)) = \exp(jkz) \cdot \exp(-j\pi\lambda z(f_x^2+f_y^2))$$

截断误差在 $$\lambda^4(f_x^2+f_y^2)^2$$ 量级。

### 5.3 实际影响评估

| 条件 | 最大 $$f$$ | $$\lambda^4 f^4$$ 项 | 影响 |
|---|---|---|---|
| N=256, dx=2.2µm | $$f_{max} \approx$$ 0.23 µm⁻¹ | ~0.0006 | ⚠️ 可察觉 |
| N=2048, dx=2.2µm | $$f_{max} \approx$$ 0.23 µm⁻¹ | ~0.0006 | ⚠️ 可察觉 |
| N=256, dx=3.45µm | $$f_{max} \approx$$ 0.14 µm⁻¹ | ~0.0001 | ✅ 可忽略 |

对于 Python 默认的 dx=2.2 µm，Nyquist 频率处菲涅尔近似有轻微偏差。对于粒子边缘的精细衍射环（高空间频率成分），重建位置可能有亚像素级的偏移。

**结论**：传播模型差异**不是导致"完全检测不到粒子"的原因**，但会对重建精度（尤其是 z 轴定位和粒径测量）产生亚像素到 1-2 像素级别的影响。

---

## 6. 中等问题：全息图形成模型的物理假设差异

### 6.1 Python 端的物理模型（Born 近似 + 二元不透明圆盘）

```python
# 每个粒子：二元掩模（0 = 透明, 1 = 不透明）
obj_plane = np.zeros((N, N), dtype=float)
obj_plane[dist_sq <= r_pixel**2] = 1.0

# 角谱传播到探测器平面
obj = fftshift(fft2(ifftshift(obj_plane)))
propagated = obj * exp(1j * k * zi * sqrt(...))
u_plane = fftshift(ifft2(ifftshift(propagated)))
u_total += u_plane

# 全息图 = |参考光 - 总散射场|²
hologram = |1 - u_total|²
```

**关键简化**：

1. **Born 近似（单次散射）**：每个粒子的散射场独立计算后线性叠加，忽略粒子间的多次散射
2. **二元不透明圆盘**：粒子内部完全挡光（振幅 = 0），粒子外部完全透明（振幅 = 1），没有部分透明、折射、相位延迟
3. **平面波照明**：参考光是理想的单位振幅平面波

### 6.2 MATLAB 端预期的实验全息图特征

MATLAB 管线（尤其是 `iDetectCandidates` 中的一系列启发式阈值）是在**实验全息图**上调优的，实验全息图通常有以下仿真缺少的特征：

| 特征 | 实验全息图 | 仿真全息图 |
|---|---|---|
| 粒子边缘 | 渐变（折射+衍射） | 锐利（二元掩模） |
| 衍射环渐晕 | 自然衰减（有限相干+像差） | 理想 Airy 图样（无穷级次） |
| 背景非均匀性 | 照明不均匀、探测器渐晕 | 完全均匀（|1|² = 1） |
| 噪声 | 散斑、暗电流、读出噪声 | 无（或加性高斯白噪声） |
| 粒子间差异 | 形状、折射率、取向 | 完全相同（仅大小不同） |
| MIP 中粒子外观 | 模糊/扩散的亮斑 | 紧凑的亮点 |

### 6.3 受影响的具体代码位置

**分水岭分割的深度参数**（`iDetectCandidates` 第 359 行）：

```matlab
mask = imextendedmin(D, 0.8);   % H-minima 深度 = 0.8
```

`0.8` 这个值是在归一化后的实验 MIP 上调出来的。仿真 MIP 的对比度特征完全不同，同样的 `0.8` 可能导致过分割或欠分割。

**边缘能量阈值**（`AS_buildAdaptiveROI` 默认 0.20）：

```matlab
edgeEnergy = (max(border) - bg) / max(peak - bg, eps);
if edgeEnergy <= 0.20 || radius >= maxRadius
    break;   % ROI 扩张停止条件
end
```

仿真粒子在焦面 patch 中的边缘更锐利，`edgeEnergy` 的计算值会与实验粒子不同，导致 ROI 扩张行为不一致。

---

## 7. 轻微问题：坐标系统与索引约定

### 7.1 0-indexed vs 1-indexed

**Python**（0-indexed）：

```python
xo = np.random.uniform(margin, N - margin)  # 范围 [margin, N-margin)
yo = np.random.uniform(margin, N - margin)
# CSV 输出为浮点像素坐标，原点在 (0, 0)
```

**MATLAB**（1-indexed）：

```matlab
% imread 读入的图像，左上角像素为 (1,1)
% WeightedCentroid 返回 1-indexed 坐标
% 真值 CSV 需要 x+1, y+1 后才能与 MATLAB 坐标对齐
```

Python 代码打印的提示已指出此问题：

```python
print(" 4. 真值 CSV 中的 (x, y) 坐标需要 +1 才是 Matlab 中的真实像素坐标。")
```

### 7.2 亚像素精度的注意事项

Python 的 `xo`、`yo` 是连续的浮点值（粒子的精确中心），CSV 中也保留了浮点精度。但 MATLAB 端：

- `iDetectCandidates` 使用 `WeightedCentroid`（亚像素精度）
- `AS_refineMeasurementROI` 用能量质心重定位（亚像素精度）

两者都有亚像素精度，对比时不应简单取整。

### 7.3 X/Y 轴方向

Python 代码：

```python
y_grid, x_grid = np.ogrid[0:N, 0:N]    # y 是第 0 维（行），x 是第 1 维（列）
dist_sq = (x_grid - xo) ** 2 + (y_grid - yo) ** 2
```

MATLAB 代码中 `regionprops` 的 `Centroid` 返回 `[column, row]` = `[x, y]`。Python CSV 中的 `x_pixel` 对应列（column），`y_pixel` 对应行（row）。**方向是一致的**，只需注意 +1 偏移。

---

## 8. 影响链路总结

```
Python 生成全息图
│
├─[P5] Born 近似 + 二元圆盘              → "过于完美"的衍射图样
├─[P4] 角谱法传播                        → 与菲涅尔重建有模型偏差
└─[P2] _to_uint8: 百分位截断+线性映射    → 物理强度信息不可逆丢失
        │
        ▼ 8-bit BMP 文件
        │
MATLAB 读取
│
├─[P2] uint8[0,255] 直接当复数场用       → 巨大 DC 分量
├─[P6] 无 1→0 索引转换                  → 坐标偏移（已记录）
│
├─ fft2(holo) → ×H(z) → ifft2 → abs      [菲涅尔反向传播]
│   │
│   ├─[P4] 与原始角谱传播的模型偏差
│   └─[P3] Amp 对比度反转 I=max-Amp      → 粒子信号可能很弱
│           │
│           ▼ MIP 累积
│           │
├─ iDetectCandidates
│   ├─[P1] 🔴 (diams>=20)&(diams<=35)    → 全部仿真粒子被丢弃！
│   ├─[P6] imextendedmin(D, 0.8)         → 对仿真 MIP 不适用
│   └─[P2] 归一化后的对比度不够          → 二值化失败
│
├─ AS_buildAdaptiveROI
│   └─[P6] edgeEnergy 阈值 0.20         → 仿真粒子的边缘能量不同
│
└─ AS_refineMeasurementROI
    └─[P3] 能量图质心定位                → 弱信号下精度下降
```

---

## 9. 修复优先级建议

### 第一优先级（必须修复，否则完全无法工作）

| 序号 | 修改 | 文件 | 方案 |
|---|---|---|---|
| 1 | **移除/参数化粒径硬编码筛选** | `Fresnel_reconstructFastStats.m` → `iDetectCandidates` | 将 `[20, 35]` 改为可配置参数 `params.minDiamPx` / `params.maxDiamPx`，或根据 `pixel` 和物理粒径范围自动计算 |
| 2 | **MATLAB 端添加归一化与 DC 去除** | `VolumeGUI_AngularSpectrum.m` → `onReconstruct` | 在 `holoPre = dataImg` 后加 `holoPre = double(holoPre) / 255; holoPre = holoPre - mean(holoPre(:));` |
| 3 | **Python 端保存浮点数据** | `generate_hologram_matlab_validation_new.py` | 除 BMP 外，同步输出 `.mat` 或 `.npy` 格式的浮点全息图，供 MATLAB 精确读取 |

### 第二优先级（改善检测率和精度）

| 序号 | 修改 | 文件 | 方案 |
|---|---|---|---|
| 4 | **验证并调整对比度反转** | `Fresnel_reconstructFastStats.m` | 在仿真数据上测试 `I = Amp`（不反转）vs `I = max-Amp` 的效果，必要时添加 `params.invertContrast` 开关 |
| 5 | **降低 `imextendedmin` H-minima 值** | `Fresnel_reconstructFastStats.m` | 仿真全息图的 MIP 噪声更少，可以用更低的 H-minima（如 0.3-0.5）来避免欠分割 |
| 6 | **对齐传播模型** | `Fresnel_reconstructFastStats.m` 或 Python 端 | 选择一端统一：要么 Python 改用菲涅尔生成，要么 MATLAB 改用角谱重建 |

### 第三优先级（锦上添花）

| 序号 | 修改 | 文件 | 方案 |
|---|---|---|---|
| 7 | **添加自动化端到端测试** | 新建 | 固定随机种子生成仿真数据 → MATLAB 重建 → 自动对比真值，量化检测率/定位误差/粒径误差 |
| 8 | **ROI 边缘能量阈值参数化** | `AS_buildAdaptiveROI.m` | `edgeEnergyThreshold` 已可配但默认 0.20 对仿真偏严格，可通过 UI 暴露 |

---

> ⚠️ **注意**：本文档**不包含**可通过 MATLAB GUI 直接调整的参数差异（波长、像元尺寸、Z 范围、Z 步长、xyThreshold、roiScale 等）。这些参数需在使用时根据仿真设定手动匹配。
