# MATLAB 版数字全息重建代码说明

> 本文档对 MATLAB 版文件夹中的全部 12 个 `.m` 文件进行详细说明，涵盖文件角色、调用关系、算法原理和实现细节。

---

## 目录

- [1. 文件清单与角色总览](#1-文件清单与角色总览)
- [2. 整体架构：两条重建链路](#2-整体架构两条重建链路)
- [3. GUI 主入口：VolumeGUI_AngularSpectrum.m](#3-gui-主入口volumegui_angularspectrumm)
- [4. 方法1 核心：Fresnel_reconstructFastStats.m](#4-方法1-核心fresnel_reconstructfaststatsm)
- [5. 方法1 备选：AS_reconstructAngularSpectrumFastStats.m](#5-方法1-备选as_reconstructangularspectrumfaststatsm)
- [6. 方法2 核心：AS_reconstructAngularSpectrum.m](#6-方法2-核心as_reconstructangularspectrumm)
- [7. 辅助模块详解](#7-辅助模块详解)
  - [7.1 AS_preprocessHologram.m — 全息图预处理](#71-as_preprocesshologramm--全息图预处理)
  - [7.2 AS_normalizeImage.m — 图像归一化](#72-as_normalizeimagem--图像归一化)
  - [7.3 AS_buildAdaptiveROI.m — 自适应 ROI 构建](#73-as_buildadaptiveroim--自适应-roi-构建)
  - [7.4 AS_refineMeasurementROI.m — 测量 ROI 细化](#74-as_refinemeasurementroim--测量-roi-细化)
  - [7.5 AS_localizeParticles3D.m — 3D 粒子定位](#75-as_localizeparticles3dm--3d-粒子定位)
  - [7.6 AS_showParticleSpheres.m — 粒子 3D 球体可视化](#76-as_showparticlespheresm--粒子-3d-球体可视化)
  - [7.7 show3d.m / vol3d.m — 3D 体渲染](#77-show3dm--vol3dm--3d-体渲染)
- [8. 关键算法详解](#8-关键算法详解)
- [9. 输出数据结构](#9-输出数据结构)
- [10. 实现细节与注意事项](#10-实现细节与注意事项)

---

## 1. 文件清单与角色总览

| 文件 | 行数 | 角色 | 活跃状态 |
|---|---|---|---|
| `VolumeGUI_AngularSpectrum.m` | 1754 | GUI 主入口，嵌套全部回调函数 | ✅ 活跃 |
| `Fresnel_reconstructFastStats.m` | 550 | **方法1 核心**（菲涅尔快速统计，当前 Python 对齐基线） | ✅ 活跃 |
| `AS_reconstructAngularSpectrumFastStats.m` | 478 | 方法1 备选（角谱版快速统计） | ⚠️ 仅参考 |
| `AS_reconstructAngularSpectrum.m` | 311 | 方法2 核心（角谱完整 3D 重建） | ⚠️ 备用 |
| `AS_preprocessHologram.m` | 39 | 全息图预处理（减背景） | ⚠️ 主路径未调用 |
| `AS_normalizeImage.m` | 49 | 图像线性归一化到 [0,1] | ✅ 活跃 |
| `AS_buildAdaptiveROI.m` | 114 | 为候选粒子构建自适应搜索 ROI | ✅ 活跃 |
| `AS_refineMeasurementROI.m` | 202 | 在焦面 patch 上细化测量窗口 | ✅ 活跃 |
| `AS_localizeParticles3D.m` | 201 | 从 3D 体中定位粒子 | ⚠️ 仅参考 |
| `AS_showParticleSpheres.m` | 63 | 3D 散点/球体粒子可视化 | ✅ 活跃 |
| `show3d.m` | 37 | 3D 体渲染封装 | ✅ 活跃 |
| `vol3d.m` | 196 | 底层 3D 纹理映射渲染引擎 | ✅ 活跃 |

---

## 2. 整体架构：两条重建链路

### 方法1（当前主链路）：菲涅尔快速统计

```
VolumeGUI_AngularSpectrum.m
  └─ onReconstruct()                                    # 重建按钮回调
       ├─ 读取 UI 参数（z 范围、波长、像元尺寸、GPU 等）
       ├─ holoPre = dataImg                             # 直接使用原图（不调预处理）
       ├─ 封装 fastStatsParams 结构体
       └─ Fresnel_reconstructFastStats()                # ★ 核心
            ├─ iFastStatsGpuCore / iFastStatsCpuCore    # GPU/CPU 双路径
            │    ├─ iFrequencyGrid()                    # 频域坐标网格
            │    ├─ phaseBase = -iπλ·(fx²+fy²)          # 菲涅尔传递函数相位
            │    ├─ k_offset = 2π/λ                      # 波数偏移
            │    ├─ fft2(holo)                           # 全息图 FFT
            │    ├─ MIP 循环：逐批反向传播               # H·z → ifft2 → 强度 → MIP
            │    ├─ gather mip2D, argmaxZMap, mean2D     # GPU→CPU 回传
            │    ├─ iDetectCandidates()                  # ★ 候选粒子检测
            │    │    ├─ AS_normalizeImage(mip2D)        # MIP 归一化到 [0,1]
            │    │    ├─ imbinarize + 分水岭分割          # 二值化 + 分离粘连
            │    │    ├─ 启发式筛选：20~35px & 圆度>0.6   # 排除噪声和异常
            │    │    └─ AS_buildAdaptiveROI()           # 自适应搜索 ROI
            │    ├─ iEstimateFocusFromMaps()             # 从 argmaxZMap 估计最佳焦面
            │    ├─ iCollectBestPatchesCpu/Gpu()         # 收集各粒子最佳焦面 patch
            │    └─ iFinalizeSummaryFromMaps()           # ★ 最终统计汇总
            │         ├─ AS_refineMeasurementROI()       # 细化测量窗口
            │         │    └─ 能量质心重定位 + 92%能量半径
            │         └─ iEstimateFocusDiameterFromPatch() # 焦面粒径精测
            └─ 返回 summary 结构体
       └─ onSimpleStatsV2()                             # 统计计算与展示
            ├─ 计算 Sauter 平均径 (D[3,2])
            ├─ 计算体积中位径 (MVD)
            ├─ 计算数浓度
            ├─ 更新直方图、表格、摘要
            └─ iShowParticleOverview()                  # 弹出粒子总览窗口
```

### 方法2（备选链路）：角谱完整 3D 重建

```
AS_reconstructAngularSpectrum()
  ├─ iReconstructGpuCore / iReconstructCpuCore
  │    ├─ iFrequencyGrid()                              # 频域坐标网格
  │    ├─ kz = sqrt(1 - λ²(fx²+fy²))                    # 纵向波数（含倏逝波截止）
  │    ├─ phaseBase = i·k·kz                            # 角谱传递函数相位
  │    ├─ fft2(holo)                                     # 全息图 FFT
  │    ├─ 分批传播：Uz = ifft2(U0f · exp(phaseBase·z))   # 角谱传播
  │    └─ 返回完整 3D 体数据 vol (Ny × Nx × Nz)          # 内存密集型
  └─ AS_localizeParticles3D()                           # 从 3D 体定位粒子
       ├─ MIP 投影 → 二值化 → regionprops
       ├─ AS_buildAdaptiveROI()                          # 自适应 ROI
       ├─ 逐粒子：轴向曲线提取 → 平滑 → 峰值检测
       ├─ AS_refineMeasurementROI()                      # 细化测量窗口
       └─ 返回 coords3D, axialCurves, validMask
```

### 两条链路的本质区别

| 特性 | 方法1（菲涅尔快速统计） | 方法2（角谱完整重建） |
|---|---|---|
| **传播模型** | 菲涅尔傍轴近似 | 角谱法（精确） |
| **传递函数** | `exp(i·k·z - iπλz·(fx²+fy²))` | `exp(i·k·z·sqrt(1-λ²(fx²+fy²)))` |
| **输出** | 统计摘要（无完整 3D 体） | 完整 3D 强度体数据 |
| **内存占用** | 低（逐层 MIP 合并） | 高（保留全 3D 体） |
| **候选筛选** | 分水岭 + 20~35px 粒径 + 圆度>0.6 | 简单阈值 + 无启发式筛选 |
| **GPU 批大小** | 最大 64 层 | 最大 32 层 |
| **每页显存** | ~36 字节/像素 | ~28 字节/像素 |

---

## 3. GUI 主入口：VolumeGUI_AngularSpectrum.m

**文件行数**：1754 行
**角色**：图形用户界面主函数，包含全部嵌套回调函数

### 3.1 界面布局

```
┌──────────────────────────────────────────────────────────────┐
│                    3D Digital Hologram Reconstruction         │
├──────────────┬──────────────────┬────────────────────────────┤
│  左列 (18%)  │   中列 (26%)      │      右列 (49%)             │
│              │                  │                            │
│ 重建参数面板  │  输入图像面板     │  粒子统计面板               │
│ ─────────── │  ─────────────── │  ────────────────────────  │
│ z最小值(mm)  │  [数据图显示]     │  [开始计算] [导出结果]       │
│ z最大值(mm)  │  [加载数据图]     │  [逆衍射逐层测试]            │
│ z步长(µm)   │                  │  状态：空闲                  │
│ 波长(nm)    │  [背景图显示]     │  ┌ 统计摘要 ──────────┐    │
│ 像素尺寸(µm) │  [加载背景图]     │  │ 方法/粒子数/平均径  │    │
│ 重建方法     │                  │  │ MVD/D[3,2]/浓度    │    │
│ [GPU加速]   │                  │  └────────────────────┘    │
│              │                  │  [粒径分布直方图]           │
│              │                  │  [粒子结果表格]             │
│              │                  │  ┌ 预处理与定位参数 ──┐    │
│              │                  │  │ 峰比值/Z边界/2D阈值 │    │
│              │                  │  │ ROI半径倍率/最小半径│    │
│              │                  │  └────────────────────┘    │
└──────────────┴──────────────────┴────────────────────────────┘
```

### 3.2 关键数据结构 `handles`

```matlab
handles = struct(
    DataImage,              % 加载的数据全息图（double 矩阵）
    BgImage,                % 加载的背景图（double 矩阵）
    VolumeData,             % 3D 重建体数据（方法2使用，方法1为空）
    PreprocessedHologram,   % 预处理缓存
    VolumeZVecUm,           % z 轴坐标向量（微米）
    FastStatsSummary,       % 快速统计结果（方法1返回的 summary）
    ParticleCoords3D,       % 粒子 3D 坐标 [N×3]
    ParticleAxialCurves,    % 粒子轴向曲线 cell array
    StatsResult,            % 最终统计结果结构体
    VolumeAlpha,            % 3D 渲染透明度（默认 0.1）
    ZoomApplied,            % 3D 缩放状态标记
    RenderTimer,            % 异步渲染定时器
    RenderRequestId,        % 渲染请求 ID（防抖）
    FastRenderEnabled,      % 快速渲染开关
    FastRenderMaxDim,       % 快速渲染最大维度（256）
    FastRenderTexture       % 纹理模式（'2D'）
);
```

### 3.3 回调函数列表

| 回调 | 触发控件 | 功能 |
|---|---|---|
| `onLoadData` | 加载数据图按钮 | 读取全息图 → 显示 → 清空旧结果 |
| `onLoadBg` | 加载背景图按钮 | 读取背景图 → 显示 → 清空旧结果 |
| `onReconstruct` | 开始计算按钮 | 参数校验 → 调用 Fresnel_reconstructFastStats → onSimpleStatsV2 |
| `onExportResults` | 导出结果按钮 | 导出粒子数据到 Excel/CSV |
| `onPlayInverseTest` | 逆衍射逐层测试按钮 | 逐帧播放角谱传播结果 |
| `onMethodChanged` | 方法下拉框 | 刷新 UI 控件状态 |
| `onSimpleStatsV2` | 菜单项 | 从 FastStatsSummary 构建统计结果并显示 |
| `onAlphaChanged` | 透明度滑块 | 调节 3D 渲染透明度 |
| `onAsyncRender` | 异步定时器 | 防抖后的 3D 渲染执行 |

### 3.4 当前状态说明

- **方法下拉框**仅有一个选项：`'方法1：反向衍射'`
- 主重建路径**不调用** `AS_preprocessHologram`，直接用 `dataImg`
- 导出按钮初始禁用，有统计结果后才启用
- 异步 3D 渲染逻辑存在但当前主流程不走（因为方法1不生成完整 3D 体）

---

## 4. 方法1 核心：Fresnel_reconstructFastStats.m

**文件行数**：550 行
**角色**：当前活跃的方法1基线，Python 重构的对齐目标

### 4.1 函数签名

```matlab
function summary = Fresnel_reconstructFastStats(holoPre, zVec, lambda, pixel, useGPU, params)
```

| 参数 | 类型 | 单位 | 说明 |
|---|---|---|---|
| `holoPre` | single 矩阵 | — | 预处理后全息图 (Ny×Nx) |
| `zVec` | double 向量 | 米 | 重建深度位置 |
| `lambda` | double 标量 | 米 | 激光波长 |
| `pixel` | double 标量 | 米 | 像素间距 |
| `useGPU` | logical | — | 是否尝试 GPU |
| `params` | struct | — | 统计参数（见下表） |

**params 字段**：

| 字段 | 默认值 | 说明 |
|---|---|---|
| `minPeakRatio` | 1.15 | 轴向峰比值阈值 |
| `edgeMargin` | 2 | Z 边界层数 |
| `xyThreshold` | 0.70 | 二维识别阈值 |
| `roiScale` | 0.50 | ROI 半径倍率 |
| `roiMinRadius` | 15 | ROI 最小半径 (px) |
| `measureRoiEnergyFraction` | 0.92 | 测量 ROI 能量占比 |
| `searchRoiGrowFactor` | 1.25 | 搜索 ROI 扩展因子 |
| `edgeEnergyThreshold` | 0.20 | 边缘能量阈值 |
| `searchMinRadiusUm` | 30 | 搜索最小半径 (µm) |

### 4.2 菲涅尔传播数学

**频域网格生成**：
```
fx = ifftshift( -⌊nx/2⌋ : ⌈nx/2⌉-1 ) / (nx · pixel)
fy = ifftshift( -⌊ny/2⌋ : ⌈ny/2⌉-1 ) / (ny · pixel)
```

**传递函数**（菲涅尔近似）：
```
H(z) = exp(i · k · z) · exp(-iπλz · (fx² + fy²))

其中 k = 2π/λ
```

**反向传播**：
```
U(z) = ifft2( fft2(holo) · H(z) )
A(z) = |U(z)|           % 振幅
I(z) = max(A) - A(z)   % 对比度反转：粒子衍射环→亮区
```

**MIP 累积**：
```
对每个像素 (x,y)：
  mip(x,y) = max_z I(x,y,z)
  argmaxZ(x,y) = argmax_z I(x,y,z)
  sum(x,y) = Σ_z I(x,y,z)
```

### 4.3 GPU/CPU 双路径

两个核心函数的差异：

| | iFastStatsGpuCore | iFastStatsCpuCore |
|---|---|---|
| 数据类型 | double（传递函数）+ single（图像） | 统一 single |
| 频域网格 | gpuArray(double) | CPU double |
| MIP 累积 | 在 GPU 上累积（gpuArray） | 在 CPU 上累积 |
| patch 收集 | 逐层 gather 回 CPU | 直接在 CPU 上处理 |
| 批大小 | 按显存的 65% 减去 0.4GB | 固定上限 1.8GB |

---

## 5. 方法1 备选：AS_reconstructAngularSpectrumFastStats.m

**文件行数**：478 行
**角色**：方法1 的角谱版实现，标注为"Reference implementation only"

### 5.1 与 Fresnel 版的关键差异

| 特性 | 菲涅尔版 | 角谱版 |
|---|---|---|
| 传递函数 | `exp(ikz - iπλz(fx²+fy²))` | `exp(ikz · sqrt(1-λ²(fx²+fy²)))` |
| 倏逝波处理 | 无 | 通带掩模 `passband = (1-λ²(fx²+fy²) ≥ 0)` |
| 精度类型 | 混合精度（double 相位 + single 图像） | 纯 single |
| 候选检测 | 分水岭 + 20~35px 筛选 | 简单阈值 + bwareaopen(1) |
| 归一化 | AS_normalizeImage → [0,1] 二值化 | 原始 MIP → levelAbs 二值化 |
| 调试代码 | 无 | **含调试代码**：循环内 figure/imshow/pause/close |

### 5.2 调试残留代码

```matlab
% 第 106-113 行（iFastStatsCore 中）
for k = 1:nz
    H = exp(phaseBase .* zVecW(k)) .* passband;
    I = abs(ifft2(U0f .* H)).^2;
    fig = figure;                                    % ⚠️ 每层弹窗
    for i=1:size(I,3)
        imshow(max(max(I(:,:,i)))-I(:,:,i),[])
        pause(0.01)
    end
    close(fig)
end
```

此代码会在每层重建时弹出一个新 figure 窗口，严重影响性能，是调试遗留。

---

## 6. 方法2 核心：AS_reconstructAngularSpectrum.m

**文件行数**：311 行
**角色**：角谱法完整 3D 重建，标注为"Reference implementation only"

### 6.1 函数签名

```matlab
function vol = AS_reconstructAngularSpectrum(holoPre, zVec, lambda, pixel, useGPU)
```

输出为完整的 `vol` 体数据（Ny × Nx × Nz, single 类型），后续可配合 `AS_localizeParticles3D` 进行粒子定位。

### 6.2 角谱传递函数

```matlab
k = 2π/λ;
kz = sqrt(1 - λ²·(fx² + fy²));      % 纵向波数分量
passband = (1 - λ²·(fx²+fy²) ≥ 0);   % 倏逝波截止
phaseBase = i · k · kz;
H(z) = exp(phaseBase · z) · passband;
```

### 6.3 GPU 显存管理

GPU 路径包含复杂的显存自适应机制：

1. **批大小选择** (`iChooseBatchSizeGPU`)：
   - 可用显存 = `52% × 总显存 - 0.8GB`（安全缓冲）
   - 每页开销 = `nx × ny × (8+8+8+4)` 字节（复数 U0f + 复数 H + 复数 Uz + 实数 I）
   - 大图（≥2048）批大小上限 24，否则 32

2. **体数据存储策略**：
   - 若显存 > `1.4×体积 + 1GB`：体数据全部存在 GPU 上
   - 否则：逐批 gather 回 CPU

3. **OOM 自动回退**：显存不足时批大小减半重试，若仍失败则整体回退 CPU

---

## 7. 辅助模块详解

### 7.1 AS_preprocessHologram.m — 全息图预处理

```matlab
function holoPre = AS_preprocessHologram(dataImg, bgImg)
```

**当前逻辑**：
1. `dataImg` 转为 double
2. 若提供 `bgImg`：尺寸对齐后执行 `holo = holo - bg`（减法，非除法）
3. 直接返回（不做反转）

**注释掉的旧逻辑**（已废弃）：
```matlab
holoPre = max(holo(:)) - holo;   % 全幅反转映射
```

**当前状态**：主重建路径 `onReconstruct` 不调用此函数，直接使用 `dataImg`。

### 7.2 AS_normalizeImage.m — 图像归一化

```matlab
function img = AS_normalizeImage(img)
```

**算法**：线性映射 `img_out = (img - min) / (max - min)`
- 正常情况：输出 [0, 1] 区间 double 型
- 常数图（max == min）：输出全零矩阵

**调用位置**（方法1链路）：
1. `iDetectCandidates`：MIP 图归一化后二值化
2. `iEstimateFocusDiameterFromPatch`：焦面 patch 归一化后粒径估计
3. `AS_localizeParticles3D`：MIP 图归一化后候选检测
4. 逆衍射逐层测试：每帧传播结果归一化显示

### 7.3 AS_buildAdaptiveROI.m — 自适应 ROI 构建

```matlab
function roi = AS_buildAdaptiveROI(img2D, candidate, params)
```

**算法步骤**：

1. **初始半径计算**：
   ```
   radius = max(roiMinRadius, bboxRadius, equivDiam × roiScale, searchMinRadiusUm/pix_um)
   radius = min(radius, roiMaxRadiusPx)
   ```

2. **边缘能量驱动扩展**：
   ```
   while edgeEnergy > threshold && radius < maxRadius:
       radius = ceil(radius × searchRoiGrowFactor)   # 1.25 倍增长
       重新裁剪 → 计算新 edgeEnergy
   ```

3. **边缘能量定义**：
   ```
   border = 四边像素值
   bg = median(patch)
   peak = max(patch)
   edgeEnergy = (max(border) - bg) / max(peak - bg, eps)
   ```
   物理含义：边缘像素有多接近内部峰值。高值表示粒子能量溢出当前 ROI。

### 7.4 AS_refineMeasurementROI.m — 测量 ROI 细化

```matlab
function out = AS_refineMeasurementROI(focusSlice, centerXY, searchBox, params)
```

**核心算法**：

1. **裁剪搜索区域**：从 `focusSlice` 中裁出 `searchBox` 指定区域

2. **背景估计**：
   ```
   border = [第一行, 最后一行, 第一列(去头尾), 最后一列(去头尾)]
   bgLevel = median(border)
   ```

3. **能量图重定位**：
   ```
   energyMap = max(patch - bgLevel, 0)
   refinedCenter = 能量图质心坐标
   ```

4. **能量分数半径**：
   ```
   按像素到质心的距离升序排序
   累计能量达到 targetFraction（默认 92%）时的距离 = 测量半径
   ```

5. **输出**：包含细化后的局部/全局坐标、测量半径、背景电平

**与搜索 ROI 的区别**：
| | 搜索 ROI | 测量 ROI |
|---|---|---|
| 目的 | 定位粒子焦面 | 精确测量粒径 |
| 大小策略 | 从候选中心向外扩展 | 从能量质心向内收缩 |
| 最小半径 | 15 px | 3 px |
| 扩展因子 | 1.25× | 不扩展，收束到 92% 能量 |

### 7.5 AS_localizeParticles3D.m — 3D 粒子定位

```matlab
function [coords3D, axialCurves, roiBoxes, bw, img2D, stats, validMask] = ...
    AS_localizeParticles3D(volData, zVecUm, minPeakRatio, edgeMargin, ...
                           xyThreshold, roiScale, roiMinRadius, roiOptions)
```

**算法步骤**（对应方法2的后处理）：

1. MIP 投影：`mip2D = max(volData, [], 3)`
2. 二值化：`bw = mip2D > xyThreshold × max(mip2D(:))`
3. 形态学：`bwareaopen(bw, 1)` + `imfill(bw, 'holes')`
4. `regionprops` 提取候选
5. 逐候选粒子：
   - `AS_buildAdaptiveROI` 构建搜索窗口
   - 沿 z 轴提取轴向曲线（ROI 内均值）
   - `smoothdata` 平滑后找峰值
   - `AS_refineMeasurementROI` 细化测量窗口
   - `iEstimateFocusDiameter` 焦面粒径精测
6. 有效性判定：`valid = 远离Z边界 AND peakRatio ≥ minPeakRatio`

### 7.6 AS_showParticleSpheres.m — 粒子 3D 球体可视化

```matlab
function hParticleFig = AS_showParticleSpheres(coords3D, diams_um, pix_um, parentAx)
```

- 将 x/y 从像素转换为微米
- 每个粒子渲染为一个彩色半透明球体（`FaceAlpha = 0.95`）
- 球体半径 = `max(baseRadius, 实际直径 × 1.2 / 2)`
- 颜色按 `lines` colormap 循环分配
- 使用 `camlight headlight` + `gouraud` 光照

### 7.7 show3d.m / vol3d.m — 3D 体渲染

**show3d.m**（37行）：
- 封装 `vol3d` 调用
- 设置坐标轴刻度、视角 `(-202°, 18°)`
- 反转 hot colormap（高值暗、低值亮，突出粒子）
- 调节 `alphamap('decrease', alpham)` 透明度

**vol3d.m**（196行，Joe Conti 2004）：
- 底层 3D 纹理映射引擎
- 支持 `'2D'`（仅最接近相机方向）和 `'3D'`（全部三个正交方向）纹理模式
- 通过堆叠半透明表面切片模拟体渲染
- 控制 `drawnow` 频率（每 `max(siz)/12` 层刷新一次）以平衡性能

---

## 8. 关键算法详解

### 8.1 菲涅尔 vs 角谱传递函数

```
菲涅尔（傍轴近似）：
  H(z) = exp(i·k·z) · exp(-iπλz·(fx²+fy²))

  条件：观测距离远大于波长，旁轴光线
  优点：计算量小，相位项是纯二次型
  缺点：大角度衍射精度下降

角谱（严格解）：
  H(z) = exp(i·k·z·sqrt(1 - λ²(fx²+fy²)))
  条件：1 - λ²(fx²+fy²) ≥ 0（传播模态）
  优点：对所有角度精确，自动滤除倏逝波
  缺点：sqrt 计算开销更大
```

### 8.2 候选粒子检测（两个版本的差异）

| 步骤 | 菲涅尔版（活跃） | 角谱版（参考） |
|---|---|---|
| 归一化 | `AS_normalizeImage(mip2D)` → [0,1] | 直接用原始 MIP |
| 二值化 | `imbinarize(img, level)` | `mip2D > levelRatio × max(mip2D(:))` |
| 分割 | 分水岭算法分离粘连粒子 | 简单阈值（不分离粘连） |
| 去噪 | `bwareaopen(bw, 3)` | `bwareaopen(bw, 1)` |
| 筛选 | **直径 20~35px + 圆度 > 0.6** | **无** |
| 中心 | WeightedCentroid（灰度加权） | WeightedCentroid，fallback Centroid |

> ⚠️ **关键差异**：方法1的"20~35px 粒径 + 圆度>0.6"启发式筛选是同一阈值在两个方法上效果不同的根本原因。这也是 Python 重构版需要特别注意的地方。

### 8.3 焦面粒径精测

```
iEstimateFocusDiameterFromPatch(roiFocus, cx, cy, xyThreshold, fallbackDiam):
  1. AS_normalizeImage(roiFocus)                 → 归一化到 [0,1]
  2. graythresh(roiFocus) 或 xyThreshold        → 自动/手动阈值
  3. level = clamp(Otsu, 0.35, 0.85)            → 限制阈值范围
  4. bwFocus = imbinarize → bwareaopen → imfill  → 二值化+清理
  5. bwconncomp → regionprops                    → 连通域分析
  6. 找到包含原候选中心的连通域                   → 选择正确粒子
  7. 返回 EquivDiameter                          → 等效直径
```

### 8.4 体积中位径 (MVD) 计算

```matlab
% 按直径排序
[diamSorted, sortIdx] = sort(diamVecUm);

% 每个粒子的体积权重 ∝ d³
volWeights = diamSorted .^ 3;
cumFrac = cumsum(volWeights) / sum(volWeights);

% 找到累计体积占比达到 50% 的直径（线性插值）
idx50 = find(cumFrac >= 0.5, 1, 'first');
mvdUm = interp(cumFrac(idx50-1:idx50), diamSorted(idx50-1:idx50), 0.5);
```

### 8.5 数浓度计算

```matlab
zSpanUm = |zMax - zMin| + zStep;                % z 方向有效跨度
reconVolumeUm3 = rows × cols × pix_um² × zSpanUm; % 重建体积 (µm³)
numConcPerMl = particleCount × 1e12 / reconVolumeUm3; % 个/mL
```

---

## 9. 输出数据结构

### 9.1 summary（方法1 返回）

```matlab
summary = struct(
    img2D,              % 归一化后的 MIP 图（double, [0,1]）
    mip2D,              % 原始 MIP 图（single）
    bw,                 % 二值化掩模
    zVecUm,             % z 轴坐标 [nz×1] (µm)
    volumeSize,         % [ny, nx, nz]
    hasVolume,          % false（方法1不保留完整3D体）
    validMask,          % 有效性掩码 [nCand×1] logical
    candidateCoords3D,  % 所有候选坐标 [nCand×3] (x_px, y_px, z_um)
    candidateRoiBoxes,  % 候选 ROI 框 [nCand×4]
    candidateStats,     % 候选 regionprops 统计
    coords3D,           % 有效粒子坐标 [nValid×3]
    axialCurves,        % 有效粒子轴向曲线 cell array
    roiBoxes,           % 有效粒子 ROI 框 [nValid×4]
    stats,              % 有效粒子统计
    diamsPx,            % 有效粒子粒径 [nValid×1] (px)
    mode,               % 运行模式：'GPU'/'CPU'
    totalSec            % 总耗时 (秒)
);
```

### 9.2 StatsResult（GUI 统计展示用）

```matlab
statsResult = struct(
    methodName,         % '方法1：反向衍射'
    coordsPx,           % 粒子坐标(像素) [N×3]
    coordsUm,           % 粒子坐标(微米) [N×3]
    diamsUm,            % 粒径(微米) [N×1]
    diamsPx,            % 粒径(像素) [N×1]
    pixUm,              % 像元尺寸(µm)
    zVecUm,             % z 坐标向量
    zStepUm,            % z 步长(µm)
    zRangeUm,           % z 范围 [min, max]
    count,              % 有效粒子数
    meanDiamUm,         % 算术平均粒径
    mvdUm,              % 体积中位径 (MVD)
    meanEffDiamUm,      % Sauter 平均径 D[3,2]
    numConcPerMl,       % 数浓度 (个/mL)
    mipImage,           % MIP 图（用于可视化）
    candidateCoordsPx,  % 候选坐标
    candidateRoiBoxes,  % 候选 ROI 框
    candidateValidMask, % 候选有效性
    reconRows,          % 重建行数
    reconCols,          % 重建列数
    paramSnapshot       % 参数快照
);
```

---

## 10. 实现细节与注意事项

### 10.1 活跃基线说明

- **`Fresnel_reconstructFastStats.m`** 是当前唯一活跃的方法1基线，文件头标注 `"Active method-1 baseline kept for Python alignment"`
- `AS_reconstructAngularSpectrumFastStats.m` 和 `AS_localizeParticles3D.m` 均标注为 `"Reference implementation only"`
- Python 重构版应以 Fresnel 版本为目标对齐

### 10.2 预处理状态

- `AS_preprocessHologram` 的旧版逻辑（`max(image)-image` 全幅反转）已被注释掉
- 新版简化为纯减背景（不反转）
- 但在 `onReconstruct` 主路径中，连新版也不调用，直接用 `dataImg`
- 只有在"逆衍射逐层测试"路径中，若没有缓存预处理结果，才会调用 `AS_preprocessHologram`

### 10.3 GPU 回退机制

所有 GPU 路径均有完善的自动回退：

```
try GPU → 失败 → warning + fallback CPU
  └─ OOM → 批大小减半 → 重试
       └─ 仍失败 → 整体回退 CPU
```

### 10.4 调试残留

`AS_reconstructAngularSpectrumFastStats.m` 第 106-113 行存在调试代码：
```matlab
fig = figure;
for i=1:size(I,3)
    imshow(max(max(I(:,:,i)))-I(:,:,i),[])
    pause(0.01)
end
close(fig)
```
**若要在生产环境使用此文件，需删除此段代码**。

### 10.5 参数默认值差异

不同函数中的 fallback 默认值不完全一致：

| 参数 | GUI 默认 | Fresnel 默认 | AS FastStats 默认 |
|---|---|---|---|
| `xyThreshold` | 0.65 (控件) / 0.85 (fallback) | 0.70 | 0.85 |
| `minPeakRatio` | 1.15 | 1.15 | 1.15 |
| `edgeMargin` | 2 | 2 | 2 |
| `roiScale` | 0.50 | 0.50 | 0.50 |
| `roiMinRadius` | 15 | 15 | 15 |

### 10.6 vol3d.m 版权

`vol3d.m` 为 Joe Conti 2004 年的第三方代码，通过正交平面 2D 纹理映射实现 3D 体渲染，是 show3d 的底层依赖。

---

*文档生成时间：2026-04-28*
*基于 MATLAB 版文件夹全部 12 个 .m 文件的完整阅读分析*
