---
name: hologram-pipeline
description: 数字全息仿真重建管线专家。掌握 Python→MATLAB 全链路、对齐问题诊断、代码修改与调试。处理全息图生成、重建、粒子检测、粒径统计等任务。
---

# 数字全息仿真与重建管线 — 超级技能

## 项目结构

```
仿真/
├── functionme/
│   ├── 仿真图像生成/          # Python：合成全息图生成
│   │   ├── generate_hologram_matlab_validation_new.py  ← 核心生成模块
│   │   ├── generate_hologram_gui.py                    ← Tkinter GUI
│   │   ├── test_generate_hologram_outputs.py            ← 单元测试
│   │   ├── 大作业仿真.py                                ← 外差干涉仿真（独立）
│   │   └── 0001~0011/          # 已生成的数据集
│   │       ├── *_hologram.bmp           # 全息图 (8-bit)
│   │       ├── *_object_reference.bmp   # 物体截面参考
│   │       ├── *_particle_ground_truth.csv  # 粒子真值
│   │       └── *_background_noise.bmp  # 噪声背景（仅 0011）
│   │
│   ├── MATLAB版/              # 活跃的 MATLAB 重建代码
│   │   ├── VolumeGUI_AngularSpectrum.m          ← GUI 主入口 (1754行)
│   │   ├── Fresnel_reconstructFastStats.m       ← 方法1 核心 (580行)
│   │   ├── AS_preprocessHologram.m              ← 预处理（支持仿真模式）
│   │   ├── AS_normalizeImage.m                  ← 归一化 [0,1]
│   │   ├── AS_buildAdaptiveROI.m                ← 自适应搜索 ROI
│   │   ├── AS_refineMeasurementROI.m            ← 测量 ROI 细化
│   │   ├── AS_showParticleSpheres.m             ← 粒子 3D 可视化
│   │   ├── AS_reconstructAngularSpectrum.m      ← 方法2（备用）
│   │   ├── AS_reconstructAngularSpectrumFastStats.m ← 方法1 备选
│   │   ├── AS_localizeParticles3D.m             ← 3D 粒子定位
│   │   ├── show3d.m / vol3d.m                   ← 3D 体渲染
│   │   └── MATLAB代码说明.md                     ← 详细文档
│   │
│   ├── archive_matlab/        # 历史 MATLAB 代码（归档）
│   │   ├── cs_pipeline/       ← 压缩感知重建管线
│   │   ├── legacy_ui_and_tools/ ← 旧版 GUI 和工具
│   │   └── experiments/       ← 实验代码
│   │
│   ├── 仿真与MATLAB重建对齐问题.md  ← P1-P6 详细分析
│   └── 工作摘要_20260611.md        ← 本次工作摘要
│
└── .claude/
    ├── settings.local.json    ← 权限配置
    └── skills/
        └── hologram-pipeline.md ← 本文件
```

## 核心数据流

```
Python 生成端                    MATLAB 重建端
─────────────                    ─────────────
参数: λ, dx, N, d, z, SNR
    ↓
创建二元粒子掩模 (0/1)           VolumeGUI_AngularSpectrum
    ↓                               ├─ 参数面板 (z, λ, pixel)
角谱传播到探测器                    ├─ Simulation mode ☑
  H = exp(1j·k·z·√(1-λ²f²))        ├─ 加载 BMP
    ↓                               ├─ AS_preprocessHologram
|1 - U_total|²                      │   /255 → -bg → -mean
    ↓                               ├─ Fresnel_reconstructFastStats
_to_uint8 → BMP                     │   ├─ fft2 → ×H_Fresnel → ifft2
    ↓                               │   ├─ MIP 累积 (max over z)
CSV 真值                            │   ├─ iDetectCandidates
                                    │   ├─ 焦面粒径精测
                                    │   └─ 统计分析 (MVD, D[3,2])
                                    └─ 结果展示
```

## 关键参数对应

| Python 形参 | MATLAB UI | 默认 (Py → MATLAB) |
|---|---|---|
| `wavelength` | 波长 (nm) | 632.8 → 638 |
| `dx`, `dy` | 像素尺寸 (µm) | 2.2 → 3.45 |
| `N` | 取决于加载图 | 256 |
| `diameter_min/max` | (无直接控件) | 10-20 µm |
| `z_min/max` | z最小值/最大值 (mm) | 2-8 → 20-35 |
| `z_step` | z步长 (µm) | 40 → 20 |

## GUI 仿真模式工作原理

仿真模式复选框 → `simulationMode = true` 沿管线传递：

1. **AS_preprocessHologram**：`double(img)/255 → [0,1] → -mean`
2. **iFillDefaultParams**：`minDiamPx=3, maxDiamPx=200, minCircularity=0.3, invertContrast=false, watershedHmin=0.3`
3. **iDetectCandidates**：用 `params.minDiamPx/maxDiamPx/minCircularity` 替代硬编码
4. **核心传播**：`if invertContrast → I=max-Amp else I=Amp`

取消勾选 → 全部恢复实验默认值，原有流程不受影响。

## 常见问题速查

### 仿真全息图重建不出粒子？
1. ✅ 勾选 Simulation mode
2. ✅ 确认像元尺寸匹配 Python 的 dx
3. ✅ Z 范围覆盖 Python 的 z_min~z_max
4. ✅ 如果粒子很小（<10 µm），检查 xyThreshold 是否过高
5. ✅ 检查 `fastStatsParams` 中 `simulationMode` 是否正确传入

### 坐标对比偏移？
- Python CSV：0-indexed
- MATLAB coords3D：1-indexed
- 真值对比：`gt_x + 1` 对齐 MATLAB 坐标

### 传播模型差异影响？
- 角谱生成 → 菲涅尔重建：高 NA 时有亚像素偏差
- 实测影响小，主要问题在预处理和筛选环节

## 修改代码时注意事项

- **不要动实验模式的默认值**：`iFillDefaultParams` 的 else 分支保持 `[20,35,0.6]`
- **AS_preprocessHologram 的 options 参数可选**：省略时默认实验模式
- **所有仿真适配通过 `simulationMode` 标志集中控制**：不要散落 if-else
- **invertContrast 参数要传到所有使用点**：GPU 核心、CPU 核心、两个 collectPatches 函数
