# 分岔图模块 (Bifdiag)

基于CUDA GPU加速的Qi系统分岔图计算模块，用于绘制非线性动力学系统的分岔图。

## 功能概述

本模块实现了：
- **CUDA加速的分岔图计算**：利用GPU并行计算大幅提升计算效率
- **Poincaré截面分析**：通过截面方法捕捉系统的周期轨道
- **分段计算**：支持大范围参数扫描的分段计算和拼接
- **自动结果管理**：自动生成带时间戳的结果文件夹和绘图脚本

## 系统要求

### 硬件要求
- NVIDIA GPU (计算能力sm_86或更高)
- 至少6GB显存

### 软件要求
- MATLAB R2022a或更高版本
- CUDA Toolkit 12.8
- Visual Studio 2019 (MSVC v142)
- Windows 10/11

## 快速开始

### 1. 编译CUDA MEX函数

```matlab
cd src/matlab/bifdiag
build_bifdiag_mex
```

### 2. 运行分岔图计算

```matlab
main
```

### 3. 查看结果

结果将自动保存在以时间戳命名的文件夹中，包含：
- `P.mat` - 原始数据
- `plot_data.m` - 绘图脚本
- `bifurcation_diagram.png` - 分岔图图像
- `bifurcation_diagram.fig` - MATLAB图形文件

## 参数说明

### 主要参数
- `parameter1Start`, `parameter1End`, `parameter1Count` - 第一个扫描参数（通常是a）的范围
- `parameter2Start`, `parameter2End`, `parameter2Count` - 第二个扫描参数（通常是c）的范围
- `parameter2` - 固定参数b
- `parameter3` - 固定参数c（系统参数）
- `whichSweep` - 扫描类型（0:扫描a，1:扫描b，2:扫描c）

### 积分参数
- `dt` - 时间步长（默认0.001）
- `N` - 总迭代次数（默认1e7）
- `stride` - 采样间隔（默认1）
- `kneadingsStart`, `kneadingsEnd` - Poincaré截面采样区间

### 初始条件
- `x0_initial`, `y0_initial`, `z0_initial` - 系统初始状态

### 分段计算参数
- `a_full_start`, `a_full_end` - 总参数范围
- `a_step` - 每段跨度（小于总范围时启用分段计算）
- `points_per_block` - 每段点数

## 算法原理

### 1. Qi系统动力学方程
```
dx/dt = a(y - x) + yz
dy/dt = cx - y - xz
dz/dt = xy - bz
```

### 2. Poincaré截面方法
通过定义截面，记录轨迹穿过截面的点，从而将连续动力系统简化为离散映射。

### 3. 龙格-库塔积分
使用4阶龙格-库塔方法进行数值积分：
```
k1 = f(y_n)
k2 = f(y_n + dt/2 * k1)
k3 = f(y_n + dt/2 * k2)
k4 = f(y_n + dt * k3)
y_{n+1} = y_n + dt/6 * (k1 + 2k2 + 2k3 + k4)
```

### 4. CUDA并行化
- 每个线程处理一个参数点
- 并行计算不同参数条件下的轨迹
- 使用共享内存优化数据访问

## 输出格式

### 数据文件 (P.mat)
- `P_all` - 分岔图数据矩阵（size: parameter2Count × total_param1_points）
- `a_all` - 对应的参数a值（或其他扫描参数）
- 所有输入参数

### 图像输出
- 散点图形式的分岔图
- 横轴：扫描参数（通常是a）
- 纵轴：第二个参数（通常是c）
- 颜色：Poincaré截面值
- 使用slanCM配色方案

## 使用示例

### 基本使用
```matlab
% 修改参数范围
parameter1Start = 30;
parameter1End = 40;
parameter1Count = 201;

% 运行计算
main;
```

### 扫描不同参数
```matlab
% 扫描参数b
whichSweep = 1;
parameter2 = 2.5;  % 固定参数a
parameter3 = 104.5; % 固定参数c
parameter1Start = 0.5;
parameter1End = 3;
parameter1Count = 101;
main;
```

### 调整精度
```matlab
% 高精度计算
dt = 0.0005;
N = 2e7;
kneadingsStart = 2000;
kneadingsEnd = kneadingsStart + 100;
main;
```

## 性能优化建议

### 1. 参数选择
- 适当减少`parameter1Count`和`parameter2Count`以提高速度
- 增大`dt`可加快计算但可能降低精度
- 减少`N`可缩短计算时间但可能无法达到稳态

### 2. 内存管理
- 大范围参数扫描建议使用分段计算（设置`a_step`）
- 监控GPU显存使用情况
- 适当调整`kneadingsEnd - kneadingsStart`控制数据量

### 3. 精度控制
- 减小`dt`提高积分精度
- 增加`N`确保达到稳态
- 调整`kneadingsStart`跳过瞬态过程

## 常见问题

### 编译错误
1. **CUDA Toolkit未找到**：检查CUDA安装路径
2. **计算能力不匹配**：更新`build_bifdiag_mex.m`中的`-arch=sm_xx`
3. **编译器配置错误**：运行`mex -setup C++`

### 运行时错误
1. **GPU内存不足**：减少参数点数或使用分段计算
2. **数值不稳定**：减小时间步长`dt`
3. **无结果输出**：检查`kneadingsStart`是否设置过大

### 结果异常
1. **分岔图空白**：调整参数范围
2. **图像噪声过大**：增加`N`或调整`kneadings`区间
3. **周期轨道不清晰**：优化Poincaré截面定义

## 算法细节

### Poincaré截面计算
代码使用了基于固定点距离的Poincaré截面：
```cuda
DerivativeCurrent = fixedPointDistance[0] * y_current[1] - fixedPointDistance[1] * y_current[0];
```
当导数符号变化时，记录截面交点，并计算对应的值。

### 固定点距离计算
```cuda
computeFixedPointDistance(fixedPointDistance, params);
```
计算平衡点与当前状态的距离，用于定义Poincaré截面。

### 线程映射
```cuda
int i = tx / parameter2Step;
int j = tx % parameter2Step;
int flat_index = i * parameter2Step + j;
```
将一维线程ID映射到二维参数网格。

## 扩展功能

### 支持其他动力系统
修改`stepper`函数中的微分方程：
```cuda
__device__ void stepper(const double* y, double* dydt, const double* params)
{
    // 替换为其他系统的方程
    // 例如Lorenz系统:
    // dydt[0] = params[0] * (y[1] - y[0]);
    // dydt[1] = y[0] * (params[1] - y[2]) - y[1];
    // dydt[2] = y[0] * y[1] - params[2] * y[2];
}
```

### 添加其他分析方法
- Lyapunov指数计算
- 功率谱分析
- 返回映射分析
- 混沌度量化

## 参考文献

1. **Qi系统分岔分析**: [Complex bifurcations and new types of structure uncovered in the Qi system](https://www.sciencedirect.com/science/article/pii/S0378475425005440)
2. **Poincaré截面方法**: Guckenheimer & Holmes, *Nonlinear Oscillations, Dynamical Systems, and Bifurcations of Vector Fields*
3. **CUDA并行计算**: NVIDIA CUDA Programming Guide

## 许可证

MIT License - 详见项目根目录LICICENSE文件

## 更新日志

### v1.0.0 (2026-03-12)
- 初始版本发布
- 基于原始bif_qi.cu代码重构
- 参照bisweep项目结构设计
- 完整的文档和示例