# Bifdiag - CUDA加速的分岔图计算模块

基于CUDA GPU加速的Qi系统分岔图计算模块，用于绘制非线性动力学系统的双参数分岔图。

## 功能特点

- **GPU并行计算**：利用CUDA实现大规模参数空间的高效扫描
- **双参数分岔图**：同时扫描两个参数，生成二维分岔图
- **Poincaré截面**：通过Poincaré截面方法分析周期轨道
- **自动分段计算**：支持大范围参数扫描的分段处理和结果拼接
- **完整工作流**：从计算、保存到绘图的一站式解决方案

## 文件结构

```
bifdiag/
├── build_bifdiag_mex.m    # 编译CUDA MEX函数的构建脚本
├── main.m                 # 主程序：分岔图计算和结果保存
├── bifdiag_mex.mexw64     # 编译后的CUDA MEX函数（Windows）
└── README.md              # 本文档

../cuda/bifdiag/
└── bifdiag_mex.cu         # CUDA内核源代码
```

## 系统要求

- **操作系统**: Windows 10/11
- **MATLAB**: R2022a 或更高版本
- **CUDA Toolkit**: 12.8
- **GPU**: NVIDIA GPU (计算能力 sm_86 或更高，建议 RTX 40/30 系列)
- **显存**: 至少 6GB，推荐 8GB+
- **编译器**: Visual Studio 2019 (MSVC v142)

## 快速开始

### 1. 编译CUDA MEX函数

确保已安装所有依赖后，在MATLAB中运行：

```matlab
cd 'F:\money\parameter_sweep\src\matlab\bifdiag'
build_bifdiag_mex
```

如果编译成功，将生成 `bifdiag_mex.mexw64` 文件。

### 2. 运行分岔图计算

```matlab
main
```

默认计算参数：
- 参数 a: 34 ~ 35，分3段计算，每段101点
- 参数 c: -50 ~ 50，101点
- 固定参数: b=1, c_sys=104.5
- 积分: dt=0.001, N=1e7, kneadings=1000~1050

### 3. 查看结果

计算完成后，会在当前目录生成类似 `bif_b1_c104.5_a34.00-35.00_20260312_154312` 的时间戳文件夹，包含：
- `P.mat` - 分岔图数据 (P_all, a_all, 参数)
- `plot_data.m` - 独立绘图脚本
- `bifurcation_diagram.png` - 分岔图图像

## 参数配置

在 `main.m` 中修改以下参数：

### 扫描参数
```matlab
% 参数a扫描（分岔参数）
a_full_start = 34;        % 起始值
a_full_end = 35;          % 结束值
a_step = 0.1;             % 分段跨度（小于总范围时分段计算）
points_per_block = 101;   % 每段点数

% 参数c扫描（Poincaré截面）
parameter2Start = -50;    % 起始值
parameter2End = 50;       % 结束值
parameter2Count = 101;    % 采样点数
```

### 系统参数
```matlab
parameter2 = 1;           % 固定参数 b
parameter3 = 104.5;       % 固定参数 c (系统参数)
whichSweep = 0;           % 扫描类型: 0=a, 1=b, 2=c
```

### 积分设置
```matlab
dt = 0.001;              % 时间步长
N = 1e7;                 % 总迭代次数
stride = 1;              % 采样间隔
kneadingsStart = 1000;   % Poincaré采样起始迭代
kneadingsEnd = 1050;     % Poincaré采样结束迭代
```

### 初始条件
```matlab
x0_initial = 0.1;
y0_initial = 0.1;
z0_initial = 0.1;
```

## 算法说明

### Qi系统方程
```
dx/dt = a(y - x) + yz
dy/dt = cx - y - xz
dz/dt = xy - bz
```

### 计算流程

1. **初始化**：根据当前参数计算平衡点距离，设置初始条件
2. **积分**：使用4阶Runge-Kutta方法进行数值积分
3. **Poincaré截面**：当轨迹穿过固定点定义的截面时记录值
4. **统计**：对每个参数点计算Poincaré截面的平均值
5. **可视化**：绘制参数-截面值散点图

### 分段计算

当 `a_step < (a_full_end - a_full_start)` 时启用分段计算，避免单次计算内存不足或超时。各段结果自动横向拼接为完整数据。

## 输出说明

### 数据文件 (P.mat)
- `P_all` - 分岔图数据矩阵 (parameter2Count × total_points)
- `a_all` - 对应的参数a值向量
- 所有输入参数（parameter2Start, parameter2End, parameter2Count 等）

### 绘图脚本 (plot_data.m)
独立脚本，可单独运行重新绘图：
```matlab
load('P.mat');
% 自动加载数据并绘制分岔图
```

### 图像文件
- `bifurcation_diagram.png` - PNG格式分岔图
- `bifurcation_diagram.fig` - MATLAB图形文件（可编辑）

## 性能优化建议

### 1. 调整参数点数
```matlab
% 降低分辨率以提高速度
parameter2Count = 51;    % 从101降低到51
points_per_block = 51;   % 减少每段点数
```

### 2. 减少积分迭代
```matlab
N = 5e6;                % 减少迭代次数（可能影响精度）
kneadingsEnd = kneadingsStart + 30;  % 减少Poincaré采样点数
```

### 3. 增大时间步长
```matlab
dt = 0.002;             % 增大时间步长（降低精度）
```

### 4. 启用分段计算
```matlab
a_step = 0.05;          % 减小分段跨度，每段计算量更小
```

## 故障排除

### 编译错误

**错误**: `CUDA源文件不存在`
- 检查 `bifdiag_mex.cu` 文件路径是否正确
- 确认 `cuda/bifdiag/` 目录存在

**错误**: `Cannot find compiler 'cl.exe'`
- 运行 `mex -setup C++` 配置MATLAB编译器
- 确保Visual Studio 2019已安装并添加到PATH

**错误**: `Unsupported GPU architecture`
- 检查GPU计算能力，修改 `build_bifdiag_mex.m` 中的 `-arch=sm_XX`
- RTX 40系列: sm_89, RTX 30系列: sm_86

### 运行时错误

**错误**: `CUDA error: out of memory`
- 减少 `parameter2Count` 或 `points_per_block`
- 使用分段计算（减小 `a_step`）

**错误**: 分岔图空白或异常
- 检查参数范围是否合理（某些参数组合可能无周期轨道）
- 增加 `N` 和 `kneadingsEnd` 以获得更稳定的统计
- 调整初始条件 `x0_initial, y0_initial, z0_initial`

**警告**: `expression has no effect`
- 来自CUDA代码的空语句，可忽略

## 示例：快速测试

创建测试脚本 `quick_test.m`：

```matlab
% 快速测试 - 低分辨率
a_full_start = 34.5;
a_full_end = 34.6;
a_step = 0.1;
points_per_block = 21;

parameter2Start = -20;
parameter2End = 20;
parameter2Count = 21;

N = 1e6;  % 减少迭代加速测试

main;
```

预期耗时：1-5分钟（取决于GPU性能）

## 高级用法

### 修改系统方程

编辑 `cuda/bifdiag/bifdiag_mex.cu` 中的 `stepper` 函数：

```cuda
__device__ void stepper(const double* y, double* dydt, const double* params)
{
    double a = params[0], b = params[1], c = params[2];
    
    // 修改为其他系统，例如 Lorenz:
    // dydt[0] = a * (y[1] - y[0]);
    // dydt[1] = y[0] * (params[1] - y[2]) - y[1];
    // dydt[2] = y[0] * y[1] - params[2] * y[2];
    
    // Qi系统（默认）:
    dydt[0] = a*y[1] - a*y[0] + y[1]*y[2];
    dydt[1] = c*y[0] - y[1] - y[0]*y[2];
    dydt[2] = y[0]*y[1] - b*y[2];
}
```

重新编译后即可计算其他系统的分岔图。

### 批量计算多个参数组合

```matlab
% 保存原始参数
orig_b = parameter2;
orig_c = parameter3;

% 循环计算不同参数
for b_val = [0.5, 1.0, 1.5, 2.0]
    parameter2 = b_val;
    main;
end

% 恢复参数
parameter2 = orig_b;
parameter3 = orig_c;
```

## 引用

如果本项目对您的研究有帮助，请引用相关论文：

- Qi系统分岔分析: [Complex bifurcations and new types of structure uncovered in the Qi system](https://www.sciencedirect.com/science/article/pii/S0378475425005440)

## 许可证

MIT License

