clc;
clear all; 

% 分岔图主程序 - 基于CUDA加速的分岔图计算
% 参照bisweep项目结构设计

addpath('../../../expake/slanCM');

% 总区间设定
a_full_start = 34;
a_full_end = 35;
a_step = 0.1;                  % 每块跨度
points_per_block = 101;     % 每块点数（横轴）

% c 参数固定
parameter2 = 1;
parameter3 = 104.5;
whichSweep = 0;

x0_initial = 0.1;
y0_initial = 0.1;
z0_initial = 0.1;

parameter2Start = -50;
parameter2End = 50;
parameter2Count = 101;

dt = 0.001;
N = 1e7;
stride = 1;

kneadingsStart = 1000;
kneadingsEnd = kneadingsStart+50;

% 准备拼接
P_all = [];
a_all = [];
tic;
% 迭代多个小区间
for a_start = a_full_start : a_step : (a_full_end - a_step)
    a_end = a_start + a_step;

    fprintf('正在计算区间 [%.2f, %.2f]...\n', a_start, a_end);
    
    
    P = bifdiag_mex( ...
        a_start, a_end, points_per_block, ...
        parameter2Start, parameter2End, parameter2Count, ...
        parameter2, ...
        parameter3, whichSweep, ...
        x0_initial, y0_initial, z0_initial, ...
        dt, N, stride, ...
        kneadingsStart, kneadingsEnd );

    
    % 拼接数据
    P_all = [P_all, P];  % 横向拼接
    a_range = linspace(a_start, a_end, points_per_block);
    a_all = [a_all, a_range];  % 保存对应的 a 参数值
end
fprintf('计算耗时: %.2f 秒\n', toc);

% 时间戳
timestamp = datestr(now,'yyyymmdd_HHMMSS');

% 生成文件夹名称并创建
folder_name = sprintf('bif_b%g_c%g_a%.2f-%.2f_%s', ...
    parameter2, parameter3, ...
    a_full_start,a_full_end, ...
    timestamp);

if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

% 保存数据到文件夹（包括 a_all 和 P_all）
save(fullfile(folder_name, 'P.mat'), 'P_all', 'a_all', 'parameter2Start', 'parameter2End', 'parameter2Count');

% 生成绘图脚本
script_name = 'plot_data.m';
script_path = fullfile(folder_name, script_name);
fid = fopen(script_path, 'w');

fprintf(fid, '%% 自动生成的绘图脚本\n');
fprintf(fid, 'clc;\n');
fprintf(fid, 'clear;\n');
fprintf(fid, '%% 添加slanCM路径\n');
fprintf(fid, 'addpath(''../../../../expake/slanCM'');\n\n');
fprintf(fid, '%% 加载数据\n');
fprintf(fid, 'load(''P.mat'');\n\n');

fprintf(fid, '%% 构造绘图网格\n');
fprintf(fid, 'c_values = linspace(%g, %g, %d);\n', parameter2Start, parameter2End, parameter2Count);
fprintf(fid, '[X, C] = meshgrid(a_all, c_values);\n\n');

fprintf(fid, '%% 绘制散点图\n');
fprintf(fid, 'figure;\n');
fprintf(fid, 'scatter(X(:), P_all(:), 1, ''.'');\n');
fprintf(fid, 'xlabel(''Parameter a'', ''FontSize'', 16);\n');
fprintf(fid, 'ylabel(''Value in Matrix P'', ''FontSize'', 16);\n');
fprintf(fid, 'title(''分岔图'', ''FontSize'', 18);\n');
fprintf(fid, 'set(gca, ''FontSize'', 14);\n');
% fprintf(fid, 'grid on;\n');
% fprintf(fid, 'box on;\n');
% fprintf(fid, 'ylim([%g %g]);\n', kneadingsStart, kneadingsEnd);  % 可调整显示范围

fclose(fid);


% 构造绘图数据
c_values = linspace(parameter2Start, parameter2End, parameter2Count);
[X, C] = meshgrid(a_all, c_values);  % 行：c，列：a

% 绘制散点图
scatter(X(:), P_all(:), 1, 'k.');  % 横坐标：a，纵坐标：P值
xlabel('Parameter a');
ylabel('Value in Matrix P');
title('分岔图');
% ylim([60 180])
