clc;
clear all; 

addpath('../../../expake/slanCM');

% 参数设置（与原始C代码main函数参数对应）

%a
parameter1Start = 34; 
parameter1End = 35; 
parameter1Count = 101; 

%c
parameter2Start = 99; 
parameter2End= 101; 
parameter2Count = 101; 

%b
parameter3 = 1;         % 固定参数3

whichSweep = 1;         % 扫描类型选择

%初始点
x0_initial=0.1;
y0_initial=0.1;
z0_initial=0.1;

dt  = 0.01;              % 时间步长
N  = 1e6;              % 迭代次数
stride = 1;             % 采样间隔

kneadingsStart = 1000;     % kneadings区间起点
kneadingsEnd = 2023;      % kneadings区间终点

% 启动CUDA核函数
tic;
P = sweep_mex(...
    parameter1Start, parameter1End, parameter1Count,...
    parameter2Start, parameter2End, parameter2Count,...
    parameter3, whichSweep,...
    x0_initial,y0_initial,z0_initial,...
    dt, N, stride,...
    kneadingsStart, kneadingsEnd);
fprintf('计算耗时: %.2f 秒\n', toc);

% 结果可视化
figure
a = linspace(parameter1Start,parameter1End, parameter1Count);
c = linspace(parameter2Start, parameter2End, parameter2Count);
imagesc(a, c, P); 
hold on;
slan = slanCM(193);
colormap(slan);
caxis([-0.99 1])
set(gca, 'LooseInset', [0,0,0,0], 'FontSize', 25);
set(gca, 'YDir', 'normal');
pbaspect([1 1 1]);   
xlabel('$a$-parameter', 'Interpreter', 'latex', 'FontSize', 25);
ylabel('$c$-parameter', 'Interpreter', 'latex', 'FontSize', 25);

% 时间戳
timestamp = datestr(now,'yyyymmdd_HHMMSS');

% 生成文件夹名称并创建
folder_name = sprintf('sweep_b%g_a%.2f-%.2f_c%.2f-%.2f_%s', ...
    parameter3, ...
    parameter1Start, parameter1End, ...
    parameter2Start, parameter2End, ...
    timestamp);

if ~exist(folder_name,'dir')
    mkdir(folder_name);
end

% 保存数据到文件夹
save(fullfile(folder_name, 'P.mat'), 'P');

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
fprintf(fid, '%% 设置参数\n');
fprintf(fid, 'parameter1Start = %g; \n',parameter1Start);
fprintf(fid, 'parameter1End = %g; \n',parameter1End);
fprintf(fid, 'parameter1Count = %g; \n\n',parameter1Count);
fprintf(fid, 'parameter2Start = %g; \n',parameter2Start);
fprintf(fid, 'parameter2End= %g; \n',parameter2End);
fprintf(fid, 'parameter2Count = %g; \n\n',parameter2Count);
fprintf(fid, 'parameter3  = %g; \n\n',parameter3 );
fprintf(fid, 'whichSweep  = %g; \n\n',whichSweep );
fprintf(fid, 'x0_initial  = %g; \n',x0_initial );
fprintf(fid, 'y0_initial  = %g; \n',y0_initial );
fprintf(fid, 'z0_initial  = %g; \n',z0_initial );
fprintf(fid, 'dt  = %g; \n',dt );
fprintf(fid, 'N  = %g; \n',N );
fprintf(fid, 'stride  = %g; \n\n',stride );
fprintf(fid, 'kneadingsStart = %g; \n',kneadingsStart );
fprintf(fid, 'kneadingsEnd = %g; \n\n',kneadingsEnd );
fprintf(fid, 'a = linspace(%g, %g, %d);\n', parameter1Start, parameter1End, parameter1Count);
fprintf(fid, 'c = linspace(%g, %g, %d);\n\n', parameter2Start, parameter2End, parameter2Count);
fprintf(fid, '%% 绘制图像\n');
fprintf(fid, 'figure;\n');
fprintf(fid, 'imagesc(a, c, P);\n');
fprintf(fid, 'set(gca, ''YDir'', ''normal'');\n');
fprintf(fid, 'colormap(slanCM(193));\n');
fprintf(fid, 'caxis([-0.99 1]);\n');
fprintf(fid, 'xlabel(''$a$-parameter'', ''Interpreter'', ''latex'', ''FontSize'', 25);\n');
fprintf(fid, 'ylabel(''$c$-parameter'', ''Interpreter'', ''latex'', ''FontSize'',25);\n');
fprintf(fid, 'pbaspect([1 1 1]);\n');
fprintf(fid, 'set(gca, ''LooseInset'', [0,0,0,0], ''FontSize'', 25);\n');
fclose(fid);