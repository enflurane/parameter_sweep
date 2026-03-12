clc;
clear all; 

% build_bifdiag_mex.m - 构建分岔图CUDA MEX函数

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);

% CUDA源文件路径
cu_file = fullfile(this_dir,'..','..','cuda','bifdiag','bifdiag_mex.cu');

% 检查CUDA文件是否存在
if ~exist(cu_file, 'file')
    error('CUDA源文件不存在: %s', cu_file);
end

fprintf('正在编译分岔图CUDA MEX函数...\n');
fprintf('CUDA源文件: %s\n', cu_file);

% 编译CUDA MEX函数
try
    mexcuda('-v', cu_file, ...
        'NVCCFLAGS=-arch=sm_86 --generate-code=arch=compute_86,code=sm_86 -Wno-deprecated-gpu-targets', ...
        '-lcudart', ...
        '-output','bifdiag_mex');
    
    fprintf('分岔图CUDA MEX函数编译成功！\n');
    fprintf('生成文件: bifdiag_mex.mexw64\n');
    
catch ME
    fprintf('编译失败！错误信息:\n');
    fprintf('%s\n', ME.message);
    
    % 提供调试建议
    fprintf('\n调试建议:\n');
    fprintf('1. 检查CUDA Toolkit是否已正确安装\n');
    fprintf('2. 检查MATLAB编译器配置: mex -setup C++\n');
    fprintf('3. 检查GPU计算能力是否匹配 (需要sm_86或更高)\n');
    fprintf('4. 检查CUDA路径: C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.8\\\n');
    
    rethrow(ME);
end