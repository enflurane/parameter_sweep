clc;
clear all; 

% build_sweep_mex.m

this_file = mfilename('fullpath');
this_dir = fileparts(this_file);

cu_file = fullfile(this_dir,'..','cuda','bisweep','sweep_mex.cu');

mexcuda('-v', cu_file, ...
    'NVCCFLAGS=-arch=sm_86 --generate-code=arch=compute_86,code=sm_86 -Wno-deprecated-gpu-targets', ...
    '-lcudart', ...
    '-output','sweep_mex');