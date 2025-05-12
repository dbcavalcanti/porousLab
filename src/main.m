%% PorousLab - FEM framework for multiphysics problems in porous media
%
% This is the main file of PorousLab.
% To run a simulation, execute this file and select an appropriate problem script.
%
clc; clearvars; close all;
addpath(genpath(pwd));
print_header;
[file_name, file_path] = uigetfile('*.m', 'Select a script to run');
if isequal(file_name, 0), return; end
addpath(file_path);
run(file_name);
rmpath(filepath);
