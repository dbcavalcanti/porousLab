function initWorkspace()
    % Clear the workspace and the command window
    close all;
    clear;
    clc;

    % PorousLab header
    fprintf('===========================================================\n');
    fprintf('                        porousLab                          \n\n');
    fprintf('     FEM solver for multiphysics problems in porous media\n\n');
    fprintf('                        - CIMNE -                          \n');
    fprintf('  International Centre for Numerical Methods in Enginnering\n');
    fprintf('  Applied Computational Physics Group\n\n');
    fprintf('  Author: Danilo Cavalcanti (dborges@cimne.upc.edu)\n');
    fprintf('  Last update: Feb/2025\n');
    fprintf('===========================================================\n\n');

    % Use all folders and subfolders
    addpath(genpath('./'));
end