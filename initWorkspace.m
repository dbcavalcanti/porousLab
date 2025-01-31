function initWorkspace()
    close all;
    % Clear the classes to avoid using a non-updated version
    clear Element RegularElement RegularElementPcPg;
    clear IntPoint;
    clear Shape Shape_CST Shape_LST Shape_ISOQ4 Shape_ISOQ8;
    clear Model IntPoint Result;
    clear Anl Anl_Linear Anl_Nonlinear Anl_Transient Anl_TransientPicard;
    clear CapillaryPressure CapillaryPressureBrooksCorey
    clear CapillaryPressureLiakopoulos  CapillaryPressureUMAT
    clear Fluid IdealGas
    clear MaterialTwoPhaseFlow
    clear PorousMedia
    clear RelativePermeability RelativePermeabilityBrooksCoreyGas
    clear RelativePermeabilityBrooksCoreyLiquid
    clear RelativePermeabilityLiakopoulosLiquid
    clear RelativePermeabilityUMAT
    clear EFEMDraw
    % Clear the workspace and the command window
    clear; clc;
    % porousLab header
    fprintf('===========================================================\n');
    fprintf('                        porousLab                          \n\n');
    fprintf('     FEM solver for multiphysics problems in porous media\n\n');
    fprintf('                        - CIMNE -                          \n');
    fprintf('  International Centre for Numerical Methods in Enginnering\n');
    fprintf('  Applied Computational Physics Group\n\n');
    fprintf('  Author: Danilo Cavalcanti\n');
    fprintf('  Last update: Jan/2025\n');
    fprintf('===========================================================\n\n');
    %Use all folders and subfolders
    addpath(genpath('./'));
end