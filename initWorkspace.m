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
    %Use all folders and subfolders
    addpath(genpath('./'));
end