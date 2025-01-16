function initWorkspace()
    close all;
    % Clear the classes to avoid using a non-updated version
    clear Element RegularElement  EnrichedElement;
    clear Fracture Fracture_ConstantJump Fracture_LinearJump;
    clear IntPoint;
    clear Shape Shape_CST Shape_LST Shape_ISOQ4 Shape_ISOQ8;
    clear Model IntPoint Result;
    clear Anl Anl_Linear Anl_Nonlinear Anl_Transient;
    clear MaterialHydro MaterialHydro_Saturated
    clear MaterialHydroInterface MaterialHydroInterface_CubicLaw
    % Clear the workspace and the command window
    clear; clc;
    %Use all folders and subfolders
    addpath(genpath('./'));
end