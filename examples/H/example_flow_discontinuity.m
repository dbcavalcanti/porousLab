%% DESCRIPTION
%
% Block crossed by a strong discontinuity.
%
% Physics:
% * Single-phase hydraulic (H)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% INITIALIZATION
close all; clear; clc;

% Path to source directory
src_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
addpath(genpath(src_dir));
print_header;

% Create model
mdl = Model_H();

%% MESH GENERATION

% Mesh geometry
Lx = 5.0;  % Horizontal dimension (m)
Ly = 3.0;  % Vertical dimension (m)
Nx = 5;    % Number of elements in the x-direction
Ny = 3;    % Number of elements in the y-direction

[mdl.NODE, mdl.ELEM] = regularMesh(Lx, Ly, Nx, Ny);

% Discontinuity generation
Dx = [1.0; 4.0];  % X-coordinates of polyline defining the fracture
Dy = [1.1; 1.9];  % Y-coordinates of polyline defining the fracture

fracture = Discontinuity([Dx, Dy], true);

%% MATERIAL CREATION

% Create fluids
water = Fluid('water');
water.rho = 1000.0;  % Density (kg/m3)
water.mu  = 1.0e-3;  % Viscosity (Pa*s)
water.K   = 2.0e9;   % Compressibility/Bulk modulus (1/Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K   = 9.8e-16;  % Intrinsic permeability (m2)
rock.phi = 0.25;     % Porosity

% Material parameters vector
mdl.mat = struct('porousMedia',rock,'fluid',water);

% Set fracture material properties
fracture.initialAperture = 1.0e-3;
fracture.fluid = water;

%% BOUNDARY CONDITIONS

% Pore pressure
CoordSupp  = [1 0.0 -1; 1 Lx -1];         
CoordLoad  = [];            
CoordPresc = [0.0 0.0 -1; 10.0 Lx -1];            
CoordInit  = [];                   
           
% Supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] =...
boundaryConditionsPressure(mdl.NODE, CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

%% PRE-PROCESS

% Create the discontinuity elements
fracture.intersectMesh(mdl);

% Add the fracture to the model
mdl.addPreExistingDiscontinuities(fracture);

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% PROCESS

% Solve the problem
anl = Anl_Linear(result);
anl.process(mdl);

%% POST-PROCESS

% Plot pressure along a segment
mdl.printResults();

mdl.plotField('Pressure'); hold on
fracture.plotIntersectedGeometry()
