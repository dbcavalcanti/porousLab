%% DESCRIPTION
%
% Block crossed by strong discontinuity.
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

%% MODEL

% Create model
mdl = Model_H();

%% MESH

% Create mesh
[node, elem] = regularMesh(200.0, 200.0, 53, 53);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.rho = 1.0e+3;  % Density (kg/m3)
water.mu  = 1.0e-3;  % Viscosity (Pa*s)
water.K   = 2.2e+9;  % Compressibility/Bulk modulus (1/Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K   = 1.0e-14;  % Intrinsic permeability (m2)
rock.phi = 0.25;     % Porosity

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY CONDITIONS

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('bottom', 0.0);
mdl.setPressureDirichletBCAtBorder('top', 1000000.0);

%% DISCONTINUITIES

% loads fracture_data
load('FractureDataReservoirCell.mat'); 

% Number of discontinuities
nd = length(FractureDataReservoirCell);

% Create the discontinuities
fractures(1, nd) = Discontinuity();
for i = 1:length(FractureDataReservoirCell)
    fractures(i) = Discontinuity(FractureDataReservoirCell{i}, true);
end

% Set fracture material properties
for i = 1:nd
    fractures(i).fluid = water;
    fractures(i).initialAperture = 1.0e-3;
end

% Add fractures to model
mdl.addPreExistingDiscontinuities(fractures);

%% PROCESS

% Analysis parameters
ti        = 1.0;    % Initial time
dt        = 1.0;    % Time step
tf        = 1000.0;  % Final time
dtmax     = 500.0;  % Maximum time step
dtmin     = 0.001;  % Minimum time step
adaptStep = true;   % Adaptive step size

% Initialize
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);

% Clear unnecessary variables
clearvars -except mdl anl

% Run analysis
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot model
mdl.plotField('Model');
colorbar off; hold on;
for i = 1:length(mdl.discontinuitySet)
    mdl.discontinuitySet(i).plotIntersectedGeometry();
end

% Plot contours
mdl.plotField('Pressure'); 
