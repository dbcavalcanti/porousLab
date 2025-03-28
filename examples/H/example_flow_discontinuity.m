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

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

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

mdl.setPressureDirichletBCAtBorder('left',0.0);
mdl.setPressureDirichletBCAtBorder('right',10.0);

%% PRE-PROCESS

% Create the discontinuity elements
fracture.intersectMesh(mdl);

% Add the fracture to the model
mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

% Solve the problem
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Plot pressure along a segment
mdl.printResults();

mdl.plotField('Pressure'); hold on
fracture.plotIntersectedGeometry()
