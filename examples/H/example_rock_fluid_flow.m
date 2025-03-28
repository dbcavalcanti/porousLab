%% DESCRIPTION
%
% ...
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

% Mesh generation
Lx = 2.0;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 2;    % Number of elements in the x-direction
Ny = 1;    % Number of elements in the y-direction

[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

%% MATERIAL CREATION

% Create fluids
water = Fluid('water');
water.rho = 1000.0;  % Density (kg/m3)
water.mu  = 1.0e-3;  % Viscosity (Pa*s)
water.K   = 2.0e9;   % Compressibility/Bulk modulus (1/Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K    = 1.0194e-14;  % Intrinsic permeability (m2)
rock.phi  = 0.3;         % Porosity
rock.Ks   = 1.0e12;      % Rock bulk modulus (Pa)
rock.biot = 0.6;         % Biot coefficient

% Material parameters vector
mdl.mat = struct('porousMedia',rock,'fluid',water);

%% BOUNDARY CONDITIONS

mdl.setPressureDirichletBCAtBorder('left',0.0);
mdl.setPressureDirichletBCAtBorder('right',10.0);

%% PROCESS

% Transient analysis parameters
tinit = 1.0;   % Initial time
dt    = 1.0;   % Time step
tf    = 500;   % Final time
dtmax = 1.0;
dtmin = 1.0;

% Solve the problem
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(tinit,dt,tf,50.0,0.001,true);
anl.setRelativeConvergenceCriteria(true);
anl.run(mdl);

%% POST-PROCESS

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , Ly/2.0];
Xf  = [Lx , Ly/2.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('Pressure')
