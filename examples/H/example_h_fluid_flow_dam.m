%% DESCRIPTION
%
% Fluid-flow through the foundation of a gravity dam.
%
% References:
% * Segura and Carol (2004). On zero-thickness interface elements for diffusion problems. Int J Numer Anal Methods Geomech, 28(9):947-962.
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
Lx = 24.0;  % Horizontal dimension (m)
Ly = 6.0;   % Vertical dimension (m)
Nx = 96/2;  % Number of elements in the x-direction
Ny = 24/2;  % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.K = 2.0e+9;  % Compressibility/Bulk modulus (1/Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K   = 1.0194e-14;  % Intrinsic permeability (m2)
rock.phi = 0.3;         % Porosity

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY CONDITIONS

% Set Dirichlet boundary conditions
for i = 1:size(mdl.NODE,1)
    if ((mdl.NODE(i,1) < 8.0) && (abs(mdl.NODE(i,2)-Ly) < 1.0e-9))
        mdl.setPressureDirichletBCAtNode(i, 120.0);
    end
    if ((mdl.NODE(i,1) > 12.0) && (abs(mdl.NODE(i,2)-Ly) < 1.0e-9))
        mdl.setPressureDirichletBCAtNode(i, 60.0);
    end
end

%% PROCESS

% Analysis parameters
ti        = 1.0;    % Initial time
dt        = 1.0;    % Time step
tf        = 500.0;  % Final time
dtmax     = 50.0;   % Maximum time step
dtmin     = 0.001;  % Minimum time step
adaptStep = true;   % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot contours
mdl.plotField('Pressure');

% Plot graphs
Xi = [0.0, Ly/2.0]; Xf = [Lx, Ly/2.0];
mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'x');
