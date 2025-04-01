%% DESCRIPTION
%
% McWhorter and Sunada problem using the Pc-Pg two-phase flow formulation.
%
% References:
% * McWhorter and Sunada (1990). Exact integral solutions for two-phase flow. Water Resour Res, 26(3):399–413.
%
% Physics:
% * Two-phase hydraulic (H2)
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
mdl = Model_H2_PcPg();

% Set model options
mdl.massLumping  = true;  % Diagonalize compressibility matrix (mass lumping)
mdl.lumpStrategy = 2;
mdl.intOrder     = 3;     % Integration rule order for the domain

%% MESH

% Create mesh
Lx = 2.6;  % Horizontal dimension (m)
Ly = 0.5;  % Vertical dimension (m)
Nx = 260;  % Number of elements in the x-direction
Ny = 1;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
gas   = Fluid('gas');

% Create porous media
rock = PorousMedia('rock');
rock.K                  = 1.0e-10;        % Intrinsic permeability (m2)
rock.phi                = 0.3;            % Porosity
rock.Slr                = 0.0;            % Residual liquid saturation
rock.Sgr                = 0.0;            % Residual gas saturation
rock.Pb                 = 5.0e+3;         % Gas-entry pressure
rock.lambda             = 2.0;            % Curve-fitting parameter
rock.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
rock.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
rock.capillaryPressure  = 'BrooksCorey';  % Saturation degree function

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Set Dirichlet boundary conditions
mdl.setCapillaryPressureDirichletBCAtBorder('left', 5.0e+3);
mdl.setGasPressureDirichletBCAtBorder('left', 200.0e+3);

% Set initial conditions
mdl.setInitialCapillaryPressureAtDomain(50.0e+3);

%% PROCESS

% Analysis parameters
ti        = 0.1;        % Initial time
dt        = 0.1;        % Time step
tf        = 50.0;       % Final time
dtmax     = 0.1;        % Maximum time step
dtmin     = 0.0000001;  % Minimum time step
adaptStep = true;       % Adaptive step size

% Run analysis
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.run(mdl);

%% POST-PROCESS

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, 0.0];
mdl.plotPressureAlongSegment(Xi, Xf, 500, 'x');
