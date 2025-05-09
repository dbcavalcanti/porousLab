%% DESCRIPTION
%
% Five-spot corner problem using the Pc-Pg two-phase flow formulation.
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

%% MESH

% Create mesh
Lx = 10.0;  % Horizontal dimension (m)
Ly = 10.0;  % Vertical dimension (m)
Nx = 50;    % Number of elements in the x-direction
Ny = 50;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water  = Fluid('water');
gas    = Fluid('gas');
gas.mu = 4.0e-3;  % Viscosity (Pa*s)

% Create porous media
rock = PorousMedia('rock');
rock.K                  = 1.0e-13;             % Intrinsic permeability (m2)
rock.phi                = 0.206;               % Porosity
rock.Pb                 = 1.0e+3;              % Gas-entry pressure
rock.lambda             = 1.0;                 % Curve-fitting parameter
rock.m                  = 2.0;                 % Expoent for the polynomial relationships
rock.liqRelPermeability = 'PolynomialLiquid';  % Liquid relative permeability
rock.gasRelPermeability = 'PolynomialGas';     % Gas relative permeability
rock.capillaryPressure  = 'BrooksCorey';       % Saturation degree function

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Set Dirichlet boundary conditions
mdl.setCapillaryPressureDirichletBCAtPoint([0.0, 0.0], 0.0);
mdl.setCapillaryPressureDirichletBCAtPoint([Lx, Ly], 105018.5554);
mdl.setGasPressureDirichletBCAtPoint([0.0, 0.0], 1.5e6);
mdl.setGasPressureDirichletBCAtPoint([Lx, Ly], 1.0e6);

% Set initial conditions
mdl.setInitialCapillaryPressureAtDomain(105018.5554);
mdl.setInitialGasPressureAtDomain(105018.5554);

%% PROCESS

% Analysis parameters
day       = 60*60*24;   % Conversion from days to seconds
ti        = 0.01*day;   % Initial time
dt        = 0.01*day;   % Time step
tf        = 10.0*day;   % Final time
dtmax     = 1.0*day;    % Maximum time step
dtmin     = 0.001*day;  % Minimum time step
adaptStep = true;       % Adaptive step size

% Run analysis
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.setRelativeConvergenceCriteria(true);
anl.maxIter = 15;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('LiquidSaturation');
mdl.plotField('GasSaturation');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, Ly];
mdl.plotFieldAlongSegment('LiquidPressure', Xi, Xf, 500, 'x');
