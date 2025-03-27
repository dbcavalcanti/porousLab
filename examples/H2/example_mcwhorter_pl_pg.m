%% DESCRIPTION
%
% McWhorter and Sunada problem using the Pl-Pg two-phase flow formulation
%
% Reference:
% D.B. McWhorter and D.K. Sunada. Exact integral solutions for
% two-phase flow. Water Resources Research, 26(3):399–413, 1990.
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

% Create model
mdl = Model_H2();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 2.6;      % Horizontal dimension (m)
Ly = 0.5;      % Vertical dimension (m)
Nx = 100;      % Number of elements in the x-direction
Ny = 1;        % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water');
gas   = Fluid('gas');

% Create the porous media
rock = PorousMedia('rock');
rock.K      = 1.0e-10;                       % Intrinsic permeability (m2)
rock.phi    = 0.3;                           % Porosity
rock.Slr    = 0.0;                           % Residual liquid saturation
rocK.Sgr    = 0.0;                           % Residual gas saturation
rock.Pb     = 5.0e3;                         % Gas-entry pressure
rock.lambda = 2.0;                           % Curve-fitting parameter
rock.liqRelPermeability = 'BrooksCorey';
rock.gasRelPermeability = 'BrooksCorey';
rock.capillaryPressure  = 'BrooksCorey';

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'liquidFluid',water,...
    'gasFluid',gas);

% --- Boundary and initial conditions -------------------------------------

% Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('left',195.0e3);
mdl.setGasPressureDirichletBCAtBorder('left',200.0e3);

% Initial conditions
mdl.setInitialPressureAtDomain(-50.0e3);

% --- Numerical model configuration ---------------------------------------

mdl.intOrder = 3;

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% PROCESS

% Transient analysis parameters
tinit = 0.1;        % Initial time
dt    = 0.1;        % Time step
tf    = 1000;       % Final time
dtmax = 10.0;       % Maximum time step
dtmin = 0.0000001;  % Minimun time step

% Solve the problem
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.maxIter = 10;
anl.process(mdl);

%% POST-PROCESS

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('CapillaryPressure');