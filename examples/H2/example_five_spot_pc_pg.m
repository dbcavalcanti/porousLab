%% DESCRIPTION
%
% Five-spot corner problem using the Pc-Pg two-phase flow formulation
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
mdl = Model_H2_PcPg();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 10.0;      % Horizontal dimension (m)
Ly = 10.0;      % Vertical dimension (m)
Nx = 50;        % Number of elements in the x-direction
Ny = 50;        % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

% --- Material properties of the domain -----------------------------------

% Create the fluids
water  = Fluid('water');
gas    = Fluid('gas');
gas.mu = 4.0e-3;

% Create the porous media
rock = PorousMedia('rock');
rock.K      = 1.0e-13;                       % Intrinsic permeability (m2)
rock.phi    = 0.206;                         % Porosity
rock.Pb     = 1.0e3;                         % Gas-entry pressure
rock.lambda = 1.0;                           % Curve-fitting parameter
rock.m      = 2.0;                           % Expoent for the polynomial relationships
rock.liqRelPermeability = 'PolynomialLiquid';
rock.gasRelPermeability = 'PolynomialGas';
rock.capillaryPressure  = 'BrooksCorey';

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'liquidFluid',water,...
    'gasFluid',gas);

% --- Boundary and initial conditions -------------------------------------

% Dirichlet boundary conditions
mdl.setCapillaryPressureDirichletBCAtPoint([0.0, 0.0], 0.0);
mdl.setCapillaryPressureDirichletBCAtPoint([Lx,   Ly], 105018.5554);
mdl.setGasPressureDirichletBCAtPoint([0.0, 0.0], 1.5e6);
mdl.setGasPressureDirichletBCAtPoint([Lx,   Ly], 1.0e6);

% Initial conditions
mdl.setInitialCapillaryPressureAtDomain(105018.5554);
mdl.setInitialGasPressureAtDomain(105018.5554);

% --- Numerical model configuration ---------------------------------------

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% PROCESS

day = 60*60*24;

% Transient analysis parameters
tinit = 0.01*day;   % Initial time
dt    = 0.01*day;   % Time step
tf    = 10*day;      % Final time
dtmax = 1.0*day;
dtmin = 0.001*day;

% Solve the problem
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.setRelativeConvergenceCriteria(true);
anl.maxIter = 15;
anl.process(mdl);

%% POST-PROCESS

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , Ly];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('LiquidSaturation');
mdl.plotField('GasSaturation');