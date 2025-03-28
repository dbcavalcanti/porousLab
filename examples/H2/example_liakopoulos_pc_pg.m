%% DESCRIPTION
%
% OGS-5 Liakopoulos problem using the Pc-Pg two-phase flow formulation
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
Lx = 0.1;     % Horizontal dimension (m)
Ly = 1.0;     % Vertical dimension (m)
Nx = 3;       % Number of elements in the x-direction
Ny = 24;      % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water');
gas = IdealGas('gas', 1.8e-5, 1.0e25);

% Create the porous media
rock = PorousMedia('rock',4.5e-13,0.2975,1.0,1.0e25,0.0,0.2,0.0,3.0,'Liakopoulos','BrooksCorey','Liakopoulos');
rock.setMinLiquidRelPermeability(0.0001);
rock.setMinGasRelPermeability(0.0001);

% Activate gravity
rock.gravityOn = true;

% Set the material to the model
mdl.setMaterial(rock, water, gas);

% --- Boundary and initial conditions -------------------------------------

% Dirichlet boundary conditions
mdl.setCapillaryPressureDirichletBCAtBorder('bottom',0.0);
mdl.setGasPressureDirichletBCAtBorder('bottom',101325.0);
mdl.setGasPressureDirichletBCAtBorder('top',101325.0);

% Initial conditions
mdl.setInitialCapillaryPressureAtDomain(0.0);
mdl.setInitialGasPressureAtDomain(101325.0);

% --- Numerical model configuration ---------------------------------------

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% PROCESS

% Conversion from days to seconds
minute = 60;

% Transient analysis parameters
tinit = 0.001*minute;          % Initial time
dt    = 0.001*minute;          % Time step
tf    = 5 * minute;          % Final time
dtmax = 1.0*minute;            % Time step
dtmin = 0.001*minute;         % Time step

% Solve the problem
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.run(mdl);

%% POST-PROCESS

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [0.0 , Ly];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotGasPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotCapillaryPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');
