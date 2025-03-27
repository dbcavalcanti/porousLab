%% DESCRIPTION
%
% Buckley-Leverett problem using the Pl-Pg two-phase flow formulation
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
Lx = 301.95;   % Horizontal dimension (m)
Ly = 100.0;    % Vertical dimension (m)
Nx = 99;       % Number of elements in the x-direction
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
rock.K      = 1.0e-7;                        % Intrinsic permeability (m2)
rock.phi    = 0.2;                           % Porosity
rock.Slr    = 0.2;                           % Residual liquid saturation
rocK.Sgr    = 0.2;                           % Residual gas saturation
rock.Pb     = 0.0;                           % Gas-entry pressure
rock.lambda = 2.0;                           % Curve-fitting parameter
rock.liqRelPermeability = 'BrooksCorey';
rock.gasRelPermeability = 'BrooksCorey';
rock.capillaryPressure  = 'UMAT';

% Set the user material capillary pressure vs. saturation law
% --------- Pc  |  Sl
SlPcUMAT = [3.0 , 0.2;
            0.0 , 0.8];
rock.setUMATCapillaryPressureCurve(SlPcUMAT);

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'liquidFluid',water,...
    'gasFluid',gas);

% --- Boundary and initial conditions -------------------------------------

% Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('left',200000.0);
mdl.setGasPressureDirichletBCAtBorder('left',200000.230769231);

% Initial conditions
mdl.setInitialPressureAtDomain(200000.0);
mdl.setInitialGasPressureAtDomain(200003.0);

% --- Numerical model configuration ---------------------------------------

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% PROCESS

% Conversion from days to seconds
day = 60*60*24;

% Transient analysis parameters
tinit = 0.1*day;          % Initial time
dt    = 0.1*day;          % Time step
tf    = 50*day;           % Final time
dtmax = 0.1*day;          % Maximum time step
dtmin = 0.1*day;          % Minimum time step

% Solve the problem
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.process(mdl);

%% POST-PROCESS

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotGasPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotCapillaryPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');