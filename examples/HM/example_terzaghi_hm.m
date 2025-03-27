%% DESCRIPTION
%
% Terzaghi consolidation problem.
%
% Physics:
% * Hydromechanical with single-phase flow (HM)
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

mdl = Model_HM();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

[node,elem] = regularMesh(0.1, 1.0, 5, 50);
mdl.setMesh(node,elem);

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water');

% Create the porous media
rock = PorousMedia('rock');
rock.K     = 1.15741e-12;   % Intrinsic permeability (m2)
rock.phi   = 0.3;           % Porosity
rock.Young = 1.0e6;         % Young modulus (Pa)
rock.nu    = 0.3;           % Poisson ratio

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluid',water);

% --- Boundary conditions -------------------------------------------------

% Displacement boundary conditions
mdl.setDisplacementDirichletBCAtBorder('left',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right', [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom',[0.0, 0.0]);

% Apply pressure at the top (Pa)
mdl.addLoadAtBorder('top', 2, -1.0e4);

% Liquid pressure boundary conditions
mdl.setPressureDirichletBCAtBorder('top',0.0);

%% RUN ANALYSIS

% Transient analysis parameters
tinit = 1.0;          % Initial time
dt    = 1.0;          % Time step
tf    = 100;          % Final time

% Solve the problem
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(tinit,dt,tf);
anl.process(mdl);

%% POST-PROCESSING

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [0.0 , 1.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotField('Pressure');