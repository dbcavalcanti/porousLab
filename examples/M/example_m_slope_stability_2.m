%% DESCRIPTION
%
% Slope stability problem.
%
% Physics:
% * Mechanical (M)
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
mdl = Model_M();

% Set model options
mdl.intOrder  = 2;     % Integration rule order for the domain
mdl.gravityOn = true;

%%  MESH

% Load nodes and elements
load('NODE');
load('ELEM');

% Set mesh to model
mdl.setMesh(NODE, ELEM);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical    = 'druckerPrager';  % Mechanical constitutive law
rock.rho           = 2.0e+3;           % Density (kg/m3)
rock.Young         = 2.0e+7;           % Young modulus (Pa)
rock.nu            = 0.49;             % Poisson ratio
rock.cohesion      = 5.0e+4;           % Cohesion (Pa)
rock.frictionAngle = 20*pi/180;        % Friction angle (rad)
rock.dilationAngle = 20*pi/180;        % Dilation angle (rad)

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
CoordSupp  = [1 0 0 -1;
              1 0 75.0 -1;
              1 1 -1 0];
CoordLoad  = [];
CoordPresc = [];                                   
           
% Define supports and loads
[mdl.SUPP_u, mdl.LOAD_u, mdl.PRESCDISPL_u] = boundaryConditionsDisplacement(mdl.NODE, CoordSupp, CoordLoad, CoordPresc, 1.0, 1.0, 100, 100);

%% PROCESS

% Analysis parameters
adapt_incr = true;    % Increment size adjustment
increment  = 0.1;     % Initial increment of load ratio
max_lratio = 4.0;     % Limit value of load ratio
max_step   = 11;      % Maximum number of steps
max_iter   = 100;     % Maximum number of iterations in each step
trg_iter   = 4;       % Desired number of iterations in each step
tol        = 1.0e-5;  % Numerical tolerance for convergence

% Create result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

% Run analysis
anl = Anl_Nonlinear(result,'LoadControl', adapt_incr, increment, max_lratio, max_step, max_iter, trg_iter, tol);
anl.run(mdl);

%% POST-PROCESS

% Plot model
mdl.plotField('Model');
mdl.plotDeformedMesh(1.0);

% Plot contours
mdl.plotField('S1');
mdl.plotField('Ux');
mdl.plotField('Sx');
mdl.plotField('Sy');
mdl.plotField('PEMAG');
