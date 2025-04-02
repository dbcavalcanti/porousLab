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
mdl.gravityOn = true;

%% MESH

% Load nodes and elements
load('MeshSlopeStability');

% Set mesh to model
mdl.setMesh(node,elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical    = 'druckerPrager';  % Mechanical constitutive law
rock.rho           = 2.0;            % Density (kg/m3)
rock.Young         = 2.0e+4;           % Young modulus (Pa)
rock.nu            = 0.49;             % Poisson ratio
rock.cohesion      = 50.0;           % Cohesion (Pa)
rock.frictionAngle = 20*pi/180;        % Friction angle (rad)
rock.dilationAngle = 20*pi/180; 

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);

%% PROCESS

% Configure analysis
anl = Anl_Nonlinear();
anl.method        = 'ArcLengthCylControl';
anl.adjustStep    = true;
anl.increment     = 0.1;
anl.max_increment = 0.1;
anl.max_lratio    = 3.0;
anl.max_step      = 100;
anl.max_iter      = 100;
anl.trg_iter      = 9;

% Node and DOF used to plot Load Factor vs Displacement
ndId = mdl.closestNodeToPoint([35.0, 40.0]);
anl.setPlotDof(ndId, 2);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('E1');
mdl.plotField('S1');
mdl.plotField('PEMAG');
