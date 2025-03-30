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
load('NODE');
load('ELEM');

% Set mesh to model
mdl.setMesh(NODE, ELEM);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.Young = 2.0e+7;       % Young modulus (Pa)
rock.nu    = 0.49;         % Poisson ratio
rock.rho   = 2039.567612;  % Density (kg/m3)

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Sxy');
