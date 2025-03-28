%% DESCRIPTION
%
% Slope stability problem
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

%% MODEL CREATION

mdl = Model_M();

% --- Mesh of continuum elements ------------------------------------------
load('NODE'); load('ELEM');
mdl.setMesh(NODE,ELEM);

% --- Material properties of the domain -----------------------------------

% Create the porous media
rock = PorousMedia('rock');
rock.Young = 2.0e7;          % Young modulus (Pa)
rock.nu    = 0.49;           % Poisson ratio
rock.rho   = 2039.567612;
rock.gravityOn = true;

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

% --- Boundary conditions -------------------------------------------------

mdl.setDisplacementDirichletBCAtBorder('left',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right', [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom',[0.0, 0.0]);

%% RUN ANALYSIS

anl = Anl_Linear();
anl.run(mdl);

%% POS-PROCESSING

mdl.plotField('Sxy');

