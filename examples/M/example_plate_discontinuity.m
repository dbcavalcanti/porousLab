%% DESCRIPTION
%
% Block crossed by a strong discontinuity
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

mdl = Model_M();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

% Generate the mesh
[node,elem] = regularMesh(2.0, 2.0, 1, 1);
mdl.setMesh(node,elem);

% Polyline that defines the fracture 
xd = [0.0 , 2.0 ];
yd = [0.25, 1.75];

% Create the discontinuity 
fracture = Discontinuity([xd', yd'],true);

% --- Material properties of the domain -----------------------------------

% Create the porous media
rock = PorousMedia('rock');   
rock.Young = 1.0e8;               % Young modulus (kPa)
rock.nu    = 0.0;                 % Poisson ratio

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

% Set the fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.initialAperture = 0.0;
fracture.shearStiffness  = 1.0;
fracture.normalStiffness = 1.0;

% --- Boundary conditions -------------------------------------------------

mdl.setDisplacementDirichletBCAtBorder('bottom',[0.0, 0.0]);
mdl.addLoadAtPoint([0.0,2.0],[-0.5,1.5]);                            

% --- Additional configurations -------------------------------------------

% Set the problem in a plane stress condition
mdl.isPlaneStress = true;

%% PRE-PROCESSING

% Create the discontinuity elements
fracture.intersectMesh(mdl);

% Add the fracture to the model
discontinuityData = struct( ...
    'addStretchingMode', false,...
    'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture,discontinuityData);

%% RUN ANALYSIS

anl = Anl_Linear();
anl.run(mdl);

%% POS-PROCESSING

% Plot pressure along a segment
mdl.printResults();

mdl.plotField('Model'); hold on
fracture.plotIntersectedGeometry()