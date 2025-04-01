%% DESCRIPTION
%
% Uniform traction on a plate with linear elastic model and quadratic finite element mesh.
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

%% MESH

% Create mesh
Lx = 0.11;  % Horizontal dimension (m)
Ly = 0.04;  % Vertical dimension (m)
Nx = 22;    % Number of elements in the x-direction
Ny = 8;     % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);
[node, elem] = convertToQuadraticMesh(node, elem);

% Set mesh to model
mdl.setMesh(node, elem);
mdl.resequenceNodes();

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical = 'elastic';  % Elastoplastic with von Mises criteria 
rock.Young      = 2.0e+10;      % Young modulus (Pa)
rock.nu         = 0.0;          % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left', [0.0, 0.0]);

% Loads
mdl.addLoadAtBorder('right', 1, 2.0e+6);

%% PROCESS

anl = Anl_Linear();
anl.run(mdl);

%% POS-PROCESS

% Print results to command window
mdl.printResults();

% Plot contours
mdl.plotField('Ux');
mdl.plotField('Sx');
