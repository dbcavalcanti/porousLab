%% DESCRIPTION
%
% Block crossed by strong discontinuity.
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
mdl.isPlaneStress = true;

%% MESH

% Create mesh
Lx = 2.0;  % Horizontal dimension (m)
Ly = 2.0;  % Vertical dimension (m)
Nx = 1;    % Number of elements in the x-direction
Ny = 1;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node,elem);

% Create discontinuities 
Dx = [0.00; 2.00];  % X-coordinates of polyline defining the fracture
Dy = [0.25; 1.75];  % Y-coordinates of polyline defining the fracture
fracture = Discontinuity([Dx, Dy], true);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 1.0e+8;  % Young modulus (kPa)
rock.nu    = 0.0;     % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.initialAperture = 0.0;
fracture.shearStiffness  = 1.0;
fracture.normalStiffness = 1.0;

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);

% Loads
mdl.addLoadAtPoint([0.0,2.0], [-0.5,1.5]);                            

%% PRE-PROCESS

% Create discontinuity elements
fracture.intersectMesh(mdl);

% Add fractures to model
discontinuityData = struct('addStretchingMode', false, 'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot model
mdl.plotField('Model');
hold on;
fracture.plotIntersectedGeometry();
