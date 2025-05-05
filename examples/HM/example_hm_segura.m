%% DESCRIPTION
%
% Segura consolidation problem.
%
% Physics:
% * Single-phase flow hydro-mechanical (HM)
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
mdl = Model_HM();

%% MESH

% Create mesh
[node, elem] = regularMesh(1.0, 1.0, 1, 1);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');

% Create porous media
rock = PorousMedia('rock');
rock.K     = 1.15741e-12;   % Intrinsic permeability (m2)
rock.phi   = 0.3;           % Porosity
rock.Young = 1.0e+6;        % Young modulus (Pa)
rock.nu    = 0.3;           % Poisson ratio

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

% Loads
mdl.addLoadAtBorder('top', 2, -1.0e4);

% Pressure
mdl.setPressureDirichletBCAtBorder('top', 0.0);

%% DISCONTINUITY

% Polyline that defines the discontinuity
%                  x     y
fracture_geom = [ 0.5 , 0.0;
                  0.5 , 1.0];

fracture = Discontinuity(fracture_geom, true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.shearStiffness  = 1.0;
fracture.normalStiffness = 1.0;
fracture.initialAperture = 1.0e-3;
fracture.fluid           = water;
fracture.leakoff         = 1.0e-19;

% Add fractures to model
discontinuityData = struct('addStretchingMode', false, 'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

% Analysis parameters
ti = 1.0;    % Initial time
dt = 1.0;    % Time step
tf = 100.0;  % Final time

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf);
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Pressure');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [0.0, 1.0];
mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'y');
mdl.plotFieldAlongSegment('Uy', Xi, Xf, 500, 'y');
