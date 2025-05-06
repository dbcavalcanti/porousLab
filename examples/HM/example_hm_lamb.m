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
[node, elem] = regularMesh(10.0, 16.0, 31, 41);

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
mdl.addLoadAtBorder('top', 2, -10.0e3);

% Pressure
mdl.setPressureDirichletBCAtBorder('top', 0.0);

%% DISCONTINUITY

% Polyline that defines the discontinuity
%                  x     y
fracture_geom = [ 0.0  , 5.0;
                  10.0 , 5.0];

fracture = Discontinuity(fracture_geom, true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.shearStiffness  = 5.0e6;
fracture.normalStiffness = 5.0e6;
fracture.initialAperture = 0.000106265856;
fracture.fluid           = water;
fracture.leakoff         = 0.0;

% Add fractures to model
discontinuityData = struct('addStretchingMode', false, 'addRelRotationMode', false);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

day = 24*60*60;

% Analysis parameters
ti = day;        % Initial time
dt = day;        % Time step
tf = 10.0*day;  % Final time

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf);
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Pressure');
hold on;
fracture.plotIntersectedGeometry();