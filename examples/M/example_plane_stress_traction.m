%% DESCRIPTION
%
% Uniform traction on a plate with an isotropic damage model
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

% Mesh properties
Lx = 0.11;     % Horizontal dimension (m)
Ly = 0.04;     % Vertical dimension (m)
Nx = 22;       % Number of elements in the x-direction
Ny = 8;        % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

% --- Material properties of the domain -----------------------------------

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'isoDamage';      % Elastoplastic with von Mises criteria 
rock.Young = 2.0e10;                % Young modulus (Pa)
rock.nu    = 0.0;                   % Poisson ratio
rock.kappa = 10.0;                  % Ratio of tensile and compressive strength
rock.DamageThreshold = 1.0e-4;
rock.FractureEnergyMode1 = 50.0;

% Set the material to the model
mdl.setMaterial(rock);

% --- Boundary conditions -------------------------------------------------

mdl.setDisplacementDirichletBCAtBorder('left',[0.0, 0.0]);

% Apply pressure at the top (Pa)
mdl.addLoadAtBorder('right', 1, 2.0e6);

%% RUN ANALYSIS

% Set the analysis
anl = Anl_Nonlinear('ArcLengthCylControl',true,0.01,200,15,100,4,1.0e-5);

% Define the node and dof that will be used to plot the load factor vs.
% displacement curve
ndId = mdl.closestNodeToPoint([Lx,0.0]);
anl.setPlotDof(ndId,1)

% Run the analysis
anl.run(mdl);

%% POS-PROCESSING
mdl.plotField('Ux');
mdl.plotField('Sx');
mdl.plotField('Sy');
mdl.plotField('Sxy');

