%% ===================== Elastic plate problem ============================
%
% Elastic traction of a elastic plate validation problem
%
% Author: Danilo Cavalcanti
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

% Mesh properties
Lx = 0.11;     % Horizontal dimension (m)
Ly = 0.04;     % Vertical dimension (m)
Nx = 22;       % Number of elements in the x-direction
Ny = 8;        % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
[node,elem] = convertToQuadraticMesh(node,elem);
mdl.setMesh(node,elem);
mdl.resequenceNodes();

% --- Material properties of the domain -----------------------------------

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'elastic';      % Elastoplastic with von Mises criteria 
rock.Young = 2.0e10;              % Young modulus (Pa)
rock.nu    = 0.0;                 % Poisson ratio

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

% --- Boundary conditions -------------------------------------------------

mdl.setDisplacementDirichletBCAtBorder('left',[0.0, 0.0]);

% Apply pressure at the top (Pa)
mdl.addLoadAtBorder('right', 1, 2.0e6);

%% RUN ANALYSIS

anl = Anl_Linear();
anl.process(mdl);

%% POS-PROCESSING

mdl.printResults();
mdl.plotField('Ux');
mdl.plotField('Sx');
