%% ===================== Elastic plate problem ============================
%
% Elastic traction of a elastic plate validation problem
%
% Author: Danilo Cavalcanti
%
%% ========================================================================
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
Ny = 8;       % Number of elements in the y-direction

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

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

% --- Boundary conditions -------------------------------------------------

mdl.setDisplacementDirichletBCAtBorder('left',[0.0, 0.0]);

% Apply pressure at the top (Pa)
mdl.addLoadAtBorder('right', 1, 2.0e6);

%% RUN ANALYSIS

% Solve the problem
anl = Anl_Nonlinear('ArcLengthCylControl',true,0.01,200,15,100,4,1.0e-5);
anl.process(mdl);

%% POS-PROCESSING
% mdl.printResults();
mdl.plotField('Ux');
mdl.plotField('Sx');
mdl.plotField('Sy');
mdl.plotField('Sxy');

