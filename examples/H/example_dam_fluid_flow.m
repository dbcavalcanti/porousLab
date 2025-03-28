%% DESCRIPTION
%
% Fluid-flow analysis through the foundation of a gravity dam.
%
% Reference:
% Segura, J. M., & Carol, I. (2004).
% On zero-thickness interface elements for diffusion problems.
% Int Journal for Numerical and Analytical Methods in Geomechanics,
% 28(9), 947-962.
% https://doi.org/10.1002/nag.358
%
% Physics:
% * Single-phase hydraulic (H)
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

% Create model
mdl = Model_H();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 24.0;      % Horizontal dimension (m)
Ly = 6.0;       % Vertical dimension (m)
Nx = 96;        % Number of elements in the x-direction
Ny = 24;        % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water');
water.K = 2.0e9;

rock = PorousMedia('rock');
rock.K     = 1.0194e-14;    % Intrinsic permeability (m2)
rock.phi   = 0.3;           % Porosity
rock.Ks    = 1.0e12;        % Rock bulk modulus (Pa)
rock.biot  = 0.6;           % Biot coefficient

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluid',water);

% --- Boundary conditions -------------------------------------------------

for i = 1:size(mdl.NODE,1)
    if ((mdl.NODE(i,1)<8.0) && (abs(mdl.NODE(i,2)-Ly) < 1.0e-9))
        mdl.setPressureDirichletBCAtNode(i,120.0);
    end
    if ((mdl.NODE(i,1)>12.0) && (abs(mdl.NODE(i,2)-Ly) < 1.0e-9))
        mdl.setPressureDirichletBCAtNode(i,60.0);
    end
end

%% RUN ANALYSIS

% Transient analysis parameters
tinit = 1.0;   % Initial time
dt    = 1.0;   % Time step
tf    = 500;  % Final time

% Solve the problem
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(tinit,dt,tf,50.0,0.001,true);
anl.run(mdl);

%% POST-PROCESSING

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , Ly/2.0];
Xf  = [Lx , Ly/2.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('Pressure')

