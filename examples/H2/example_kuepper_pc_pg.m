%% DESCRIPTION
%
% Kueper and Frind problem using the Pc-Pg two-phase flow formulation
%
% Physics:
% * Two-phase hydraulic (H2)
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
mdl = Model_H2_PcPg();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 0.7;       % Horizontal dimension (m)
Ly = 0.5;       % Vertical dimension (m)
Nx = 56;        % Number of elements in the x-direction
Ny = 40;        % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

% Compute the centroids of the elements
nelem = Nx*Ny;
Xc = zeros(nelem,2);
for el = 1:nelem
    xcoord = mdl.NODE(mdl.ELEM(el,:),1);
    ycoord = mdl.NODE(mdl.ELEM(el,:),2);
    xcentr = sum(xcoord) / 4;            % Considering linear quad elements
    ycentr = sum(ycoord) / 4;            % Considering linear quad elements
    Xc(el,:) = [xcentr , ycentr];
end

% Setting material ids to each element
mdl.matID = ones(nelem,1);

% Sand 2 
reg = isInsideRectangle(Xc,[0.10,0.15],[0.60,0.20]); mdl.matID(reg==1) = 2;

% Sand 3 
reg = isInsideRectangle(Xc,[0.10,0.20],[0.25,0.30]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc,[0.35,0.20],[0.60,0.25]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc,[0.20,0.35],[0.50,0.40]); mdl.matID(reg==1) = 3;

% Sand 4 
reg = isInsideRectangle(Xc,[0.00,0.00],[0.70,0.05]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.05,0.05],[0.20,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.50,0.05],[0.65,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.05,0.15],[0.10,0.40]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.20,0.10],[0.45,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.60,0.15],[0.65,0.40]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.35,0.25],[0.60,0.30]); mdl.matID(reg==1) = 4;

% --- Material properties of the domain -----------------------------------

% Fluid properties
water = Fluid('water');
gas   = Fluid('gas');
gas.rho = 1630.0;
gas.mu  = 0.9e-3;

% Porous media properties
% -----------------------  |   K(m2) | phi | biot |  Ks   | Slr | Sgr |   Pb  | lambda | LiqRelPerm |  GasRelPerm  |  capPressure
sand1 = PorousMedia('sand1', 5.04e-10, 0.40, 1.0,  1.0e25, 0.078, 0.0, 369.73,  3.86,  'BrooksCorey', 'BrooksCorey','BrooksCorey');
sand2 = PorousMedia('sand2', 2.05e-10, 0.39, 1.0,  1.0e25, 0.069, 0.0, 434.45,  3.51,  'BrooksCorey', 'BrooksCorey','BrooksCorey');
sand3 = PorousMedia('sand3', 5.26e-11, 0.39, 1.0,  1.0e25, 0.098, 0.0, 1323.95, 2.49,  'BrooksCorey', 'BrooksCorey','BrooksCorey');
sand4 = PorousMedia('sand4', 8.19e-12, 0.41, 1.0,  1.0e25, 0.189, 0.0, 3246.15, 3.30,  'BrooksCorey', 'BrooksCorey','BrooksCorey');

% Set the material to the model
mdl.setMaterial([sand1, sand2, sand3, sand4], water, gas);

% Activate the gravity
mdl.gravityOn = true;

% --- Boundary and initial conditions -------------------------------------

% Prescribed pressures at the top corners
mdl.setCapillaryPressureDirichletBCAtPoint([0.0, Ly],369.73);
mdl.setCapillaryPressureDirichletBCAtPoint([Lx , Ly],369.73);
mdl.setGasPressureDirichletBCAtPoint([0.0, Ly],369.73);
mdl.setGasPressureDirichletBCAtPoint([Lx , Ly],369.73);

% Add prescribed gas pressure at the infiltration zone
tol = 1.0e-4;
reg = find(isInsideRectangle(mdl.NODE,[0.3-tol,0.5-tol],[0.4+tol,0.5+tol]));
for i = 1:length(reg)
    mdl.setGasPressureDirichletBCAtNode(reg(i), 639.35);
end

% Initial conditions
mdl.setInitialCapillaryPressureAtDomain(369.73);
mdl.setInitialGasPressureAtDomain(369.73);

% --- Order of the integration rule for the domain ------------------------

mdl.intOrder = 2;

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% PROCESS

% Transient analysis parameters
tinit = 1.0;   % Initial time
dt    = 1.0;   % Time step
tf    = 100;    % Final time

% Solve the problem
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(tinit,dt,tf,1.0,0.001,true);
anl.setRelativeConvergenceCriteria(true);
anl.run(mdl);

%% POST-PROCESS

mdl.plotField('GasSaturation',[0.0, 1.0]);