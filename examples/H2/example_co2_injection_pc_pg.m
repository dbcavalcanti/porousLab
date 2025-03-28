%% DESCRIPTION
%
% Dense non-liquid phase infiltration problem using the Pc-Pg two-phase 
% flow formulation
%
% References:
%
% Benisch, K., Graupner, B., & Bauer, S. (2013). 
% The coupled OpenGeoSys-eclipse simulator for simulation of CO2
% storage–code comparison for fluid flow and geomechanical processes.
% Energy Procedia, 37, 3663-3671.
%
% Graupner, B. J., Li, D., & Bauer, S. (2011). The coupled simulator 
% ECLIPSE–OpenGeoSys for the simulation of CO2 storage in saline 
% formations. Energy Procedia, 4, 3794-3800.
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
Lx = 200.0;     % Horizontal dimension (m)
Ly = 6.0;       % Vertical dimension (m)
Nx = 30;        % Number of elements in the x-direction
Ny = 60;        % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny,[],[],'ISOQ4',true,false);
mdl.setMesh(node,elem);

mdl.isAxisSymmetric = true;

% --- Material properties of the domain -----------------------------------

% Create the fluids
brine     = Fluid('brine');
brine.rho = 1173.0;
brine.mu  = 1.252e-3;
co2       = Fluid('co2');
co2.rho   = 848.0;
co2.mu    = 8.1e-5;

% Porous media properties
% --------------------------   |   K(m2) | phi | biot |  Ks   | Slr | Sgr |   Pb  | lambda | LiqRelPerm |  GasRelPerm  |  capPressure
aquifer = PorousMedia('aquifer', 3.0e-12 , 0.26 , 1.0 , 1.0e25 , 0.35 , 0.0 , 1.0e4 , 2.0 , 'BrooksCorey', 'BrooksCorey','BrooksCorey');

% Set the material to the model
mdl.setMaterial(aquifer, brine, co2);

% Activate the gravity
mdl.gravityOn = true;

% --- Boundary and initial conditions -------------------------------------

% Capillary pressure conditions
mdl.setCapillaryPressureDirichletBCAtBorder('right', 1.0e4);
mdl.setInitialCapillaryPressureAtDomain(1.0e4);

% Gas pressure conditions
% Hydrostatic profile
depth = 1500.0;          % Depth of the aquifer (top border)
depth0 = depth + Ly;     % Depth of the aquifer (bottom border)
grav   = 9.806;          % Gravity acceleration (m/s2)
for i = 1:mdl.nnodes
    h = depth0 - mdl.NODE(i,2);
    pl = grav * brine.rho * h;
    pc = 1.0e4;
    mdl.setInitialGasPressureAtNode(i,pc+pl);
    % Fix the gas pressure at the right borders
    if ((abs(mdl.NODE(i,1) - Lx))<1.0e-12)
        mdl.setGasPressureDirichletBCAtNode(i, pc+pl);
    end
end

% Add prescribed gas pressure at the infiltration zone
day = 60 * 60 * 24;
qinj = (1.0 / day) * Ly / (Ny + 1);
mdl.setGasPressureNeumannBCAtBorder('left',qinj);

% --- Numerical model configuration ---------------------------------------

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% PROCESS

% Transient analysis parameters
tinit = 0.05*day;   % Initial time
dt    = 0.05*day;   % Time step
tf    = 1*day;      % Final time
dtmax = 0.1*day;    % Maximum time step
dtmin = 1.0e-3*day; % Minimum time step

% Solve the problem
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.setRelativeConvergenceCriteria(true);
anl.run(mdl);

%% POST-PROCESS

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');
