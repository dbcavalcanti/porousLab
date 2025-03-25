%% DESCRIPTION
%
% ...
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

%% MESH GENERATION

% Mesh generation
Lx = 2.0;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 2;    % Number of elements in the x-direction
Ny = 1;    % Number of elements in the y-direction

[mdl.NODE, mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Element type
mdl.type = 'ISOQ4';

% Element thickness (m)
mdl.t = 1.0;

%% MATERIAL CREATION

% Create fluids
water = Fluid('water');
water.rho = 1000.0;  % Density (kg/m3)
water.mu  = 1.0e-3;  % Viscosity (Pa*s)
water.K   = 2.0e9;   % Compressibility/Bulk modulus (1/Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K    = 1.0194e-14;  % Intrinsic permeability (m2)
rock.phi  = 0.3;         % Porosity
rock.Ks   = 1.0e12;      % Rock bulk modulus (Pa)
rock.biot = 0.6;         % Biot coefficient

% Material parameters vector
mdl.mat = struct('porousMedia',rock,'fluid',water);

%% BOUNDARY CONDITIONS
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% pore pressure boundary conditions
CoordSupp  = [1 0.0 -1;1 Lx -1];         
CoordLoad  = [];            
CoordPresc = [0.0 0.0 -1;10.0 Lx -1];            
CoordInit  = [];                   
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 2;

% Diagonalize compressibility matrix
mdl.massLumping = false;
mdl.lumpStrategy = 2;

%% ========================= INITIALIZATION ===============================

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Transient analysis parameters
tinit = 1.0;   % Initial time
dt    = 1.0;   % Time step
tf    = 500;  % Final time
dtmax = 1.0;
dtmin = 1.0;

% Solve the problem
anl = Anl_Transient(result,"Picard");
anl.setUpTransientSolver(tinit,dt,tf,50.0,0.001,true);
anl.setRelativeConvergenceCriteria(true);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , Ly/2.0];
Xf  = [Lx , Ly/2.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('Pressure')

