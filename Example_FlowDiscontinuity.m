%% ========================================================================
%
% Block crossed by a strong discontinuity
%
% Author: Danilo Cavalcanti
%
%% ========================================================================
%
% Initialize workspace
clear
initWorkspace; 
%
%% ============================== MESH  ===================================

mdl = Model_H();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 5.0;     % Horizontal dimension (m)
Ly = 3.0;     % Vertical dimension (m)
Nx = 51;       % Number of elements in the x-direction
Ny = 31;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Polyline that defines the fracture 
xd = [1.0, 4.0];
yd = [1.1, 1.9];

% Create the discontinuity 
fracture = Discontinuity([xd', yd'],true);

%% ============================= MATERIAL =================================

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,2.0e9);

% Create the porous media
rock = PorousMedia('rock');
rock.K     = 9.8e-16;        % Intrinsic permeability (m2)
rock.phi   = 0.25;          % Porosity

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluid',water);

% Set the fracture material properties
fracture.initialAperture = 1.0e-3;
fracture.fluid = water;

%% ======================= BOUNDARY CONDITIONS ============================

% Pore pressure boundary conditions
CoordSupp  = [1 0.0 -1;1 Lx -1];         
CoordLoad  = [];            
CoordPresc = [0.0 0.0 -1;10.0 Lx -1];            
CoordInit  = [];                   
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

%% ===================== MODEL CONFIGURATION ==============================

% Using Gauss quadrature
mdl.intOrder = 2;

%% ========================= INITIALIZATION ===============================

% Create the discontinuity elements
fracture.intersectMesh(mdl);

% Add the fracture to the model
mdl.addPreExistingDiscontinuities(fracture);

% Perform the basic pre-computations associated to the model
mdl.preComputations();

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Transient analysis parameters
tinit = 1.0;   % Initial time
dt    = 1.0;   % Time step
tf    = 400;  % Final time

% Solve the problem
anl = Anl_Transient(result,"Newton");
anl.setUpTransientSolver(tinit,dt,tf,50.0,0.001,true);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
mdl.printResults();

mdl.plotField('Pressure'); hold on
fracture.plotIntersectedGeometry()