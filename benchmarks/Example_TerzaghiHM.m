%% ================ Terzaghi consolidation problem ========================
%
% Hydromechanical with single-phase flow validation problem
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

mdl = Model_HM();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 0.1;     % Horizontal dimension (m)
Ly = 1.0;     % Vertical dimension (m)
Nx = 5;       % Number of elements in the x-direction
Ny = 50;      % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

%% ============================= MATERIAL =================================

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,1.0e25);

% Create the porous media
rock = PorousMedia('rock');
rock.K     = 1.15741e-12;   % Intrinsic permeability (m2)
rock.phi   = 0.3;           % Porosity
rock.Young = 1.0e6;         % Young modulus (Pa)
rock.nu    = 0.3;           % Poisson ratio

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluid',water);

%% ======================= BOUNDARY CONDITIONS ============================
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Displacement boundary conditions
CoordSupp  = [1 0 0 -1;
              1 0 Lx -1
              1 1 -1 0.0;];
CoordLoad  = [];
CoordPresc = [];                                   
           
% Define supports and loads
[mdl.SUPP_u, mdl.LOAD_u, mdl.PRESCDISPL_u] = boundaryConditionsDisplacement(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, Lx, Ly, Nx, Ny);

% Apply pressure at the top (Pa)
[mdl.LOAD_u] = pressureLoad(-1.0e4,[Lx, Ly],2,mdl.NODE,mdl.ELEM,mdl.LOAD_u);

% Liquid pressure boundary conditions
CoordSupp  = [1 -1 Ly];                              
CoordLoad  = [];                      
CoordPresc = [];                    
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

%% ===================== MODEL CONFIGURATION ==============================

% Using Gauss quadrature
mdl.intOrder = 2;

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
tinit = 1.0;          % Initial time
dt    = 1.0;          % Time step
tf    = 100;          % Final time
dtmax = 1.0;          % Maximum time step
dtmin = 0.001;        % Minimum time step

% Solve the problem
anl = Anl_Transient(result);
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [0.0 , Ly];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotField('Pressure');

