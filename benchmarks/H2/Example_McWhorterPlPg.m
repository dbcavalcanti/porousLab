%% ================ Two-Phase flow in porous media ========================
% 
% This scripts solves the McWhorter and Sunada (1991) problem
%
% Reference: D.B. McWhorter and D.K. Sunada. Exact integral solutions for
% two-phase flow. Water Resources Research, 26(3):399–413, 1990.
%
% Author: Danilo Cavalcanti
%
%% ========================================================================
%
% Initialize workspace
clear
initWorkspace; 
%
%% ========================== MODEL CREATION ==============================

mdl = Model_H2();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 2.6;      % Horizontal dimension (m)
Ly = 0.5;      % Vertical dimension (m)
Nx = 260;      % Number of elements in the x-direction
Ny = 5;        % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% --- Material properties of the domain -----------------------------------

% Create the fluid components
water = Fluid('water',1000.0,1.0e-6,1.0e25);
gas   = Fluid('gas'  ,1000.0,1.0e-6,1.0e25);

% Porous media properties
% -----------------------  |   K(m2) | phi | biot |  Ks   | Slr | Sgr | Pb | lambda | LiqRelPerm |  GasRelPerm  |  capPressure
rock = PorousMedia('rock'  , 1.0e-10 , 0.3 , 1.0 , 1.0e25 , 0.0 , 0.0 , 5.0 , 2.0 ,'BrooksCorey','BrooksCorey','BrooksCorey');

rock.setMinLiquidRelPermeability(1.0e-5);
rock.setMinGasRelPermeability(1.0e-5);

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'liquidFluid',water,...
    'gasFluid',gas);

% --- Boundary conditions -------------------------------------------------
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Capillary pressure boundary conditions
CoordSupp  = [1 0 -1];                     
CoordLoad  = [];              
CoordPresc = [195.0 0 -1];
CoordInit  = [];
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = -50.0*ones(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
CoordSupp  = [1 0 -1];                % [r cx cy] If cx,cy<0, line              
CoordLoad  = [];                      % [q cx cy] If cx,cy<0, line [m3/s]
CoordPresc = [200.0 0 -1];            % [p cx cy] If cx,cy<0, line [kPa]
CoordInit  = [];                      % [p cx cy] If cx,cy<0, line [kPa]
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 2;

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
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
tinit = 0.1;        % Initial time
dt    = 0.1;        % Time step
tf    = 1000;       % Final time
dtmax = 10.0;       % Maximum time step
dtmin = 0.0000001;  % Minimun time step

% Solve the problem
anl = Anl_Transient0(result,"Picard");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('CapillaryPressure');
