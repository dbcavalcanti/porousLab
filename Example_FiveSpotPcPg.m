%% ================ Two-Phase flow in porous media ====================
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

mdl = Model();

% --- Physics -------------------------------------------------------------
mdl.physics = 'hydraulicTwoPhasePcPg';

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 10.0;      % Horizontal dimension (m)
Ly = 10.0;      % Vertical dimension (m)
Nx = 50;        % Number of elements in the x-direction
Ny = 50;        % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,1.0e25);
gas   = Fluid('gas'  ,1000.0,4.0e-3,1.0e25);

rock = PorousMedia('rock',1.0e-13,0.206,1.0,1.0e25,0.0,0.0,1.0e3,1.0,'PolynomialLiquid','PolynomialGas','BrooksCorey');

rock.m = 2;
rock.setMinLiquidRelPermeability(1.0e-9);
rock.setMinGasRelPermeability(1.0e-9);

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'liquidFluid',water,...
    'gasFluid',gas);

% --- Boundary conditions -------------------------------------------------
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Initial pressure
pi = 105018.5554;
% pi = 1.0e3;

% Capillary pressure boundary conditions
CoordSupp  = [1 0 0;
              1 Lx Ly];                            
CoordLoad  = [];                      
CoordPresc = [pi Lx Ly];            
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = pi*ones(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
CoordSupp  = [1 0 0;
              1 Lx Ly];                             
CoordLoad  = [];                      
CoordPresc = [1.5e6 0 0;
              1.0e6 Lx Ly];             
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_pg = pi*ones(size(mdl.INITCOND_pg,1),1);

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 2;

% Diagonalize compressibility matrix
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% ========================= INITIALIZATION ===============================

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();

% Plot the mesh with the supports
% mdl.plotMeshWithBC();

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

day = 60*60*24;

% Transient analysis parameters
tinit = 0.01*day;   % Initial time
dt    = 0.01*day;   % Time step
tf    = 10*day;      % Final time
dtmax = 1.0*day;
dtmin = 0.001*day;

% Solve the problem
anl = Anl_TransientPicard(result);
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.setPicardRelaxation();
anl.useRelativeError = true;
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , Ly];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
% mdl.plotField('CapillaryPressure');
% mdl.plotField('GasPressure');
mdl.plotField('LiquidSaturation');
mdl.plotField('GasSaturation');