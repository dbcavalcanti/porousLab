%% ================ Two-Phase flow in porous media ====================
%
% Reference:
% Zhou, Q., J. T. Birkholzer, C.-F. Tsang, and J. Rutqvist (2008), 
% A method for quick assessment of CO 2 storage capacity in closed
% and semi-closed saline formations. Int. J. Greenh. Gas Control, 2(4),
% 626–639, doi:10.1016/j.ijggc.2008.02.004
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
Lx = 20.0e3;      % Horizontal dimension (m)
Ly = 30.0;       % Vertical dimension (m)
Nx = 200;      % Number of elements in the x-direction
Ny = 20;        % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,1.0e25);
gas   = Fluid('gas'  ,1000.0,1.0e-3,1.0e25);

% Parameters of the reservoir
permeability     = 20 * 9.869233 * 1.0e-16; 
porosity         = 0.15;
biotCoeff        = 1.0;
compressibility  = 1.0e25;
resLiqSaturation = 0.0;
resGasSaturation = 0.0;
gasEntryPressure = 1.0;
lambda           = 1.0;

rock = PorousMedia('rock',permeability,porosity,biotCoeff, ...
    compressibility,resLiqSaturation,resGasSaturation,gasEntryPressure, ...
    lambda,'BrooksCorey','BrooksCorey','BrooksCorey');
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

% Capillary pressure boundary conditions
CoordSupp  = [1 0 -1];                % [r cx cy] If cx,cy<0, line              
CoordLoad  = [];                      % [q cx cy] If cx,cy<0, line [m3/s]
CoordPresc = [5.0e3 0 -1];            % [p cx cy] If cx,cy<0, line [kPa]
CoordInit  = [];                      % [p cx cy] If cx,cy<0, line [kPa]
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = 50.0e3*ones(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
CoordSupp  = [1 0 -1];                % [r cx cy] If cx,cy<0, line              
CoordLoad  = [];                      % [q cx cy] If cx,cy<0, line [m3/s]
CoordPresc = [200.0e3 0 -1];          % [p cx cy] If cx,cy<0, line [kPa]
CoordInit  = [];                      % [p cx cy] If cx,cy<0, line [kPa]
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 3;

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

% Transient analysis parameters
tinit = 0.1;   % Initial time
dt    = 0.1;   % Time step
tf    = 50;  % Final time

% Solve the problem
anl = Anl_TransientPicard(result);
anl.setUpTransientSolver(tinit,dt,tf,0.1,0.0000001,true);
% anl.setPicardRelaxation();
anl.useRelativeError = false;
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')

