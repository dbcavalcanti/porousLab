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
Lx = 2.6;      % Horizontal dimension (m)
Ly = 0.5;      % Vertical dimension (m)
Nx = 26;      % Number of elements in the x-direction
Ny = 1;        % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

% --- Material properties of the domain -----------------------------------

% Gas compressibility
R = 8.3144598;      % Universal gas constant (J/(mol*K)
T = 293.15;         % Temperature (K)
Mg = 0.0289;        % Molar mass (kg/mol)
rhog = 1000.0;      % Density (kg/m3)
Kg = (rhog * T * R) /  Mg;

fluids = [Fluid('water',1000.0,1.0e-3,1.0e25),...
          Fluid('CO2',  rhog,1.0e-3,Kg)];

rock = PorousMedia('rock',1.0e-10,0.3,1.0,1.0e25,0.0,5.0e3,2.0,'BrooksCorey','BrooksCorey');
rock.setMinLiquidRelPermeability(1.0e-9);
rock.setMinGasRelPermeability(1.0e-9);

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluids',fluids);

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
mdl.intOrder = 2;

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
tf    = 100;  % Final time

% Solve the problem
anl = Anl_TransientPicard(result);
anl.setUpTransientSolver(tinit,dt,tf,5.0,0.0000001,true);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
