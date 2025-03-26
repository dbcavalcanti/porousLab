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

mdl = Model_H2_PcPg();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 301.95;   % Horizontal dimension (m)
Ly = 100.0;    % Vertical dimension (m)
Nx = 99;       % Number of elements in the x-direction
Ny = 1;        % Number of elements in the y-direction

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

% Create the porous media
rock = PorousMedia('rock',1.0e-7,0.2,1.0,1.0e25,0.2,0.2,0.0,2.0,'BrooksCorey','BrooksCorey','UMAT');
rock.setMinLiquidRelPermeability(1.0e-9);
rock.setMinGasRelPermeability(1.0e-9);

% Set the user material capillary pressure vs. saturation law
% --------- Pc  |  Sl
SlPcUMAT = [3.0 , 0.2;
            0.0 , 0.8];
rock.setUMATCapillaryPressureCurve(SlPcUMAT);

% Material parameters vector
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
CoordPresc = [0.230769231 0 -1];      % [p cx cy] If cx,cy<0, line [kPa]
CoordInit  = [];                      % [p cx cy] If cx,cy<0, line [kPa]
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = 3.0*ones(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
CoordSupp  = [1 0 -1];                % [r cx cy] If cx,cy<0, line              
CoordLoad  = [];                      % [q cx cy] If cx,cy<0, line [m3/s]
CoordPresc = [200000.230769231 0 -1]; % [p cx cy] If cx,cy<0, line [kPa]
CoordInit  = [];                      % [p cx cy] If cx,cy<0, line [kPa]
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_pg = 200003*ones(size(mdl.INITCOND_pg,1),1);

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 2;

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% ========================= INITIALIZATION ===============================

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Conversion from days to seconds
day = 60*60*24;

% Transient analysis parameters
tinit = 0.1*day;          % Initial time
dt    = 0.1*day;          % Time step
tf    = 50*day;           % Final time
dtmax = 0.1*day;          % Maximum time step
dtmin = 0.1*day;          % Minimum time step

% Solve the problem
anl = Anl_Transient(result,"Picard");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
% mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotGasPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotCapillaryPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');
