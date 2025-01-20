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
mdl.physics = 'H2_PcPg';

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 0.1;     % Horizontal dimension (m)
Ly = 1.0;     % Vertical dimension (m)
Nx = 3;       % Number of elements in the x-direction
Ny = 24;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,1.0e25);
gas = IdealGas('gas', 1.8e-5, 1.0e25);

% Create the porous media
rock = PorousMedia('rock',4.5e-13,0.2975,1.0,1.0e25,0.0,0.2,0.0,3.0,'Liakopoulos','BrooksCorey','Liakopoulos');
rock.setMinLiquidRelPermeability(0.0001);
rock.setMinGasRelPermeability(0.0001);

% Set the user material capillary pressure vs. saturation law (NOT BEING USED)
% --------- Pc  |  Sl
SlPcUMAT = [
    9938.064, 0.900000007;
    9516.018, 0.910000009;
    9065.393, 0.920000003;
    8589.271, 0.929821599;
    8052.432, 0.939999996;
    7469.886, 0.950000006;
    6813.948, 0.960000005;
    6052.562, 0.97;
    5121.663, 0.980000004;
    4847.584, 0.982500005;
    4549.372, 0.984999999;
    4220.252, 0.9875;
    3849.668, 0.989999997;
    3419.508, 0.9925;
    2893.579, 0.995000002;
    2174.942, 0.997499999;
    1000.000, 0.99962099;
       0.000, 1.0];
rock.setUMATCapillaryPressureCurve(SlPcUMAT);

% Set the user material capillary pressure vs. saturation law (NOT BEING USED)
% -- Sl  |  klr
krlUMAT = [
    0.575, 0.00000;
    0.600, 0.12693;
    0.625, 0.18214;
    0.650, 0.23730;
    0.675, 0.29242;
    0.700, 0.34748;
    0.725, 0.40248;
    0.750, 0.45743;
    0.775, 0.51231;
    0.800, 0.56711;
    0.825, 0.62184;
    0.850, 0.67646;
    0.875, 0.73098;
    0.900, 0.78536;
    0.925, 0.83958;
    0.950, 0.89358;
    0.975, 0.94723;
    1.000, 1.00000;
];
rock.setUMATLiquidRelPermCurve(krlUMAT);

% Activate gravity
rock.gravityOn = true;

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'liquidFluid',water,...
    'gasFluid',gas);

% --- Boundary conditions -------------------------------------------------
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Capillary pressure boundary conditions
CoordSupp  = [1 -1 0];                              
CoordLoad  = [];                      
CoordPresc = [];                     
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = 0.0*ones(size(mdl.INITCOND_p,1),1);  % Based on OGS6

% Gas pressure boundary conditions
CoordSupp  = [1 -1 0;
              1 -1 Ly];                        
CoordLoad  = [];                     
CoordPresc = [101325 -1 0;
              101325 -1 Ly];   
CoordInit  = []; 
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_pg = 101325*ones(size(mdl.INITCOND_pg,1),1);

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

% Conversion from days to seconds
minute = 60;

% Transient analysis parameters
tinit = 0.001*minute;          % Initial time
dt    = 0.001*minute;          % Time step
tf    = 5 * minute;          % Final time
dtmax = 0.001*minute;            % Time step
dtmin = 0.001*minute;         % Time step

% Solve the problem
anl = Anl_TransientPicard(result);
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.setPicardRelaxation();
anl.useRelativeError = true;
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
% mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [0.0 , Ly];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotGasPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotCapillaryPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');
