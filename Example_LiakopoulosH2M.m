%% ======================== Liakopoulos problem ===========================
%
% Consists in a drainage test of a sand column.
%
% References:
%
% Schrefler, B. A., and Z. Xiaoyong (1993).
% A fully coupled model for water flow and airflow in deformable
% porous media, Water Resour. Res., 29(1), 155–167,
% doi:10.1029/92WR01737.
%
% Schrefler, B. A., & Scotta, R. (2001).
% A fully coupled dynamic model for two-phase fluid flow in deformable
% porous media. Computer methods in applied mechanics and engineering,
% 190(24-25), 3223-3246.
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

mdl = Model_H2M();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 0.1;     % Horizontal dimension (m)
Ly = 1.0;     % Vertical dimension (m)
Nx = 4;       % Number of elements in the x-direction
Ny = 40;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

%% ============================= MATERIAL =================================

% Create the fluids
water = Fluid('water',1000.0, 1.0e-3, 2.0e9);
gas   = Fluid('gas'  ,1.20  , 1.8e-5, 1.0e5);

% Create the porous media
rock = PorousMedia('rock',4.5e-13,0.2975,1.0,1.0e12,0.2,0.0,0.0,3.0,'Liakopoulos','BrooksCorey','Liakopoulos');
rock.setMechanicalProperties(1.3e6,0.4);
rock.setDensity(2000.0);
rock.setMinLiquidRelPermeability(1.0e-9);
rock.setMinGasRelPermeability(1.0e-9);

% Activate gravity
rock.gravityOn = true;

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'liquidFluid',water,...
    'gasFluid',gas);

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

% Capillary pressure boundary conditions
CoordSupp  = [1 -1 0];                              
CoordLoad  = [];                      
CoordPresc = [101325.0 -1 0];                    
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = 101325.0*ones(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
CoordSupp  = [1 -1  0;
              1  0 -1;
              1  0 -1;
              1 -1 Ly];                        
CoordLoad  = [];                     
CoordPresc = [101325.0 -1  0;
              101325.0  0 -1;
              101325.0  0 -1;
              101325.0 -1 Ly];   
CoordInit  = []; 
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_pg = 101325.0*ones(size(mdl.INITCOND_pg,1),1);

%% ===================== MODEL CONFIGURATION ==============================

% Using Gauss quadrature
mdl.intOrder = 2;

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = false;
mdl.lumpStrategy = 2;

%% ========================= INITIALIZATION ===============================

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();

% Initialize the stress tensor
for el = 1:mdl.nelem
    for i = 1:mdl.element(el).type.nIntPoints
        mdl.element(el).type.intPoint(i).stressOld = [101325.0;101325.0;101325.0;0.0];
    end
end

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Transient analysis parameters
tinit = 5.0;          % Initial time
dt    = 5.0;          % Time step
tf    = 60.0;           % Final time
dtmax = 5.0;          % Time step
dtmin = 5.0;        % Time step

% Solve the problem
% anl = Anl_TransientPicard(result);
% anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
% anl.setPicardRelaxation();
% anl.useRelativeError = true;

anl = Anl_Transient(result);
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
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
