%% ================= Single-Phase flow in porous media ====================
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

mdl = Model_H();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 24.0;      % Horizontal dimension (m)
Ly = 6.0;       % Vertical dimension (m)
Nx = 96;        % Number of elements in the x-direction
Ny = 24;        % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,2.0e9);


rock = PorousMedia('rock');
rock.K     = 1.0194e-14;    % Intrinsic permeability (m2)
rock.phi   = 0.3;           % Porosity
rock.Ks    = 1.0e12;        % Rock bulk modulus (Pa)
rock.biot  = 0.6;           % Biot coefficient

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluid',water);

% --- Boundary conditions -------------------------------------------------
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% pore pressure boundary conditions
CoordSupp  = [];         
CoordLoad  = [];            
CoordPresc = [];            
CoordInit  = [];                   
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

% Apply hydraulic pressure head at the left
for i = 1:size(mdl.NODE,1)
    if ((mdl.NODE(i,1)<8.0) && (abs(mdl.NODE(i,2)-Ly) < 1.0e-9))
        mdl.PRESCDISPL_p(i) = 120.0;
        mdl.SUPP_p(i) = 1;
    end

    if ((mdl.NODE(i,1)>12.0) && (abs(mdl.NODE(i,2)-Ly) < 1.0e-9))
        mdl.PRESCDISPL_p(i) = 60.0;
        mdl.SUPP_p(i) = 1;
    end
end

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

% Solve the problem
anl = Anl_Transient0(result,"Newton");
anl.setUpTransientSolver(tinit,dt,tf,50.0,0.001,true);
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

