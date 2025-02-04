%% ====================== Slope stability problem =========================
%
% Elastoplastic pressurized cylinder example.
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

mdl = Model_M();

% --- Mesh of continuum elements ------------------------------------------
load('NODE'); load('ELEM');
mdl.NODE = NODE;
mdl.ELEM = ELEM;

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

%% ============================= MATERIAL =================================

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'elastic'; % Elastoplastic with von Mises criteria 
rock.Young = 2.0e7;           % Young modulus (Pa)
rock.nu    = 0.49;             % Poisson ratio
rock.sy0   = 2.40e8;          % Initial yield stress (Pa)
rock.Kp    = 0.0;             % Plastic modulus (Pa)
rock.rho   = 2039.567612;

rock.gravityOn = true;

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

%% ======================= BOUNDARY CONDITIONS ============================
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Displacement boundary conditions
CoordSupp  = [1 0 0 -1;
              1 0 75.0 -1;
              1 1 -1 0];
CoordLoad  = [];
CoordPresc = [];                                   
           
% Define supports and loads
[mdl.SUPP_u, mdl.LOAD_u, mdl.PRESCDISPL_u] = boundaryConditionsDisplacement(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, 1.0, 1.0, 100, 100);

%% ===================== MODEL CONFIGURATION ==============================

% Using Gauss quadrature
mdl.intOrder = 2;

%% ========================= INITIALIZATION ===============================

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();
mdl.plotField('Model');

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Solve the problem
anl = Anl_Linear(result);
% anl = Anl_Nonlinear(result,'ArcLengthCylControl',true,0.01,2.0,100,100,4,1.0e-5);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
Xi  = [0.0 , min(mdl.NODE(:,2))];
Xf  = [0.0 , max(mdl.NODE(:,2))];
npts = 500;
% mdl.plotDeformedMesh(1.0);
% mdl.plotField('S1');
% mdl.plotField('Ux');
% mdl.plotField('Sx');
% mdl.plotField('Sy');
mdl.plotField('Sxy');

