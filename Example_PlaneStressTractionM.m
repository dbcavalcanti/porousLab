%% ===================== Elastic plate problem ============================
%
% Elastic traction of a elastic plate validation problem
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

% Mesh properties
Lx = 2;       % Horizontal dimension (m)
Ly = 1.0;     % Vertical dimension (m)
Nx = 20;       % Number of elements in the x-direction
Ny = 10;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

%% ============================= MATERIAL =================================

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'vonMises'; % Elastoplastic with von Mises criteria 
rock.Young = 1.0e6;          % Young modulus (Pa)
rock.nu    = 0.0;             % Poisson ratio
rock.sy0   = 1.0e6;          % Initial yield stress (Pa)
rock.Kp    = 0.0;             % Plastic modulus (Pa)

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

%% ======================= BOUNDARY CONDITIONS ============================
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Displacement boundary conditions
CoordSupp  = [1 0 0 -1;
              1 1 0 0];
CoordLoad  = [];
CoordPresc = [];                                   
           
% Define supports and loads
[mdl.SUPP_u, mdl.LOAD_u, mdl.PRESCDISPL_u] = boundaryConditionsDisplacement(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, Lx, Ly, Nx, Ny);

% Apply pressure at the top (Pa)
[mdl.LOAD_u] = pressureLoad(2.0e6,[Lx, Ly],1,mdl.NODE,mdl.ELEM,mdl.LOAD_u);

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

% Solve the problem
anl = Anl_Nonlinear(result,'ArcLengthCylControl',true,0.01,1.0,10,100,4,1.0e-5);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [0.0 , Ly];
npts = 500;
mdl.plotDeformedMesh(1.0);
mdl.plotField('Ux');
mdl.plotField('Sx');
mdl.plotField('Sy');
mdl.plotField('Sxy');

