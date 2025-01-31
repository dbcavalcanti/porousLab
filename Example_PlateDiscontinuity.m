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
Lx = 2.0;     % Horizontal dimension (m)
Ly = 2.0;     % Vertical dimension (m)
Nx = 20;       % Number of elements in the x-direction
Ny = 20;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

[mdl.NODE, mdl.ELEM] = convertToQuadraticMesh(mdl.NODE, mdl.ELEM);

mdl.resequenceNodes();

mdl.type = 'ISOQ8';

xd = [0.0, 2.0];
yd = [0.5, 1.5];

% xd = linspace(0, 2, 100);
% yd = 1.0 + 0.5 * sin(0.5 * pi * xd);
% 
% fracture = Discontinuity([xd', yd'],true);
% 
% fracture.setRepelTol(0.1);
% fracture.setSavePerturbNodes(true);
% 
% % Perform intersection and repel process
% fracture.intersectMesh(mdl);
% 
% % Add the fracture to the model
% mdl.addPreExistingDiscontinuities(fracture);

%% ============================= MATERIAL =================================

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'elastic';      % Elastoplastic with von Mises criteria 
rock.Young = 2.0e10;              % Young modulus (Pa)
rock.nu    = 0.0;                 % Poisson ratio

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

%% ======================= BOUNDARY CONDITIONS ============================
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Displacement boundary conditions
CoordSupp  = [1 1 0 -1];
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
anl = Anl_Linear(result);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [0.0 , Ly];
npts = 500;
mdl.plotDeformedMesh(1.0);
mdl.plotField('Ux');
mdl.plotField('Sx');

% mdl.plotField('Model'); hold on
% fracture.plotOriginalGeometry()
% fracture.plotIntersectedGeometry()
% fracture.plotPerturbNodes()