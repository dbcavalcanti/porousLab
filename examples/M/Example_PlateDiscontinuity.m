%% ========================================================================
%
% Block crossed by a strong discontinuity
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
Nx = 1;       % Number of elements in the x-direction
Ny = 1;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Polyline that defines the fracture 
xd = [0.0, 2.0];
yd = [0.25, 1.75];

% Create the discontinuity 
fracture = Discontinuity([xd', yd'],true);

%% ============================= MATERIAL =================================

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'elastic';     
rock.Young = 1.0e8;               % Young modulus (kPa)
rock.nu    = 0.0;                 % Poisson ratio

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

% Set the fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.initialAperture = 0.0;
fracture.shearStiffness  = 1.0;
fracture.normalStiffness = 1.0;

%% ======================= BOUNDARY CONDITIONS ============================

% Displacement boundary conditions
CoordSupp  = [1 1 -1 0.0];
CoordLoad  = [-0.5 1.5 0.0 Ly];
CoordPresc = [];                                   
           
% Define supports and loads
[mdl.SUPP_u, mdl.LOAD_u, mdl.PRESCDISPL_u] = boundaryConditionsDisplacement(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, Lx, Ly, Nx, Ny);

%% ===================== MODEL CONFIGURATION ==============================

% Set the problem in a plane stress condition
mdl.isPlaneStress = true;

% Using Gauss quadrature
mdl.intOrder = 2;

%% ========================= INITIALIZATION ===============================

% Create the discontinuity elements
fracture.intersectMesh(mdl);

% Add the fracture to the model
discontinuityData = struct( ...
    'addStretchingMode', false,...
    'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture,discontinuityData);

% Perform the basic pre-computations associated to the model
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
mdl.printResults();

mdl.plotField('Model'); hold on
fracture.plotIntersectedGeometry()