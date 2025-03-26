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
Lx = 0.11;     % Horizontal dimension (m)
Ly = 0.04;     % Vertical dimension (m)
Nx = 22;       % Number of elements in the x-direction
Ny = 8;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE, mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);
[mdl.NODE, mdl.ELEM] = convertToQuadraticMesh(mdl.NODE, mdl.ELEM);

mdl.resequenceNodes();

% Type of elements
mdl.type = 'ISOQ8';

% Thickness (m)
mdl.t = 1.0;

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

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Solve the problem
anl = Anl_Linear(result);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

mdl.plotField('Ux');
mdl.plotField('Sx');
