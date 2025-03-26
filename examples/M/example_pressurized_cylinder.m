%% ================= Pressurized cylinder problem =========================
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

% Mesh properties
ri = 0.1;      % Internal radius of the cylinder
re = 0.2;      % External radius of the cylinder
Nx = 10;       % Number of elements in the x-direction
Ny = 10;       % Number of elements in the y-direction

% Generate the mesh
[NODE,mdl.ELEM] = regularMeshY(1.0, 1.0, Nx, Ny);

% Transform to cylindrical coordinates
r = ri + NODE(:, 1) * (re - ri);
theta = NODE(:, 2) * (pi / 2);

% Set the nodes with the transformed coordinates
mdl.NODE = [r .* cos(theta), r .* sin(theta)];

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

%% ============================= MATERIAL =================================

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'vonMises'; % Elastoplastic with von Mises criteria 
rock.Young = 2.1e11;          % Young modulus (Pa)
rock.nu    = 0.3;             % Poisson ratio
rock.sy0   = 2.40e8;          % Initial yield stress (Pa)
rock.Kp    = 0.0;             % Plastic modulus (Pa)

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

%% ======================= BOUNDARY CONDITIONS ============================
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Displacement boundary conditions
CoordSupp  = [1 0 0 -1;
              0 1 -1 0];
CoordLoad  = [];
CoordPresc = [];                                   
           
% Define supports and loads
[mdl.SUPP_u, mdl.LOAD_u, mdl.PRESCDISPL_u] = boundaryConditionsDisplacement(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, 1.0, 1.0, Nx, Ny);

% Internal pressure (Pa)
pint = 192.0905814164710e06;
% pint = 500;

% Compute the radius of each node wrt to the center at (0,0)
r  = sqrt((mdl.NODE(:,1)).^2 + (mdl.NODE(:,2)).^2);
sn = mdl.NODE(:,2) ./ r;
cs = mdl.NODE(:,1) ./ r;

% Loaded nodes
internalNodes = (r-ri)<eps;

% Number of loaded nodes
nInternalNodes = sum(internalNodes);

% Force magnitude
F0 = (pint*mdl.t*ri*pi/2.0)/(nInternalNodes-1)/2.0;

% Count occurrences of each node
nodeCount = histcounts(mdl.ELEM(:), 1:(size(mdl.NODE,1)+1))';

% Apply the internal pressure to the nodes located at the internal face
mdl.LOAD_u(internalNodes,:) = F0 * nodeCount(internalNodes) .* [cs(internalNodes) , sn(internalNodes)];

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
% anl = Anl_Linear(result);
anl = Anl_Nonlinear(result,'ArcLengthCylControl',true,0.01,2.0,100,100,4,1.0e-5);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
Xi  = [0.0 , min(mdl.NODE(:,2))];
Xf  = [0.0 , max(mdl.NODE(:,2))];
npts = 500;
% mdl.plotDeformedMesh(1.0);
mdl.plotField('S1');
mdl.plotField('Sr');
% mdl.plotField('Ux');
mdl.plotField('Sx');
% mdl.plotField('Sy');
% mdl.plotField('Sxy');

