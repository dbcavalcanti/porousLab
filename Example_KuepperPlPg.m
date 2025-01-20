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
mdl.physics = 'H2';

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 0.7;       % Horizontal dimension (m)
Ly = 0.5;       % Vertical dimension (m)
Nx = 56;        % Number of elements in the x-direction
Ny = 40;        % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

% Thickness (m)
mdl.t = 1.0;

% Compute the centroids of the elements
nelem = Nx*Ny;
Xc = zeros(nelem,2);
for el = 1:nelem
    xcoord = mdl.NODE(mdl.ELEM(el,:),1);
    ycoord = mdl.NODE(mdl.ELEM(el,:),2);
    xcentr = sum(xcoord) / 4;            % Considering linear quad elements
    ycentr = sum(ycoord) / 4;            % Considering linear quad elements
    Xc(el,:) = [xcentr , ycentr];
end

% Setting material ids to each element
mdl.matID = ones(nelem,1);

% Sand 2 
reg = isInsideRectangle(Xc,[0.10,0.15],[0.60,0.20]); mdl.matID(reg==1) = 2;

% Sand 3 
reg = isInsideRectangle(Xc,[0.10,0.20],[0.25,0.30]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc,[0.35,0.20],[0.60,0.25]); mdl.matID(reg==1) = 3;
reg = isInsideRectangle(Xc,[0.20,0.35],[0.50,0.40]); mdl.matID(reg==1) = 3;

% Sand 4 
reg = isInsideRectangle(Xc,[0.00,0.00],[0.70,0.05]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.05,0.05],[0.20,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.50,0.05],[0.65,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.05,0.15],[0.10,0.40]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.20,0.10],[0.45,0.15]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.60,0.15],[0.65,0.40]); mdl.matID(reg==1) = 4;
reg = isInsideRectangle(Xc,[0.35,0.25],[0.60,0.30]); mdl.matID(reg==1) = 4;

% --- Material properties of the domain -----------------------------------

% Fluid properties
fluids = [Fluid('water', 1000.0, 1.0e-3, 1.0e25),...
          Fluid('DNAPL', 1630.0, 0.9e-3, 1.0e25)];

% Porous media properties
% -----------------------  |   K(m2)  | phi | biot|  Ks   | Slr |   Pb   | lambda |  relPerm   |  capPressure
sand1 = PorousMedia('sand1', 5.04e-10, 0.40,  1.0, 1.0e25, 0.078, 369.73,   3.86, 'BrooksCorey','BrooksCorey');
sand2 = PorousMedia('sand2', 2.05e-10, 0.39,  1.0, 1.0e25, 0.069, 434.45,   3.51, 'BrooksCorey','BrooksCorey');
sand3 = PorousMedia('sand3', 5.26e-11, 0.39,  1.0, 1.0e25, 0.098, 1323.95,  2.49, 'BrooksCorey','BrooksCorey');
sand4 = PorousMedia('sand4', 8.19e-12, 0.41,  1.0, 1.0e25, 0.189, 3246.15,  3.30, 'BrooksCorey','BrooksCorey');

% Activate gravity
sand1.gravityOn = true;
sand2.gravityOn = true;
sand3.gravityOn = true;
sand4.gravityOn = true;

rock = [sand1, sand2, sand3, sand4];

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluids',fluids);

% --- Boundary conditions -------------------------------------------------
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Liquid pressure boundary conditions
CoordSupp  = [1 0 Ly;1 Lx Ly];                
CoordLoad  = [];                     
CoordPresc = [];
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = zeros(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
CoordSupp  = [1 0 Ly;1 Lx Ly];                
CoordLoad  = [];                     
CoordPresc = [369.73 0 Ly;
              369.73 Lx Ly];
CoordInit  = [];    
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_pg = 369.73*ones(size(mdl.INITCOND_pg,1),1);

% Add prescribed gas pressure at the infiltration zone
tol = 1.0e-4;
reg = isInsideRectangle(mdl.NODE,[0.3-tol,0.5-tol],[0.4+tol,0.5+tol]);
mdl.SUPP_pg(reg == 1) = 1;
mdl.PRESCDISPL_pg(reg == 1) = 639.35;

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 3;

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% ========================= INITIALIZATION ===============================

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();

% Plot the mesh with the supports
mdl.plotMeshWithMatId();

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Transient analysis parameters
tinit = 1.0;   % Initial time
dt    = 1.0;   % Time step
tf    = 184;  % Final time

% Solve the problem
anl = Anl_TransientPicard(result);
anl.setUpTransientSolver(tinit,dt,tf,1.0,0.001,true);
anl.setPicardRelaxation();
anl.useRelativeError = true;
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');
