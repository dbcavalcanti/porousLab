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
mdl.physics = 'hydraulicTwoPhasePcPg';

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

% Gas compressibility
R = 8.3144598;      % Universal gas constant (J/(mol*K)
T = 293.15;         % Temperature (K)
Mg = 0.0289;        % Molar mass (kg/mol)
rhog = 1630.0;      % Density (kg/m3)
Kg = (rhog * T * R) /  Mg;

fluids = [Fluid('water', 1000.0, 1.0e-3, 1.0e25),...
          Fluid('DNAPL', 1630.0, 0.9e-3, Kg)];

% Porous media properties
%                          |   K(m2)   | phi | biot|  Ks   | Slr |   Pb   | lambda |  relPerm   |  capPressure
sand1 = PorousMedia('sand1', 5.041e-10, 0.40,  1.0, 1.0e25, 0.078, 369.73,   3.86, 'BrooksCorey','BrooksCorey');
sand2 = PorousMedia('sand2', 2.051e-10, 0.39,  1.0, 1.0e25, 0.069, 434.45,   3.51, 'BrooksCorey','BrooksCorey');
sand3 = PorousMedia('sand3', 5.261e-11, 0.39,  1.0, 1.0e25, 0.098, 1323.95,  2.49, 'BrooksCorey','BrooksCorey');
sand4 = PorousMedia('sand4', 8.191e-12, 0.41,  1.0, 1.0e25, 0.189, 3246.15,  3.30, 'BrooksCorey','BrooksCorey');

rock = [sand1, sand2, sand3, sand4];

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluids',fluids);

% --- Boundary conditions -------------------------------------------------
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Capillary pressure boundary conditions
CoordSupp  = [1 0 Ly;1 Lx Ly];                
CoordLoad  = [];                     
CoordPresc = [369.73 0 Ly;
              369.73 Lx Ly];
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = 369.73*ones(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
CoordSupp  = [1 0 Ly;1 Lx Ly];                
CoordLoad  = [];                     
CoordPresc = [369.73 0 Ly;
              369.73 Lx Ly];
CoordInit  = [];    
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = 369.73*ones(size(mdl.INITCOND_p,1),1);

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 2;

%% ========================= INITIALIZATION ===============================

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();

% Plot the mesh with the supports
mdl.plotMeshWithBC();
return

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Transient analysis parameters
tinit = 0.1;   % Initial time
dt    = 0.1;   % Time step
tf    = 1000;  % Final time

% Solve the problem
anl = Anl_TransientPicard(result);
anl.setUpTransientSolver(tinit,dt,tf,5.0,0.0000001,true);
anl.setPicardRelaxation();
anl.useRelativeError = false;
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
