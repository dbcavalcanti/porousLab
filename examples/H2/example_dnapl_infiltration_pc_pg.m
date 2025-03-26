%% ================ Two-Phase flow in porous media ========================
%
% References:
%
% Wang, W., Fischer, T., Zehner, B. et al.
% A parallel finite element method for two-phase flow processes in
% porous media: OpenGeoSys with PETSc. 
% Environ Earth Sci 73, 2269–2285 (2015).
% https://doi.org/10.1007/s12665-014-3576-z
%
% Bastian, P., Ippisch, O., Rezanezhad, F., Vogel, H.J., Roth, K. (2007).
% Numerical Simulation and Experimental Studies of Unsaturated Water Flow
% in Heterogeneous Systems. 
% In: Jäger, W., Rannacher, R., Warnatz, J. (eds) Reactive Flows,
% Diffusion and Transport. Springer, Berlin, Heidelberg.
% https://doi.org/10.1007/978-3-540-28396-6_22
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

mdl = Model_H2_PcPg();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 0.9;       % Horizontal dimension (m)
Ly = 0.65;      % Vertical dimension (m)
Nx = 60;        % Number of elements in the x-direction
Ny = 40;        % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Type of elements
mdl.type = 'ISOQ4';

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
reg = isInsideRectangle(Xc,[0.3,0.325],[0.45,0.4875]); mdl.matID(reg==1) = 2;

% --- Material properties of the domain -----------------------------------

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,1.0e25);
dnapl = Fluid('dnapl',1460.0,0.9e-3,1.0e25);

mat1 = PorousMedia('rock',6.64e-11,0.4,1.0,1.0e25,0.1,0.0,755.0,2.7,'BrooksCorey','BrooksCorey','BrooksCorey');
mat2 = PorousMedia('rock',6.64e-13,0.4,1.0,1.0e25,0.1,0.0,755.0,2.7,'BrooksCorey','BrooksCorey','BrooksCorey');

% Activate gravity
mat1.gravityOn = true;
mat2.gravityOn = true;

rock = [mat1 , mat2];

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'liquidFluid',water,...
    'gasFluid',dnapl);

% --- Boundary conditions -------------------------------------------------
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Capillary pressure boundary conditions
CoordSupp  = [1 0 -1;       % Left border
              1 Lx -1;      % Right border
              1 -1 0.0];    % Bottom border                        
CoordLoad  = [];                      
CoordPresc = [755.0 0 -1;
              755.0 Lx -1;
              755.0 -1 0.0];            
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);
mdl.INITCOND_p = 755.0*ones(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
CoordSupp  = [1 0 -1;       % Left border
              1 Lx -1];     % Right border                          
CoordLoad  = [];                      
CoordPresc = [];             
CoordInit  = [];                      
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

% Set hydrostatic pressure 
for i = 1:size(mdl.NODE,1)
    pg = 7886.5 - 9810.0*mdl.NODE(i,2);
    mdl.INITCOND_pg(i) = pg;
    if (mdl.SUPP_pg(i,1) == 1)
        mdl.PRESCDISPL_pg(i) = pg;
    end
end

% Add prescribed gas pressure at the infiltration zone
qginj = 0.3 * 0.075 / dnapl.rho;
tol = 1.0e-4;
reg = isInsideRectangle(mdl.NODE,[0.3-tol,Ly-tol],[0.6+tol,Ly+tol]);
mdl.LOAD_pg(reg == 1) = qginj/sum(reg);

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 2;

% Diagonalize compressibility matrix
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% ========================= INITIALIZATION ===============================

% Plot the mesh with the supports
mdl.plotMeshWithMatId();

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Transient analysis parameters
tinit = 1.0;       % Initial time
dt    = 1.0;       % Time step
tf    = 800.0;      % Final time
dtmax = 20.0;
dtmin = 0.01;

% Solve the problem
anl = Anl_TransientPicard(result);
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.setPicardRelaxation();
anl.useRelativeError = true;
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , Ly];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
% mdl.plotField('CapillaryPressure');
% mdl.plotField('GasPressure');
mdl.plotField('LiquidSaturation');
mdl.plotField('GasSaturation');