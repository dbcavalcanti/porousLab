%% DESCRIPTION
%
% Dense non-liquid phase infiltration problem using the Pc-Pg two-phase 
% flow formulation
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
% Physics:
% * Two-phase hydraulic (H2)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% INITIALIZATION
close all; clear; clc;

% Path to source directory
src_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
addpath(genpath(src_dir));
print_header;

% Create model
mdl = Model_H2_PcPg();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 0.9;       % Horizontal dimension (m)
Ly = 0.65;      % Vertical dimension (m)
Nx = 60;        % Number of elements in the x-direction
Ny = 40;        % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

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
reg = isInsideRectangle(Xc,[0.3,0.325],[0.45,0.4875]);
mdl.matID(reg==1) = 2;

% --- Material properties of the domain -----------------------------------

% Create the fluids
water     = Fluid('water');
dnapl     = Fluid('dnapl');
dnapl.rho = 1460.0;
dnapl.mu  = 0.9e-3;

% Create the porous media
mat1 = PorousMedia('rock',6.64e-11,0.4,1.0,1.0e25,0.1,0.0,755.0,2.7,'BrooksCorey','BrooksCorey','BrooksCorey');
mat2 = PorousMedia('rock',6.64e-13,0.4,1.0,1.0e25,0.1,0.0,755.0,2.7,'BrooksCorey','BrooksCorey','BrooksCorey');

% Activate gravity
mat1.gravityOn = true;
mat2.gravityOn = true;

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',[mat1 , mat2], ...
    'liquidFluid',water,...
    'gasFluid',dnapl);

% --- Boundary and initial conditions -------------------------------------

% Capillary pressure conditions
mdl.setCapillaryPressureDirichletBCAtBorder('left'  ,755.0);
mdl.setCapillaryPressureDirichletBCAtBorder('right' ,755.0);
mdl.setCapillaryPressureDirichletBCAtBorder('bottom',755.0);
mdl.setInitialCapillaryPressureAtDomain(755.0);

% Gas pressure conditions
% It follows an hydrostatic profile
for i = 1:mdl.nnodes
    pg = 7886.5 - 9810.0*mdl.NODE(i,2);
    mdl.setInitialGasPressureAtNode(i,pg);
    % Fix the gas pressure at the lateral borders
    if ((abs(mdl.NODE(i,1))<1.0e-12) || ((abs(mdl.NODE(i,1) - Lx))<1.0e-12))
        mdl.setGasPressureDirichletBCAtNode(i, pg);
    end
end

% Add prescribed gas pressure at the infiltration zone
qginj = 0.3 * 0.075 / dnapl.rho;
tol = 1.0e-4;
reg = find(isInsideRectangle(mdl.NODE,[0.3-tol,Ly-tol],[0.6+tol,Ly+tol]));
for i = 1:length(reg)
    mdl.setGasPressureNeumannBCAtNode(reg(i), qginj/length(reg));
end

% --- Numerical model configuration ---------------------------------------

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% PROCESS

% Transient analysis parameters
tinit = 1.0;       % Initial time
dt    = 1.0;       % Time step
tf    = 800.0;      % Final time
dtmax = 20.0;
dtmin = 0.01;

% Solve the problem
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.setRelativeConvergenceCriteria(true);
anl.run(mdl);

%% POST-PROCESS

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , Ly];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('LiquidSaturation');
mdl.plotField('GasSaturation');