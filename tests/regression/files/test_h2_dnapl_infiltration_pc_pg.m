%% DESCRIPTION
%
% Dense non-liquid phase infiltration problem using the Pc-Pg two-phase flow formulation.
%
% References:
% * Wang et al (2015). A parallel finite element method for two-phase flow processes in porous media: OpenGeoSys with PETSc. Environ Earth Sci, 73:2269–2285.
% * Bastian et al (2007). Numerical simulation and experimental studies of unsaturated water flow in heterogeneous systems. Reactive Flows, Diffusion and Transport. https://doi.org/10.1007/978-3-540-28396-6_22
%
% Physics:
% * Two-phase hydraulic (H2)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_H2_PcPg();

% Set model options
mdl.massLumping  = true;  % Diagonalize compressibility matrix (mass lumping)
mdl.lumpStrategy = 2;
mdl.gravityOn    = true;

%% MESH

% Create mesh
[node, elem] = regularMesh(0.9, 0.65, 60, 40);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water     = Fluid('water');
dnapl     = Fluid('dnapl');
dnapl.rho = 1.460e+3;  % Density (kg/m3)
dnapl.mu  = 0.900e-3;  % Viscosity (Pa*s)

% Create porous media
mat1 = PorousMedia('rock');
mat1.K                  = 6.64e-11;       % Intrinsic permeability (m2)
mat1.phi                = 0.4;            % Porosity
mat1.biot               = 1.0;            % Biot's coefficient
mat1.Ks                 = 1.0e+25;        % Solid bulk modulus (Pa)
mat1.Slr                = 0.1;            % Residual liquid saturation
mat1.Sgr                = 0.0;            % Residual gas saturation
mat1.Pb                 = 755.0;          % Gas-entry pressure
mat1.lambda             = 2.7;            % Curve-fitting parameter
mat1.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
mat1.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
mat1.capillaryPressure  = 'BrooksCorey';  % Saturation degree function

mat2 = PorousMedia('rock');
mat2.K                  = 6.64e-13;       % Intrinsic permeability (m2)
mat2.phi                = 0.4;            % Porosity
mat2.biot               = 1.0;            % Biot's coefficient
mat2.Ks                 = 1.0e+25;        % Solid bulk modulus (Pa)
mat2.Slr                = 0.1;            % Residual liquid saturation
mat2.Sgr                = 0.0;            % Residual gas saturation
mat2.Pb                 = 755.0;          % Gas-entry pressure
mat2.lambda             = 2.7;            % Curve-fitting parameter
mat2.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
mat2.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
mat2.capillaryPressure  = 'BrooksCorey';  % Saturation degree function

% Compute element centroids to set material IDs
nelem = 60 * 40;
Xc = zeros(nelem, 2);
for i = 1:nelem
    xcoord = mdl.NODE(mdl.ELEM{i}, 1);
    ycoord = mdl.NODE(mdl.ELEM{i}, 2);
    xcentr = sum(xcoord) / 4;  % Considering linear quad elements
    ycentr = sum(ycoord) / 4;  % Considering linear quad elements
    Xc(i, :) = [xcentr, ycentr];
end

% Porous media material 1
mdl.matID = ones(nelem, 1);

% Porous media material 2
reg = isInsideRectangle(Xc, [0.3,0.325], [0.45,0.4875]);
mdl.matID(reg==1) = 2;

% Set materials to model
mdl.setMaterial([mat1, mat2], water, dnapl);

%% BOUNDARY AND INITIAL CONDITIONS

% Capillary pressure
mdl.setCapillaryPressureDirichletBCAtBorder('left'  , 755.0);
mdl.setCapillaryPressureDirichletBCAtBorder('right' , 755.0);
mdl.setCapillaryPressureDirichletBCAtBorder('bottom', 755.0);
mdl.setInitialCapillaryPressureAtDomain(755.0);

% Gas pressure (hydrostatic profile)
for i = 1:mdl.nnodes
    pg = 7886.5 - 9810.0 * mdl.NODE(i,2);
    mdl.setInitialGasPressureAtNode(i, pg);

    % Fix gas pressure at lateral borders
    if ((abs(mdl.NODE(i,1)) < 1.0e-12) || ((abs(mdl.NODE(i,1) - 0.9)) < 1.0e-12))
        mdl.setGasPressureDirichletBCAtNode(i, pg);
    end
end

% Set prescribed gas pressure to infiltration zone
qginj = 0.3 * 0.075 / dnapl.rho;
tol = 1.0e-4;
reg = find(isInsideRectangle(mdl.NODE, [0.3-tol,0.65-tol], [0.6+tol,0.65+tol]));
for i = 1:length(reg)
    mdl.setGasPressureNeumannBCAtNode(reg(i), qginj/length(reg));
end

%% PROCESS

% Analysis parameters
ti        = 1.0;    % Initial time
dt        = 1.0;    % Time step
tf        = 5.0;  % Final time
dtmax     = 20.0;   % Maximum time step
dtmin     = 0.01;   % Minimum time step
adaptStep = true;   % Adaptive step size

% Run analysis
anl = Anl_Transient("Picard");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.setRelativeConvergenceCriteria(true);
anl.echo = false;
anl.run(mdl);