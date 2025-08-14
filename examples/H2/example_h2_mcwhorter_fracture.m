%% DESCRIPTION
%
% McWhorter and Sunada problem using the Pl-Pg two-phase flow formulation
%
% References:
% * McWhorter and Sunada (1990). Exact integral solutions for two-phase flow. Water Resour Res, 26(3):399–413.
%
% Physics:
% * Two-phase hydraulic (H2)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_H2();
%% MESH

% Create mesh
Lx = 2.6;  % Horizontal dimension (m)
Ly = 2.0;  % Vertical dimension (m)
Nx = 30;  % Number of elements in the x-direction
Ny = 35;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
gas   = Fluid('gas');
gas.mu = 0.005;

% Create porous media
rock = PorousMedia('rock');
rock.K                  = 1.0e-10;        % Intrinsic permeability (m2)
rock.phi                = 0.15;           % Porosity
rock.Slr                = 0.02;           % Residual liquid saturation
rock.Sgr                = 0.001;          % Residual gas saturation
rock.Pb                 = 5.0e+3;         % Gas-entry pressure
rock.lambda             = 3.0;            % Curve-fitting parameter
rock.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
rock.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
rock.capillaryPressure  = 'BrooksCorey';  % Saturation degree function

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Brooks and Corey capillary pressure
pc = @(Sl) (rock.Pb * ((Sl-rock.Slr)/(1.0-rock.Slr-rock.Sgr))^(-1/rock.lambda));

% Initial capillary pressure
pc_ini = pc(0.05);

% Capillary pressure at the left border
pc_left = pc(0.8);

% Initial gas pressure
pg_ini = 100000.0;

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('left', pg_ini - pc_left);
mdl.setGasPressureDirichletBCAtBorder('left', pg_ini);

% Set initial conditions
mdl.setInitialPressureAtDomain(pg_ini - pc_ini);
mdl.setInitialGasPressureAtDomain(pg_ini);

%% DISCONTINUITIES

% Create discontinuities
Dx = [0.0; Lx];  % X-coordinates of polyline defining the fracture
Dy = [Ly/2; Ly/2];  % Y-coordinates of polyline defining the fracture
fracture = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fracture.porousMedia = rock;
fracture.liquidFluid = water;
fracture.gasFluid = gas;
fracture.initialAperture = 9.0e-3;

% Add fractures to model
mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

% Analysis parameters
ti        = 0.001;      % Initial time
dt        = 0.001;      % Time step
tf        = 1000;       % Final time
dtmax     = 5.0;        % Maximum time step
dtmin     = 0.0000001;  % Minimum time step
adaptStep = true;       % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 20;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('LiquidSaturation');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, 0.0];
mdl.plotFieldAlongSegment('LiquidSaturation', Xi, Xf, 500, 'x');
