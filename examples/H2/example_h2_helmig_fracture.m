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

% Set model options
mdl.gravityOn  = true;

%% MESH

% Create mesh
Lx = 1.0;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 31;  % Number of elements in the x-direction
Ny = 31;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
gas   = IdealGas('gas');
gas.mu = 1.65e-5;
gas.T  = 293.2016;

% Create porous media
rock = PorousMedia('rock');
rock.K                  = 1.0e-12;        % Intrinsic permeability (m2)
rock.phi                = 0.1;            % Porosity
rock.Slr                = 0.0;            % Residual liquid saturation
rock.Sgr                = 0.0;            % Residual gas saturation
rock.Pb                 = 2.0e+3;         % Gas-entry pressure
rock.lambda             = 2.0;            % Curve-fitting parameter
rock.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
rock.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
rock.capillaryPressure  = 'BrooksCorey';  % Saturation degree function
rock.setMinGasRelPermeability(1.0e-4);
rock.setMinLiquidRelPermeability(1.0e-4);

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('top', 0.0);
mdl.setGasPressureDirichletBCAtBorder('top', 1000.0);

% Set Neumann boundary conditions
mdl.setGasPressureNeumannBCAtBorder('bottom', 0.00001);

% Initial conditions
mdl.setInitialGasPressureAtDomain(1000.0);

%% DISCONTINUITIES

% Create discontinuities
Xd = [0.5, 0.0;
      0.5, 1.0]; 
fracture = Discontinuity(Xd, true);

% Set fracture material properties
fracture.porousMedia = rock;
fracture.liquidFluid = water;
fracture.gasFluid = gas;
fracture.porosity = 0.3;
fracture.initialAperture = 5.0e-3;

% Add fractures to model
mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

% Analysis parameters
ti        = 0.00001;      % Initial time
dt        = 0.00001;      % Time step
tf        = 100;       % Final time
dtmax     = 2.0;        % Maximum time step
dtmin     = 0.0000001;  % Minimum time step
adaptStep = true;       % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('GasSaturation');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, 0.0];
mdl.plotFieldAlongSegment('LiquidSaturation', Xi, Xf, 500, 'x');
