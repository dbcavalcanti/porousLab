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
Lx = 1.0;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 40;  % Number of elements in the x-direction
Ny = 40;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.mu = 0.8e-3;
gas   = Fluid('gas');
gas.mu = 0.45e-3;
gas.rho = 660.0;

% Saturation law
pc_curve = 101325 * [1.0, 0.6, 0.4 ,0.24, 0.18, 0.14,0.1,0.07,0.04,0.02,0.0];
sl_curve = linspace(0.0, 1.0, 11);
PcSl = [pc_curve', sl_curve'];

% Create porous media
rock = PorousMedia('rock');
rock.K                  = 9.87e-16;       % Intrinsic permeability (m2)
rock.phi                = 0.2;            % Porosity
rock.liqRelPermeability = 'PolynomialLiquid';  % Liquid relative permeability
rock.gasRelPermeability = 'PolynomialGas';  % Gas relative permeability
rock.capillaryPressure  = 'UMAT';           % Saturation degree function
rock.SlPc_umat = PcSl;
rock.setMinGasRelPermeability(1.0e-5);
rock.setMinLiquidRelPermeability(1.0e-5);

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Set boundary conditions
mdl.setPressureNeumannBCAtPoint([0.0, 0.0], 2.3148e-8);  %m3/s
% mdl.setPressureNeumannBCAtPoint([1.0, 1.0], -2.3148e-8);
mdl.setPressureDirichletBCAtPoint([1.0, 1.0], 0.0);
mdl.setGasPressureDirichletBCAtPoint([1.0, 1.0], 101225.0);

% Initial conditions
mdl.setInitialGasPressureAtDomain(101225.0);

%% DISCONTINUITIES

% Create discontinuities
Xd = [0.2, 0.2;
      0.8, 0.8]; 
fracture = Discontinuity(Xd, true);

% Set fracture material properties
fracture.porousMedia = rock;
fracture.liquidFluid = water;
fracture.gasFluid = gas;
fracture.porosity = 1.0;
fracture.initialAperture = 1.0e-2;

% Add fractures to model
mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

days = 60*60*24;

% Analysis parameters
ti        = 0.0;      % Initial time
dt        = 0.001*days;      % Time step
tf        = 40*days;      % Final time
dtmax     = 1*days;        % Maximum time step
dtmin     = 0.0000001;  % Minimum time step
adaptStep = true;       % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('LiquidSaturation');
mdl.plotField('GasSaturation');
hold on;
fracture.plotIntersectedGeometry();

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, Ly];
mdl.plotFieldAlongSegment('LiquidSaturation', Xi, Xf, 500, 'x');
mdl.plotFieldAlongSegment('GasSaturation', Xi, Xf, 500, 'x');
