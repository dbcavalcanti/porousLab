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
Nx = 51;  % Number of elements in the x-direction
Ny = 51;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
gas   = Fluid('gas');
gas.mu = 5.7e-7;
gas.rho = 1460.0;

% Create porous media
rock = PorousMedia('rock');
rock.K                  = 1.0e-12;        % Intrinsic permeability (m2)
rock.phi                = 0.2;            % Porosity
rock.Slr                = 0.0;            % Residual liquid saturation
rock.Sgr                = 0.0;            % Residual gas saturation
rock.Pb                 = 1.0e+4;         % Gas-entry pressure
rock.lambda             = 2.0;            % Curve-fitting parameter
rock.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
rock.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
rock.capillaryPressure  = 'BrooksCorey';  % Saturation degree function
rock.setMinGasRelPermeability(1.0e-4);
rock.setMinLiquidRelPermeability(1.0e-4);

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Brooks and Corey capillary pressure
pc = @(Sl) (rock.Pb * ((Sl-rock.Slr)/(1.0-rock.Slr-rock.Sgr))^(-1/rock.lambda));

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('left',  2.0e5, [0.0, 0.1]);
mdl.setPressureDirichletBCAtBorder('right', 1.0e5, [0.9, 1.0]);
mdl.setGasPressureDirichletBCAtBorder('left',  2.0e5+pc(0.001), [0.0, 0.1]);
mdl.setGasPressureDirichletBCAtBorder('right',  1.0e5+pc(0.99), [0.9, 1.0]);

% Initial conditions
mdl.setInitialPressureAtDomain(1.0e5);
mdl.setInitialGasPressureAtDomain(1.0e5+pc(0.99));

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
% mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

% Analysis parameters
ti        = 0.0;      % Initial time
dt        = 0.01;      % Time step
tf        = 1000;      % Final time
dtmax     = 10;        % Maximum time step
dtmin     = 0.0000001;  % Minimum time step
adaptStep = true;       % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 5;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('GasSaturation');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, Ly];
mdl.plotFieldAlongSegment('GasSaturation', Xi, Xf, 500, 'x');
