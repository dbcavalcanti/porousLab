%% DESCRIPTION
%
% OGS-5 Liakopoulos problem using the Pc-Pg two-phase flow formulation.
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
Lx = 0.1;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 3;    % Number of elements in the x-direction
Ny = 24;   % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water  = Fluid('water');
gas    = IdealGas('gas');
gas.mu = 1.8e-5;  % Viscosity (Pa*s)
gas.K  = 1.0e25;  % Compressibility/Bulk modulus (1/Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K                  = 4.5e-13;        % Intrinsic permeability (m2)
rock.phi                = 0.2975;         % Porosity
rock.biot               = 1.0;            % Biot's coefficient
rock.Ks                 = 1.0e25;         % Solid bulk modulus (Pa)
rock.Slr                = 0.0;            % Residual liquid saturation
rock.Sgr                = 0.2;            % Residual gas saturation
rock.Pb                 = 0.0;            % Gas-entry pressure
rock.lambda             = 3.0;            % Curve-fitting parameter
rock.liqRelPermeability = 'Liakopoulos';  % Liquid relative permeability
rock.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
rock.capillaryPressure  = 'Liakopoulos';  % Saturation degree function
rock.setMinLiquidRelPermeability(0.0001);
rock.setMinGasRelPermeability(0.0001);

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Set Dirichlet boundary conditions
mdl.setCapillaryPressureDirichletBCAtBorder('bottom', 0.0);
mdl.setGasPressureDirichletBCAtBorder('bottom', 101325.0);
mdl.setGasPressureDirichletBCAtBorder('top', 101325.0);

% Set initial conditions
mdl.setInitialCapillaryPressureAtDomain(0.0);
mdl.setInitialGasPressureAtDomain(101325.0);

%% PROCESS

% Analysis parameters
min       = 60;         % Conversion from days to seconds
ti        = 0.001*min;  % Initial time
dt        = 0.001*min;  % Time step
tf        = 5.0*min;    % Final time
dtmax     = 1.0*min;    % Maximum time step
dtmin     = 0.001*min;  % Minimum time step
adaptStep = true;       % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [0.0, Ly];
mdl.plotFieldAlongSegment('LiquidPressure', Xi, Xf, 500, 'y');
mdl.plotFieldAlongSegment('CapillaryPressure', Xi, Xf, 500, 'y');
mdl.plotFieldAlongSegment('GasPressure', Xi, Xf, 500, 'y');
