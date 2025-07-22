%% DESCRIPTION
%
% Buckley-Leverett problem using the Pl-Pg two-phase flow formulation.
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
mdl.massLumping  = true;  % Diagonalize compressibility matrix (mass lumping)
mdl.lumpStrategy = 2;

%% MESH

% Create mesh
Lx = 301.95;  % Horizontal dimension (m)
Ly = 100.00;  % Vertical dimension (m)
Nx = 99;      % Number of elements in the x-direction
Ny = 1;       % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
gas   = Fluid('gas');

% Create porous media
rock = PorousMedia('rock');
rock.K                  = 1.0e-7;         % Intrinsic permeability (m2)
rock.phi                = 0.2;            % Porosity
rock.Slr                = 0.2;            % Residual liquid saturation
rock.Sgr                = 0.2;            % Residual gas saturation
rock.Pb                 = 0.0;            % Gas-entry pressure
rock.lambda             = 2.0;            % Curve-fitting parameter
rock.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
rock.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
rock.capillaryPressure  = 'UMAT';         % Saturation degree function

% Set user material values for capillary pressure vs. saturation law
%           Pc | Sl
SlPcUMAT = [3.0, 0.2;
            0.0, 0.8];
rock.setUMATCapillaryPressureCurve(SlPcUMAT);

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('left', 200000.0);
mdl.setGasPressureDirichletBCAtBorder('left', 200000.230769231);

% Set initial conditions
mdl.setInitialPressureAtDomain(200000.0);
mdl.setInitialGasPressureAtDomain(200003.0);

%% PROCESS

% Analysis parameters
day       = 60*60*24;  % Conversion from days to seconds
ti        = 0.1*day;   % Initial time
dt        = 0.1*day;   % Time step
tf        = 50.0*day;  % Final time
dtmax     = 0.1*day;   % Maximum time step
dtmin     = 0.1*day;   % Minimum time step
adaptStep = true;      % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [Lx, 0.0];
mdl.plotFieldAlongSegment('LiquidPressure', Xi, Xf, 500, 'x');
mdl.plotFieldAlongSegment('CapillaryPressure', Xi, Xf, 500, 'x');
mdl.plotFieldAlongSegment('GasPressure', Xi, Xf, 500, 'x');
