%% DESCRIPTION
%
% Liakopoulos problem using the Pc-Pg two-phase flow formulation.
%
% References:
% * Schrefler and Xiaoyong (1993). A fully coupled model for water flow and airflow in deformable porous media. Water Resour Res, 29(1):155–167
% * Schrefler and Scotta (2001). A fully coupled dynamic model for two-phase fluid flow in deformable porous media. CMAME, 190(24-25):3223-3246.
%
% Physics:
% * Two-phase flow hydro-mechanical (H2M)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_H2M();

mdl.gravityOn = true;
mdl.intOrder = 3;

%% MESH

% Create mesh
Lx = 0.1;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 10;    % Number of elements in the x-direction
Ny = 60;   % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);
[node, elem] = convertToQuadraticMesh(node, elem);

% Set mesh to model
mdl.setMesh(node, elem);
mdl.resequenceNodes();

%% MATERIALS

% Create fluids
water = Fluid('water');
water.K = 0.43e+13;  % Compressibility/Bulk modulus (1/Pa)

gas = Fluid('gas');
gas.mu  = 1.0e-3;  % Viscosity (Pa*s)
gas.rho = 1.22;
gas.K = 0.1e6;

% Create the porous media
rock = PorousMedia('rock');
rock.K                  = 0.46e-11;       % Intrinsic permeability (m2)
rock.phi                = 0.3;            % Porosity
rock.Ks                 = 0.14e+10;       % Solid bulk modulus (Pa)
rock.Slr                = 0.3966;         % Residual liquid saturation 
rock.lambda             = 3.0;            % Curve-fitting parameter
rock.Pb                 = 2.25e+5;        % Gas-entry pressure (Pa)
rock.Young              = 6.0e+6;         % Young modulus (Pa)
rock.nu                 = 0.4;            % Poisson ratio
rock.rho                = 2000.0;         % Density (kg/m3)
rock.liqRelPermeability = 'BrooksCorey';  % Liquid relative permeability
rock.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
rock.capillaryPressure  = 'BrooksCorey';  % Saturation degree function
rock.setMinLiquidRelPermeability(1.0e-4);
rock.setMinGasRelPermeability(1.0e-4);

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

% Loads
mdl.addLoadAtBorder('top', 2, -1.0e3);

% Initial capillary pressure
pci = 3.819e+5;

% Initial gas pressur
pgi = pci - 280.0e3;

% Liquid pressure
mdl.setPressureDirichletBCAtBorder('top', -420.0e3);
mdl.setInitialPressureAtDomain(-280.0e3);

% Gas pressure
mdl.setGasPressureDirichletBCAtBorder('top', 102.0e3);
% mdl.setGasPressureDirichletBCAtBorder('bottom', 0.0);
mdl.setInitialGasPressureAtDomain(pgi);

%% PROCESS

days = 60*60*24;

% Analysis parameters
ti        = 0.0;     % Initial time
dt        = 0.01*days;     % Time step
tf        = 0.01*days;     % Final time
dtmax     = 0.01 * days;     % Maximum time step
dtmin     = 0.0001;  % Minimum time step
adaptStep = true;    % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.setRelativeConvergenceCriteria(true);
anl.maxIter = 15;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('GasSaturation');
mdl.plotField('GasPressure');

% Plot graphs
Xi = [Lx/2.0, 0.0]; Xf = [Lx/2.0, Ly];
mdl.plotFieldAlongSegment('LiquidPressure', Xi, Xf, 500, 'y');
% mdl.plotFieldAlongSegment('CapillaryPressure', Xi, Xf, 500, 'y');
% mdl.plotFieldAlongSegment('GasPressure', Xi, Xf, 500, 'y');
