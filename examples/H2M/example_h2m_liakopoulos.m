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

% Set model options
mdl.gravityOn    = true;

%% MESH

% Create mesh
Lx = 0.1;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 4;    % Number of elements in the x-direction
Ny = 40;   % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.K = 2.0e+9;  % Compressibility/Bulk modulus (1/Pa)

gas = IdealGas('gas');
gas.mu  = 1.8e-5;  % Viscosity (Pa*s)

% Create the porous media
rock = PorousMedia('rock');
rock.K                  = 4.5e-13;        % Intrinsic permeability (m2)
rock.phi                = 0.2975;         % Porosity
rock.Ks                 = 1.0e+12;        % Solid bulk modulus (Pa)
rock.Slr                = 0.2;            % Residual liquid saturation 
rock.lambda             = 3.0;            % Curve-fitting parameter
rock.Young              = 1.3e+6;         % Young modulus (Pa)
rock.nu                 = 0.4;            % Poisson ratio
rock.rho                = 2000.0;         % Density (kg/m3)
rock.liqRelPermeability = 'Liakopoulos';  % Liquid relative permeability
rock.gasRelPermeability = 'BrooksCorey';  % Gas relative permeability
rock.capillaryPressure  = 'Liakopoulos';  % Saturation degree function
rock.setMinLiquidRelPermeability(1.0e-4);
rock.setMinGasRelPermeability(1.0e-4);

% Set materials to model
mdl.setMaterial(rock, water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

% Liquid pressure
mdl.setPressureDirichletBCAtBorder('bottom', 101025.0);
mdl.setInitialPressureAtDomain(101025.0);

% Gas pressure
mdl.setGasPressureDirichletBCAtBorder('top',    101325.0);
mdl.setGasPressureDirichletBCAtBorder('bottom', 101325.0);
mdl.setInitialGasPressureAtDomain(101325.0);

%% PRE-PROCESS

mdl.preComputations();

% Initialize stress tensor
for i = 1:mdl.nelem
    for j = 1:mdl.element(i).type.nIntPoints
        mdl.element(i).type.intPoint(j).stressOld = [101325.0; 101325.0; 101325.0; 0.0];
    end
end

%% PROCESS

% Analysis parameters
ti        = 0.1;     % Initial time
dt        = 0.1;     % Time step
tf        = 1.0;     % Final time
dtmax     = 1.0;     % Maximum time step
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
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');

% Plot graphs
Xi = [0.0, 0.0]; Xf = [0.0, Ly];
mdl.plotFieldAlongSegment('LiquidPressure', Xi, Xf, 500, 'y');
mdl.plotFieldAlongSegment('CapillaryPressure', Xi, Xf, 500, 'y');
mdl.plotFieldAlongSegment('GasPressure', Xi, Xf, 500, 'y');
