%% DESCRIPTION
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
[node, elem] = regularMesh(500, 270, 85, 45);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.mu = 1.0e-3;
water.K  = 2.15e9;
gas      = Fluid('gas');
gas.mu   = 0.45e-3;
gas.rho  = 660.0;
gas.K    = 1.5e9;

% Create porous media
rock_perm = PorousMedia('rock');
rock_perm.K                  = 9.86e-14;            % Intrinsic permeability (m2)
rock_perm.phi                = 0.2;                 % Porosity
rock_perm.Slr                = 0.01;                % Residual liquid saturation
rock_perm.Sgr                = 0.01;                % Residual gas saturation
rock_perm.Pb                 = 5.0e6;               % Gas-entry pressure
rock_perm.lambda                  = 2.0;                 % Curve-fitting parameter
rock_perm.m                  = 2.0;                 % Curve-fitting parameter
rock_perm.liqRelPermeability = 'PolynomialLiquid';  % Liquid relative permeability
rock_perm.gasRelPermeability = 'PolynomialGas';     % Gas relative permeability
rock.capillaryPressure  = 'BrooksCorey';  % Saturation degree function
% rock_perm.capillaryPressure  = 'Log';               % Saturation degree function

% Impermeable layer
rock_imp = copy(rock_perm);
rock_imp.K = rock_perm.K / 100;
rock_imp.Pb = 5.0e5;

% Compute element centroids to set material IDs
Xc = zeros(mdl.nelem, 2);
for i = 1:mdl.nelem
    xcoord = mdl.NODE(mdl.ELEM{i}, 1);
    ycoord = mdl.NODE(mdl.ELEM{i}, 2);
    xcentr = sum(xcoord) / 4;  % Considering linear quad elements
    ycentr = sum(ycoord) / 4;  % Considering linear quad elements
    Xc(i, :) = [xcentr, ycentr];
end

% Layer thickness
nlayers = 9;
layer_h = 270.0 / nlayers;
mdl.matID = ones(mdl.nelem,1);
for i = 1:nlayers
    if mod(i,2) == 0
        reg = isInsideRectangle(Xc, [0.0,layer_h * (i-1)], [500.0,layer_h * i]);
        mdl.matID(reg==1) = 2;
    end
end

% Set materials to model
mdl.setMaterial([rock_perm, rock_imp], water, gas);

%% BOUNDARY AND INITIAL CONDITIONS

% Capillary pressure function
% pc = @(Sl) -rock_perm.Pb*log(Sl);
pc = @(Sl) (rock_perm.Pb * ((Sl-rock_perm.Slr)/(1.0-rock_perm.Slr-rock_perm.Sgr))^(-1/rock_perm.lambda));

% Initial capillary pressure
pc_ini = pc(0.05);

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtBorder('right', 0.0);
mdl.setGasPressureDirichletBCAtBorder('right', pc_ini);

% Set Neumann boundary conditions
mdl.setPressureNeumannBCAtBorder('left', 9.42e-5); % m/s3

% Set initial conditions
mdl.setInitialGasPressureAtDomain(pc_ini);

%% PROCESS

day = 60*60*24;
year = 365*day;

% Analysis parameters
ti        = 0.0;        % Initial time
dt        = 0.001*day;      % Time step
tf        = 10*day;       % Final time
dtmax     = 5*day;        % Maximum time step
dtmin     = 0.01;  % Minimum time step
adaptStep = true;       % Adaptive step size

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.maxIter = 10;
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('LiquidSaturation');