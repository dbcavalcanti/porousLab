%% DESCRIPTION
%
% Fluid-flow through the foundation of a gravity dam.
%
% References:
% * Segura and Carol (2004). On zero-thickness interface elements for diffusion problems. Int J Numer Anal Methods Geomech, 28(9):947-962.
%
% Physics:
% * Single-phase hydraulic (H)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_H();

mdl.intOrder = 3;

%% MESH

% Create mesh
Lx = 24.0;  % Horizontal dimension (m)
Ly = 6.0;   % Vertical dimension (m)
Nx = 90;  % Number of elements in the x-direction
Ny = 35;  % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.K = 2.2e+9;  % Compressibility/Bulk modulus (1/Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K   = 1.0194e-14;  % Intrinsic permeability (m2)
rock.phi = 0.3;         % Porosity

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY CONDITIONS

% Set Dirichlet boundary conditions
for i = 1:size(mdl.NODE,1)
    if ((mdl.NODE(i,1) < 8.0) && (abs(mdl.NODE(i,2)-Ly) < 1.0e-9))
        mdl.setPressureDirichletBCAtNode(i, 120.0e3);
    end
    if ((mdl.NODE(i,1) > 12.0) && (abs(mdl.NODE(i,2)-Ly) < 1.0e-9))
        mdl.setPressureDirichletBCAtNode(i, 60.0e3);
    end
end

%% DISCONTINUITIES

% Generate fractures
% FractureDataDamFoundation = generateRandomFractures(Lx, Ly, 20, 2*Lx/Nx, 'Seed', 123);
FractureDataDamFoundation = generateParallelFractures(Lx, Ly, -pi/4.0, 1.0);
% FractureDataDamFoundation2 = generateParallelFractures(Lx, Ly,  pi/3.0, 3.0);
% FractureDataDamFoundation = [FractureDataDamFoundation1; FractureDataDamFoundation2];

% Number of discontinuities
nd = length(FractureDataDamFoundation);

% Create discontinuities
fractures(1,nd) = Discontinuity();
for i = 1:nd
    fractures(i) = Discontinuity(FractureDataDamFoundation{i}, false);
end

% Set fracture material properties
rng(123,'twister');
for i = 1:nd
    fractures(i).liquidFluid = water;
    % fractures(i).initialAperture = 1e-5 + (5e-4 - 1e-5) .* rand;
    fractures(i).initialAperture = 1e-4;
end

% Add fractures to model
mdl.addPreExistingDiscontinuities(fractures);

%% PROCESS

% Analysis parameters
ti        = 1.0;    % Initial time
dt        = 1.0;    % Time step
tf        = 500.0;  % Final time
dtmax     = 50.0;   % Maximum time step
dtmin     = 0.001;  % Minimum time step
adaptStep = true;   % Adaptive step size

% Run analysis
% anl = Anl_Linear();
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
mdl.printResults();

% Plot model
mdl.plotField('Model');
colorbar off; hold on;
for i = 1:length(mdl.discontinuitySet)
    mdl.discontinuitySet(i).plotIntersectedGeometry();
end

% Plot contours
mdl.plotField('Pressure',[50e3,120e3]);
hold on;
for i = 1:length(mdl.discontinuitySet)
    mdl.discontinuitySet(i).plotIntersectedGeometry();
end

% Plot graphs
Xi = [0.0, Ly/2.0]; Xf = [Lx, Ly/2.0];
mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'x');
