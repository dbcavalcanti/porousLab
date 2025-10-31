%% DESCRIPTION
%
% Block crossed by strong discontinuity.
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

%% MESH

% Create mesh
Lx = 5.0;  % Horizontal dimension (m)
Ly = 3.0;  % Vertical dimension (m)
Nx = 15*3;    % Number of elements in the x-direction
Ny = 9*3;    % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');

% Create porous media
rock = PorousMedia('rock');
rock.K   = 9.8e-16;  % Intrinsic permeability (m2)
rock.phi = 0.25;     % Porosity

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY CONDITIONS

% Set Dirichlet boundary conditions
mdl.setPressureDirichletBCAtPoint([0.0,Ly], 10.0);
% mdl.setPressureDirichletBCAtBorder('left', 10.0);
mdl.setPressureDirichletBCAtBorder('right', 0.0);


%% DISCONTINUITIES

fractures(1,2) = Discontinuity();

% Create discontinuities
Dx = [2.5; 2.5];  % X-coordinates of polyline defining the fracture
Dy = [0.0; Ly];   % Y-coordinates of polyline defining the fracture
fractures(1) = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fractures(1).fluid = water;
fractures(1).initialAperture = 1.0e-3;
fractures(1).transversalPerm = 0.0;
% fracture.leakoff = 1.0e-15;

% Create discontinuities
Dx = [0.0; Lx];  % X-coordinates of polyline defining the fracture
Dy = [1.5; 1.5];   % Y-coordinates of polyline defining the fracture
fractures(2) = Discontinuity([Dx, Dy], true);

% Set fracture material properties
fractures(2).fluid = water;
fractures(2).initialAperture = 1.0e-3;
fractures(2).transversalPerm = 0.0;
% fracture.leakoff = 1.0e-15;

% Add fractures to model
mdl.addPreExistingDiscontinuities(fractures);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Print results to command window
% mdl.printResults();

% Plot contours
mdl.plotField('Pressure');
hold on;
fractures(1).plotIntersectedGeometry();
fractures(2).plotIntersectedGeometry();
