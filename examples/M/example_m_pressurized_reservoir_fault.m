%% DESCRIPTION
%
% Uniform traction on a plate with isotropic damage model.
%
% Physics:
% * Mechanical (M)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_M();

% Set model options
mdl.condenseEnrDofs   = false;
mdl.addPorePressure   = true;
mdl.subDivIntegration = true;
mdl.symmetricSDAEFEM  = false;

%% MESH

% Create mesh
Lx = 4000.0;  % Horizontal dimension (m)
Ly = 1000.0;  % Vertical dimension (m)
Nx = 51;      % Number of elements in the x-direction
Ny = 25;      % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny, [], [], 'ISOQ4', 0.5, 0.7, 0.5, 0.7);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Earth pressure coefficient
K0 = 0.65;

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 30.0e+9;        % Young modulus (Pa)
rock.nu    = K0/(1+K0);      % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% DISCONTINUITIES

% Create discontinuities 
Xd = [ 0.5*Lx , 0.0;
       0.5*Lx , Ly ];
fault = Discontinuity(Xd, true);

% Set fracture material properties
fault.cohesiveLaw     = 'elastic';
fault.shearStiffness  = 1.0e15;       % Pa/m
fault.normalStiffness = 1.0e15;       % Pa/m

% Add fractures to model
discontinuityData = struct('addStretchingMode', false, 'addRelRotationMode', false);
mdl.addPreExistingDiscontinuities(fault, discontinuityData);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

% Pressure load
mdl.addLoadAtBorder('top', 2, -70.0e6);

% Set external pore-pressure
P = 35.0e6 * ones(mdl.nnodes,1);
mdl.setPorePressureField(P);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

% Update the pressure at the left side
tol = 1.0e-5;
reservoir = isInsideRectangle(mdl.NODE, [0.0-tol,350.0-tol], [0.5*Lx+tol,650.0+tol]);
P(reservoir == 1) = P(reservoir == 1) + 20.0e6;

% Update pressure at the right side
offset = 0.0;
reservoir = isInsideRectangle(mdl.NODE, [0.5*Lx,350.0+offset-tol], [Lx+tol,650.0+offset+tol]);
P(reservoir == 1) = P(reservoir == 1) + 20.0e6;
mdl.setPorePressureField(P);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

% Plot stresses
mdl.plotField('Sx');
mdl.plotField('Sy');
mdl.plotField('Sxy');
hold on;
fault.plotIntersectedGeometry();

mdl.plotFieldAlongDiscontinuiy('Sn',1,'y');
mdl.plotFieldAlongSegment('Sx',[0.5*Lx,0.0],[0.5*Lx,Ly],100,'y')
