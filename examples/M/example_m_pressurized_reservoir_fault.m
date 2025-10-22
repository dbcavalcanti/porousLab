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
mdl.symmetricSDAEFEM  = false;

%% MESH

% Create mesh
Lx = 4000.0;  % Horizontal dimension (m)
Ly = 1000.0;  % Vertical dimension (m)
Nx = 60;      % Number of elements in the x-direction
Ny = 60;      % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny, [], [], 'ISOQ4', 0.5, 0.7, 0.5, 0.8);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Earth pressure coefficient
% K0 = 0.60;

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 11868.0e+6;        % Young modulus (Pa)
% rock.nu    = K0/(1+K0);         % Poisson ratio
rock.nu    = 0.29;         % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% DISCONTINUITIES

% Fault center
Xdc = 0.5 * [Lx, Ly];

% Fault dip
dip = 60.0;
dxdy = 0.0;
if (dip-90)<1.0e-8
    dxdy = 1.0/tand(dip);
end

% Create discontinuities 
Xd = [ Xdc(1)+(0.0 - Xdc(2))*dxdy , 0.0;
       Xdc(1)+(Ly - Xdc(2))*dxdy , Ly ];
fault = Discontinuity(Xd, true);

% Set fracture material properties
fault.cohesiveLaw     = 'elastic';
fault.shearStiffness  = 1.0e15;       % Pa/m
fault.normalStiffness = 1.0e15;       % Pa/m

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', false, 'addNormalStretchingMode', false, 'addRelRotationMode', false);
mdl.addPreExistingDiscontinuities(fault);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

% Pressure load
mdl.addLoadAtBorder('top', 2, -70.0e6);

% Set external pore-pressure
P0 = 35.0e6;
P = P0 * ones(mdl.nnodes,1);
mdl.setPorePressureField(P);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

mdl.resetDisplacements();

% Parameters
tol    = 1.0e-5;

% Update pressure at the reservoir
reservoir = isInsideRectangle(mdl.NODE, [0.0-tol,350.0-tol], [Lx+tol,650.0+tol]);
P(reservoir == 1) = P0 + 20.0e6;

% Update values
mdl.setPorePressureField(P);

% Analysis parameters
clear anl

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(0.0, 1.0, 1.0);
anl.run(mdl);

%% POST-PROCESS

mdl.plotField('PressureExt');
hold on;
fault.plotIntersectedGeometry();

% % Plot stresses
% mdl.plotField('Sx');
% hold on;
% fault.plotIntersectedGeometry();
% mdl.plotField('Sy');
% hold on;
% fault.plotIntersectedGeometry();
% mdl.plotField('Sxy');
% hold on;
% fault.plotIntersectedGeometry();
mdl.plotFieldAlongSegment('Sy',[0.5*Lx,0.0],[0.5*Lx,Ly],100,'y');
mdl.plotFieldAlongSegment('Sx',[0.5*Lx,0.0],[0.5*Lx,Ly],100,'y');
mdl.plotFieldAlongDiscontinuiy('St',1,'y');
mdl.plotFieldAlongDiscontinuiy('Sn',1,'y');
mdl.plotField('Uy');
mdl.plotField('Ux');

