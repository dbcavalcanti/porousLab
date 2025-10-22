%% DESCRIPTION
%
% Strip footing problem.
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
mdl.isPlaneStress = true;
mdl.symmetricSDAEFEM  = false;

%% MESH

% Load mesh
[node, elem] = regularMesh(1.0, 1.0, 51, 51);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.Young         = 1.0e6;            % Young modulus (kPa)
rock.nu            = 0.3;              % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('top', [NaN, -0.1]);
mdl.setDisplacementDirichletBCAtPoint([1.0, 0.0], [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtPoint([1.0, 1.0], [0.0, NaN]);

% Pressure load
% mdl.addLoadAtBorder('top', 2, -1.0e5);
% mdl.addLoadAtBorder('left', 1, 0.5);

%% DISCONTINUITIES

% Create discontinuities 
Xd = [0.0 , 0.4586;
      1.0 , 0.6586];
fracture = Discontinuity(Xd, true);

% Set fracture material properties
fracture.cohesiveLaw     = 'mohrCoulomb';
fracture.initialAperture = 0.0;
fracture.shearStiffness  = 1.0e13;
fracture.normalStiffness = 1.0e13;
fracture.frictionAngle   = atan(0.19);
fracture.dilationAngle   = atan(0.19);
fracture.cohesion        = 10.0;
fracture.tensionCutOff   = 0.0;

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', false, 'addNormalStretchingMode', false, 'addRelRotationMode', false);
mdl.addPreExistingDiscontinuities(fracture, discontinuityData);

%% PROCESS

% % Setup analysis
% anl = Anl_NonlinearQuasiStatic('GeneralizedDisplacement');
% anl.adjustStep    = true;
% anl.increment     = 1.0;
% anl.max_increment = 10.0;
% anl.max_lratio    = 1.0;
% anl.max_step      = 100;
% anl.max_iter      = 100;
% anl.trg_iter      = 3;
% 
% % Node and DOF used to plot Load Factor vs Displacement
% ndId = mdl.closestNodeToPoint([0.0, 1.0]);
% anl.setPlotDof(ndId, 2);
% 
% % Run analysis
% anl.run(mdl);

% Analysis parameters
ti = 0.0;    % Initial time
dt = 1.0;    % Time step
tf = 1.0;  % Final time

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf);
anl.run(mdl);

%% POST-PROCESS

% Plot Load Factor vs Displacement
% anl.plotCurves();

% mdl.plotDeformedMesh(1.0);

% Plot contours
mdl.plotField('Sy');
hold on;
fracture.plotIntersectedGeometry();
mdl.plotField('Sx');
mdl.plotField('Sxy');
mdl.plotField('Ux');
mdl.plotFieldAlongSegment('Sy',[0.5,0.0],[0.5,1.0],100,'y');
mdl.plotFieldAlongSegment('Ux',[0.5,0.0],[0.5,1.0],100,'y');
mdl.plotFieldAlongSegment('Uy',[0.0,1.0],[1.0,1.0],100,'y');
mdl.plotFieldAlongDiscontinuiy('St',1,'x');
mdl.plotFieldAlongDiscontinuiy('Sn',1,'x');


% Expected stresses along the discontinuity
sxx = -0.0;
syy = -1.0;
stress = [sxx , 0.0; 0.0, syy];
theta = atan(0.2);
cs  = cos(theta);
sn  = sin(theta);
R = [cs , sn; -sn , cs];
tnm = R * stress * R'
