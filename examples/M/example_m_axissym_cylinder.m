%% DESCRIPTION
%
% Elastic pressurized cilynder using axisymmetric model.
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

mdl.isAxisSymmetric = true;

%% MESH

% Create mesh
[node, elem] = regularMesh(0.1, 0.5, 10, 50);
node(:,1) = node(:,1) + 0.1;

[node, elem] = convertToQuadraticMesh(node, elem);

% Set mesh to model
mdl.setMesh(node, elem);
mdl.resequenceNodes();

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical = 'elastic';  % Constitutive law
rock.Young      = 30.0e6;     % Young modulus (Pa)
rock.nu         = 0.3;        % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('top'  ,  [NaN, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);
% mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

% Loads
pint = 1.0;
mdl.addLoadAtBorder('left', 1, pint);

%% PROCESS

% Run analysis
anl = Anl_Linear();
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Sx');
mdl.plotFieldAlongSegment('Sx',[0.1,0.0],[0.2,0.0])
