%% DESCRIPTION
%
% Wijesinghe problem.
%
% Physics:
% * Single-phase flow hydro-mechanical (HM)
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% MODEL

% Create model
mdl = Model_HM();

%% MESH

% Create mesh
[node, elem] = regularMesh(25.0, 1.0, 50, 21, [], [], 'ISOQ4', 0.0, 0.8, 0.5, 0.95);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');

% Create porous media
rock = PorousMedia('rock');
rock.K     = 1.0e-21;        % Intrinsic permeability (m2)
rock.phi   = 0.001;         % Porosity
rock.Young = 60.0e+9;       % Young modulus (Pa)
rock.nu    = 0.0;           % Poisson ratio

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [0.0, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

% Loads
mdl.addLoadAtBorder('top', 2, -50.0e6);

% Pressure
initialPressure = 11.0e6;  % Pa
inletPressure = 11.9e6;    % Pa, applied after equilibration
mdl.setPressureDirichletBCAtDomain(initialPressure);

%% DISCONTINUITY

% Polyline that defines the discontinuity
Xd = [ 0.0  , 0.5;
       25.0 , 0.5];
fracture = Discontinuity(Xd, true);

% Set fracture material properties
fracture.cohesiveLaw     = 'elastic';
fracture.shearStiffness  = 100.0e9;
fracture.normalStiffness = 100.0e9;
fracture.initialAperture = 1.0e-5;
fracture.liquidFluid     = water;

% Add fractures to model
discontinuityData = struct('addTangentialStretchingMode', false, ...
                           'addRelRotationMode', true);
mdl.addPreExistingDiscontinuities(fracture,discontinuityData);

%% PROCESS

% Analysis parameters
ti = 0.0;    % Initial time
dt = 1.0;    % Time step
tf = 1.0;   % Final time

% Set transient analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf);

% Begin time measurement
caseTimer = tic;
anl.run(mdl);
initialIterations = anl.totalIterations;

% mdl.plotFieldAlongDiscontinuiy('Aperture',1,'x');
% mdl.plotFieldAlongDiscontinuiy('Pressure',1,'x');

% mdl.resetDisplacements();
mdl.resetPressureDirichletBC();
mdl.setUpdateAperture(true);

for i = 1:mdl.nnodes
    if (mdl.NODE(i,1)<1.0e-12 && (abs(mdl.NODE(i,2)-0.5))<0.03)
        mdl.setPressureDirichletBCAtNode(i,inletPressure);
    end
end
mdl.setPressureDirichletBCAtBorder('right', initialPressure);
mdl.updateDirichletBC();

% clear anl
% anl = Anl_Transient("Newton");
anl.setUpTransientSolver(0.0, 0.01, 2000.0, 1, 0.0001, true);

anl.run(mdl);

% Time measurement
elapsedTime = toc(caseTimer);
fprintf('Elapsed time: %.3f seconds\n', elapsedTime);
fprintf('Total iterations: %d\n', initialIterations + anl.totalIterations);

%% POST-PROCESS

% Plot mesh and fracture
mdl.plotField('Model'); hold on;
fracture.plotIntersectedGeometry();

% Plot graphs
Xi = [0.0, 0.0]; Xf = [0.0, 1.0];
% mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'x');
mdl.plotFieldAlongDiscontinuiy('Aperture',1,'x');
mdl.plotFieldAlongDiscontinuiy('Pressure',1,'x');

% RMSE at the final time, measured from the inlet pressure step.
% Run only after the transient analysis has converged to anl.tf.
comparisonTime = anl.tf - anl.ti;
disc = mdl.discontinuitySet(1);
nSegments = disc.getNumberOfDiscontinuitySegments();
xError = []; pressureNum = []; apertureNum_mm = [];
for j = 1:nSegments
    seg = disc.segment(j);
    host = mdl.element(disc.elemID(j)).type;
    % Sample both fields at the same fracture integration points.
    [Xip, bj_mm] = seg.getAperture();
    pj = seg.getPorePressure(Xip, host);
    xError = [xError; Xip(:,1)-disc.X(1,1)]; %#ok<AGROW>
    pressureNum = [pressureNum; pj]; %#ok<AGROW>
    apertureNum_mm = [apertureNum_mm; bj_mm]; %#ok<AGROW>
end
% Each segment integration point is one sample, including shared endpoints.
% Aperture updates were disabled during equilibration, so bi is unchanged.
[apertureRef_mm, pressureRef] = wijesingheReference(xError, comparisonTime, ...
    water.mu, fracture.normalStiffness, fracture.initialAperture, ...
    initialPressure, inletPressure);

% RMSE in MPa for pressure and mm for aperture, as in the manuscript.
pressureRMSE = 1.0e-6*sqrt(mean((pressureNum-pressureRef).^2));
apertureRMSE = sqrt(mean((apertureNum_mm-apertureRef_mm).^2));
finalTimeErrors = table(comparisonTime, pressureRMSE, apertureRMSE, ...
    'VariableNames', {'Time_s','PressureRMSE_MPa','ApertureRMSE_mm'});
fprintf('\nRMSE against the semi-analytical solution:\n');
disp(finalTimeErrors);
