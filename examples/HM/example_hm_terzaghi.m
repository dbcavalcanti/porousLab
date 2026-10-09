%% DESCRIPTION
%
% Terzaghi consolidation problem.
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
Lx = 0.1;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 5;    % Number of elements in the x-direction
Ny = 50;   % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.K = 2.2e9;  % Fluid bulk modulus (Pa)

% Create porous media
rock = PorousMedia('rock');
rock.K     = 1.15741e-16;   % Intrinsic permeability (m2)
rock.phi   = 0.3;           % Porosity
rock.Young = 8.64e9;        % Young modulus (Pa)
rock.nu    = 0.3;           % Poisson ratio
rock.Ks    = 36.0e9;        % Solid grain bulk modulus (Pa)

% Biot coefficient from the drained and grain bulk modulus
Kd = rock.Young/(3*(1-2*rock.nu));
rock.biot = 1-Kd/rock.Ks;  % 0.8

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);

% Loads
q = 1.0e4;  % Applied compressive traction (Pa)
mdl.addLoadAtBorder('top', 2, -q);

% Pressure
mdl.setPressureDirichletBCAtBorder('top', 0.0);

%% PROCESS

% Analysis parameters
ti = 0.0;    % Initial time
dt = 1.0;    % Time step
tf = 100.0;  % Final time

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(ti, dt, tf);
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Pressure');

% Final numerical profiles at the left-edge nodes, from bottom to top
comparisonNodes = mdl.getNodesAtBorder('left');
[y, order] = sort(mdl.NODE(comparisonNodes,2));
comparisonNodes = comparisonNodes(order);
pressureNum = mdl.U(mdl.ID(comparisonNodes,3));
displacementNum = mdl.U(mdl.ID(comparisonNodes,2));

% Analytical solution at the same time and nodes
comparisonTime = tf-ti;
alpha = rock.biot;
Sm = (alpha-rock.phi)/rock.Ks + rock.phi/water.K;
[pressureRef, displacementRef] = terzaghiReference(y, comparisonTime, ...
    Ly, q, rock.Young, rock.nu, rock.K, water.mu, alpha, Sm);

% RMSE in Pa for pressure and mm for vertical displacement
pressureRMSE = sqrt(mean((pressureNum-pressureRef).^2));
displacementRMSE = 1000*sqrt(mean((displacementNum-displacementRef).^2));
fprintf('\nFinal-time comparison at t = %g s:\n', comparisonTime);
fprintf('Pressure RMSE: %.6g Pa\n', pressureRMSE);
fprintf('Vertical displacement RMSE: %.6g mm\n', displacementRMSE);

% Plot pressure profiles
Xi = [0.0, 0.0]; Xf = [0.0, Ly];
mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'y');
plot(pressureRef, y, 'ok');
xlabel('Pore pressure (Pa)');
ylabel('Height, y (m)');
title(['t = ',num2str(comparisonTime), ' s']);
legend('PorousLab', 'Analytical', 'Location', 'best');

% Plot vertical displacement profiles
mdl.plotFieldAlongSegment('Uy', Xi, Xf, 500, 'y');
plot(displacementRef, y, 'ok');
xlabel('Vertical displacement, u_y (m)');
ylabel('Height, y (m)');
title(['t = ',num2str(comparisonTime), ' s']);
legend('PorousLab', 'Analytical', 'Location', 'best');
