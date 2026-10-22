%% DESCRIPTION
%
% Benchmark 1 from Flemisch et al. (2018): regular network with conductive fractures.
% DOI: 10.1016/j.advwatres.2017.10.036 (Section 4.1, Figures 8 and 11).
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
mdl.gravityOn = false;
mdl.condenseDiscontinuityPressure = true;

%% MESH

% Create mesh (N must be a multiple of 8 to retain the network junctions)
Lx = 1.0;  % Horizontal dimension
Ly = 1.0;  % Vertical dimension
N  = 24;   % Refinement parameter; 2*(N+1)^2 triangular elements
assert(N > 0 && mod(N,8) == 0, 'N must be a positive multiple of 8.');
[node, elem] = regularMesh(Lx, Ly, N+1, N+1, [], [], 'CST');

% Shift the interior grid by half a cell in both directions. Junctions lie
% on diagonal mesh edges, while the boundary strips have half-width cells.
column = round(node(:,1)*(N+1)/Lx);
row = round(node(:,2)*(N+1)/Ly);
coord = [0, ((0:N-1)+0.5)/N, 1];
node(:,1) = Lx*coord(column+1);
node(:,2) = Ly*coord(row+1);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluid (normalized benchmark coefficients, not SI water properties)
water = Fluid('water');
water.rho = 1.0;
water.mu  = 1.0;

% Create porous media
rock = PorousMedia('rock');
rock.K   = 1.0;
rock.phi = 0.25;

% Set materials to model
mdl.setMaterial(rock, water);

%% DISCONTINUITIES

% Fracture polylines from Figure 8
fractureData = { [0.000, 0.500; 1.000, 0.500];
                 [0.500, 0.000; 0.500, 1.000];
                 [0.500, 0.750; 1.000, 0.750];
                 [0.750, 0.500; 0.750, 1.000];
                 [0.500, 0.625; 0.750, 0.625];
                 [0.625, 0.500; 0.625, 0.750] };

% Physical aperture and equal tangential/normal fracture permeabilities
b = 1.0e-4;
kf = 1.0e4 * rock.K;

% Equivalent cubic-law aperture: w^3/(12*mu) = b*kf/mu.
% This computational aperture is distinct from the physical aperture b.
aperture = (12.0 * b * kf)^(1/3);

% Create discontinuities
nd = length(fractureData);
fractures(1,nd) = Discontinuity();
for i = 1:nd
    fractures(i) = Discontinuity(fractureData{i}, false);
    fractures(i).liquidFluid = water;
    fractures(i).initialAperture = aperture;
    % Transfer through half the physical aperture on each face
    fractures(i).leakoff = 2*kf/(water.mu*b);
end

% Add fractures to model
mdl.addPreExistingDiscontinuities(fractures);

%% BOUNDARY CONDITIONS

% Prescribed pressure on the right; top and bottom are impermeable
mdl.setPressureDirichletBCAtBorder('right', 1.0);

% Unit inward flux on the left (the paper's outward flux is -1).
% Neumann setters take nodal loads: integrate over the boundary edges.
leftNodes = mdl.getNodesAtBorder('left');
[y, order] = sort(mdl.NODE(leftNodes,2));
leftNodes = leftNodes(order);
dy = diff(y);
nodalFlow = 0.5 * ([dy; 0.0] + [0.0; dy]) * mdl.t;
for i = 1:length(leftNodes)
    mdl.setPressureNeumannBCAtNode(leftNodes(i), nodalFlow(i));
end

%% PROCESS

anl = Anl_Linear();

% Begin time measurement
caseTimer = tic;

% Run analysis
anl.run(mdl);

% Time measurement
elapsedTime = toc(caseTimer);
fprintf('Elapsed time: %.3f seconds\n', elapsedTime);

%% POST-PROCESS

% Plot model and embedded fractures
mdl.plotField('Model');
colorbar off; hold on;
for i = 1:nd
    fractures(i).plotIntersectedGeometry();
end

% Pressure field
mdl.plotField('Pressure');

% Pressure profile along y = 0.7
Xi = [0.0,0.7]; Xf = [1.0,0.7];
mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'x');

% Pressure profile along x = 0.5
Xi = [0.5,0.0]; Xf = [0.5,1.0];
mdl.plotFieldAlongSegment('Pressure', Xi, Xf, 500, 'x');