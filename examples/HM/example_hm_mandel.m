%% DESCRIPTION
%
% Mandel consolidation problem.
%
% Reference:
% * Keilegavlen et al. (2021). PorePy: an open-source software for
%   simulation of multiphysics processes in fractured porous media.
% * Mandel (1953). Consolidation des sols. Géotechnique, 3(7):287–299.
% * Cheng and Detournay (1988). A direct boundary element method for plane
%   strain poroelasticity. Int J Numer Anal Methods Geomech, 12:551–572.
% * Mikelić et al. (2014). Numerical convergence study of iterative
%   coupling for coupled flow and geomechanics. Comput Geosci, 18:325–341.
% * Walker et al. (2023). Multiphysics modelling in PyLith: poroelasticity.
%   Geophys J Int, 235:2442–2475.
% * GEOS Mandel benchmark documentation. 
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

% Problem dimensions (m)
a = 1.0;
b = 1.0;

% Create mesh
[node, elem] = regularMesh(a, b, 30, 30);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');
water.K = 2.2727e+09;

% Bulk modulus and Shear modulus (Pa)
Kb = 6.6667e7;
G  = 4.0e7;

% Create porous media
rock = PorousMedia('rock');
rock.K     = 1.0e-12;                        % Intrinsic permeability (m2)
rock.phi   = 0.375;                          % Porosity
rock.Young = 9*Kb*G/(3*Kb + G);              % Young modulus (Pa)
rock.nu    = (3*Kb - 2*G)/(2*(3*Kb + G));    % Poisson ratio

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

% Mandel load parameter
F = 1.0e4;

% Load the at top
mdl.addLoadAtBorder('top',2,-F/a);

% Set master-slave displacement to impose the rigid plate condition
nodeTopLeft = mdl.closestNodeToPoint([0,b]);
nodesTop = setdiff(mdl.getNodesAtBorder("top"), nodeTopLeft);
mdl.setMasterSlaveDisplacements(nodeTopLeft, nodesTop, 2)

% Pressure
mdl.setPressureDirichletBCAtBorder('right', 0.0);

%% PROCESS

anl = Anl_Transient("Newton");

% Analysis parameters
ti        = 0.0;         % Initial time
dt        = 0.0001;      % Time step
tf        = 0.5;         % Final time
dtmax     = 0.1;         % Maximum time step
dtmin     = 0.00001;     % Minimum time step
adaptStep = true;        % Adaptive step size
anl.setUpTransientSolver(ti, dt, tf, dtmax, dtmin, adaptStep);

% Run analysis
anl.run(mdl);

%% ANALYTICAL SOLUTION

% Number of expansion terms
nRoots = 200;

% Spatial discretization
x = 0:0.1:a;     % Horizontal coordinate in the quarter domain
y = 0:0.1:b;     % Vertical coordinate in the quarter domain
X = x/a;         % Normalized horizontal coordinate
Y = y/b;         % Normalized vertical coordinate

% Biot modulus (Pa)
M = 1.0 / (rock.phi/water.K + (rock.biot - rock.phi)/rock.Ks); 

% Undrained bulk modulus (Pa)
Ku = Kb + rock.biot * rock.biot * M; 

% Skempton coefficient (-)
B = rock.biot*M / Ku;             

% Undrained Poisson's ratio (-)
nu_u = (3*Ku - 2*G) / (2*(3*Ku + G));       

% Hydraulic coefficient (m^2/(Pa*s))
kh = rock.K/water.mu;

% Diffusivity coefficient (m^2/s)
c = 2*kh*B*B*G*(1 - rock.nu)*(1 + nu_u)^2 / ...
    (9*(1 - nu_u)*(nu_u - rock.nu)); 

% INITIAL UNDRAINED RESPONSE ----------------------------------------------

% Initial pore pressure
p0 = B * (1 + nu_u)*F/(3*a);

% Initial horizontal displacement along y = 0
ux0 = F * nu_u/(2*G) * (x/a);

% Initial vertical displacement along x = 0
uy0 = -F*(1 - nu_u)/(2*G) * (y/a);

% Reference values for normalization
ux0_right = F * nu_u/(2*G);            % u_x at x = a and t = 0+
uy0_top   = -F * b * (1 - nu_u)/(2*G*a); % u_y at y = b and t = 0+

% ROOTS OF THE MANDEL EQUATION --------------------------------------------

% Equation: tan(alpha_i) = ((1 - nu)/(nu_u - nu)) * alpha_i

A = (1 - rock.nu)/(nu_u - rock.nu);
alpha_i = zeros(nRoots,1);

for i = 1:nRoots
    % Limits of the search. The tangent function has vertical asymptotes
    alphaLeft  = (i - 1.0)*pi + 1.0e-10;
    alphaRight = (i - 0.5)*pi - 1.0e-10;
    
    % Objective function
    fRoot = @(alpha) tan(alpha) - A*alpha;
    
    % Solution of the equation
    alpha_i(i) = fzero(fRoot, [alphaLeft, alphaRight]);
end

% MANDEL SOLUTION ---------------------------------------------------------

% Dimensionless time
T = c*tf/(a^2);

% Auxiliary sums
sumPressure = zeros(size(x));
sumUx       = zeros(size(x));
sumUy       = 0.0;

% Compute the series term
for j = 1:nRoots
    ai = alpha_i(j);

    den     = ai - sin(ai)*cos(ai);
    expTerm = exp(-ai*ai*T);

    coefP  = sin(ai)/den;
    coefUx = cos(ai)/den;
    coefUy = sin(ai)*cos(ai)/den;

    % Pressure auxiliary sum
    sumPressure = sumPressure + coefP * (cos(ai*X) - cos(ai)) * expTerm;

    % Horizontal displacement auxiliary sum
    sumUx = sumUx + coefUx * sin(ai*X) * expTerm;

    % Vertical displacement auxiliary sum
    sumUy = sumUy + coefUy * expTerm;
end

% Solution
p  = 2*p0*sumPressure;
ux = F/G*sumUx + (F*rock.nu/(2*G*a) - F*nu_u/(G*a)*sumUy)*x;
uy = (-F*(1 - rock.nu)/(2*G*a) + F*(1 - nu_u)/(G*a)*sumUy)*y;

%% POST-PROCESS

% Pore pressure
mdl.plotFieldAlongSegment('Pressure', [0.0, b], [a, b], 500, 'x');
plot(X, p,'xr', 'LineWidth', 2, 'DisplayName', 'Analytical');
legend('PorousLab', 'Analytical');

% Horizontal displacement
mdl.plotFieldAlongSegment('Ux', [0.0, b], [a, b], 500, 'x');
plot(X, ux, 'xr', 'LineWidth', 2, 'DisplayName', 'Analytical');
legend('PorousLab', 'Analytical');

% Vertical displacement
mdl.plotFieldAlongSegment('Uy', [0.0, 0.0], [0.0, b], 500, 'x');
plot(X, uy, 'xr', 'LineWidth', 2, 'DisplayName', 'Analytical');
legend('PorousLab', 'Analytical');

% Plot field
mdl.plotField('Pressure');
mdl.plotField('Ux');