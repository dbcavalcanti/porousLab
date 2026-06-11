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

% Set model options
mdl.updatePorosity = true;

%% MESH

% Create mesh
Lx = 1.0;  % Horizontal dimension (m)
Ly = 1.0;  % Vertical dimension (m)
Nx = 1;    % Number of elements in the x-direction
Ny = 1;   % Number of elements in the y-direction
[node, elem] = regularMesh(Lx, Ly, Nx, Ny);

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create fluids
water = Fluid('water');

% Create porous media
rock = PorousMedia('rock');
rock.K     = 1.15741e-12;                   % Intrinsic permeability (m2)
rock.permeabilityModel = 'KozenyCarman';    % Permeability model
rock.phi   = 0.3;                           % Porosity
rock.Young = 1.0e+6;                        % Young modulus (Pa)
rock.nu    = 0.3;                           % Poisson ratio

% Set materials to model
mdl.setMaterial(rock, water);

%% BOUNDARY AND INITIAL CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('top',    [NaN, 0.0]);

% Pressure
Pleft = 1000.0; Pright = 0.0;  %(Pa)
mdl.setPressureDirichletBCAtBorder('left',  Pleft);
mdl.setPressureDirichletBCAtBorder('right', Pright);

%% EXPECTED SOLUTION - INCREMENTAL POROSITY UPDATE

% Pressure gradient
gradP = (Pleft - Pright) / Lx;

% Initial intrinsic permeability
K0 = rock.K;

% Initial porosity
phi0 = rock.phi;

% Initial Darcy flux
qx0 = K0 / water.mu * gradP;

% Kozeny-Carman permeability law
kozenyCarman = @(phi) K0 * ((phi / phi0)^3) * ((1.0 - phi0) / (1.0 - phi))^2;

% Prescribed vertical displacement magnitudes
uv = -[0.005; 0.01; 0.02; 0.05];

% Initial values
evOld  = 0.0;
phiOld = phi0;

% Preallocate
evRef  = zeros(length(uv),1);
phiRef = zeros(length(uv),1);
kRef   = zeros(length(uv),1);
qxRef  = zeros(length(uv),1);

for i = 1:length(uv)

    % Positive uv means imposed compression
    % Compression is negative in the strain convention
    evNew = uv(i) / Ly;

    % Volumetric strain increment
    dEv = evNew - evOld;

    % Incremental porosity update
    phiNew = 1.0 - (1.0 - phiOld) * exp(-dEv);

    % Updated permeability
    kNew = kozenyCarman(phiNew);

    % Updated Darcy flux
    qxNew = kNew / water.mu * gradP;

    % Store
    evRef(i)  = evNew;
    phiRef(i) = phiNew;
    kRef(i)   = kNew;
    qxRef(i)  = qxNew;

    % Update history variables
    evOld  = evNew;
    phiOld = phiNew;
end

ReferenceSolution = table(uv, evRef, phiRef, kRef, qxRef);

%% PROCESS

% Run analysis
anl = Anl_Transient("Newton");
anl.setUpTransientSolver(0.0, 1.0, 1.0);
anl.run(mdl);

% Initialize solution vector
dx = zeros(mdl.ndof);
XOld = mdl.U;

phiVector = zeros(length(uv)+1,1);
phiVector(1) = mdl.element(1).type.porosity;

kVector = zeros(length(uv)+1,1);
kVector(1) = mdl.element(1).type.intPoint(1).constitutiveMdl.intrinsicPermeability(phiVector(1));

for i = 1:length(uv)
    mdl.setDisplacementDirichletBCAtBorder('top', [NaN, uv(i)]);
    mdl.updateDirichletBC();
    
    % Compute model global matrices
    [A,b] = mdl.getLinearSystem(mdl.U,XOld,anl.nlscheme,1.0);
    XOld = mdl.U;

    mdl.updateStateVar();
    phiVector(1+i) = mdl.element(1).type.porosity;
    kVector(1+i) = mdl.element(1).type.intPoint(1).constitutiveMdl.intrinsicPermeability(phiVector(1+i));
end

%% POST-PROCESS

figure; hold on; title('Porosity');
plot(uv,phiRef,'-k','DisplayName','Analytical');
plot([0;uv],phiVector,'ro','DisplayName','Numerical');
xlabel('Volumetric strain');
ylabel('Porosity');

figure; hold on; title('Intrinsic permeability');
plot(phiRef, kRef,'-k','DisplayName','Analytical');
plot(phiVector,kVector,'ro','DisplayName','Numerical');
xlabel('Porosity');
ylabel('Intrinsic permeability');