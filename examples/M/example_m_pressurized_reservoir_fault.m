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
mdl.subDivIntegration = false;
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
K0 = 0.60;

% Create porous media
rock = PorousMedia('rock');   
rock.Young = 14.95e+9;        % Young modulus (Pa)
rock.nu    = K0/(1+K0);         % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% DISCONTINUITIES

% Fault center
Xdc = 0.5 * [Lx, Ly];

% Fault dip
dip = 70.0;
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
mdl.addPreExistingDiscontinuities(fault, discontinuityData);

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

% Parameters
tol    = 1.0e-5;
offset = 0.0;
DP     = 20.0e6;

region = [0.0 , 350.0;
          Xdc(1)+(350.0 - Xdc(2))*dxdy , 350.0;
          Xdc(1)+(650.0 - Xdc(2))*dxdy , 650.0;
          0.0 , 650.0;];

% Update the pressure at the left side
reservoir = inpoly(mdl.NODE,region);
P(reservoir == 1) = P0 + DP;

% Update pressure at the right side
% reservoir = isInsideRectangle(mdl.NODE, [0.5*Lx,350.0+offset-tol], [Lx+tol,650.0+offset+tol]);
% P(reservoir == 1) = P0 + DP;

% Set pressure jump at the discontinuities
nDiscontinuities = mdl.getNumberOfDiscontinuities();
for i = 1:nDiscontinuities
    nDiscontinuitySeg = mdl.discontinuitySet(i).getNumberOfDiscontinuitySegments();
    for j = 1:nDiscontinuitySeg
        Xr = mdl.discontinuitySet(i).segment(j).referencePoint();
        if (Xr(2) > 350.0) && (Xr(2) < 650.0)
            mdl.discontinuitySet(i).segment(j).DP = DP;
            econnect = mdl.ELEM{mdl.discontinuitySet(i).elemID(j)};
            P(econnect) = P0 + DP;
        end
    end
end

% Update values
mdl.setPorePressureField(P);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

mdl.plotField('PressureExt');
hold on;
fault.plotIntersectedGeometry();

% Plot stresses
mdl.plotField('Sx');
hold on;
fault.plotIntersectedGeometry();
mdl.plotField('Sy');
hold on;
fault.plotIntersectedGeometry();
mdl.plotField('Sxy');
hold on;
fault.plotIntersectedGeometry();

mdl.plotFieldAlongSegment('Sx',[0.5*Lx,0.0],[0.5*Lx,Ly],100,'y');
mdl.plotFieldAlongDiscontinuiy('St',1,'y');
mdl.plotFieldAlongDiscontinuiy('Sn',1,'y');


