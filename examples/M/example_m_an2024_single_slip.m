%% DESCRIPTION
%
% An, Shovkun and Tchelepi (2024), Section 3.2.1
% DOI: 10.1016/j.geoen.2024.213200.
%
% Physics:
% * Mechanical (M), plane strain
%
%% MODEL

% Create model
mdl = Model_M();

% Set model options
mdl.isPlaneStress    = false;
mdl.condenseEnrDofs  = false;
mdl.useNodalEnrDofs  = false;
mdl.symmetricSDAEFEM = false;
mdl.subDivIntegration = false;
mdl.t = 1.0;          % Out-of-plane thickness (m)

%% MESH

% Create a graded Q4 mesh in a 200 m by 200 m square.
% Start with the opening-example mesh and add vertical lines at the tips.
% The fracture crosses element interiors; the mesh is not rotated.
Lx           = 200.0;   % Horizontal dimension (m)
Ly           = 200.0;   % Vertical dimension (m)
Lf           = 10.0;    % Fracture length (m)
nSegments    = 128;     % Divisions along the horizontal reference fracture
h            = Lf/nSegments;
growthFactor = 1.35;
alpha        = pi/6.0/2*3;  % Fracture inclination (rad)

[node,elem] = singleFractureMesh(Lx,Lf,nSegments,alpha,growthFactor);

% Set mesh to model
mdl.setMesh(node,elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.Young = 1.0e10;   % Young modulus (Pa)
rock.nu    = 0.25;     % Poisson ratio

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Dirichlet boundary conditions
mdl.setDisplacementDirichletBCAtPoint([-Lx/2.0, 0.0],[NaN, 0.0]);
mdl.setDisplacementDirichletBCAtPoint([ Lx/2.0, 0.0],[NaN, 0.0]);
mdl.setDisplacementDirichletBCAtPoint([0.0, -Ly/2.0],[0.0, NaN]);
mdl.setDisplacementDirichletBCAtPoint([0.0,  Ly/2.0],[0.0, NaN]);

% Horizontal compressive traction
sigma = 1.0e7;        % Applied stress magnitude (Pa)
mdl.addLoadAtBorder('left',  1,  sigma);
mdl.addLoadAtBorder('right', 1, -sigma);

%% DISCONTINUITIES

% Geometry
Xd = [-Lf/2.0*cos(alpha), -Lf/2.0*sin(alpha); ...
       Lf/2.0*cos(alpha),  Lf/2.0*sin(alpha)];
fracture = Discontinuity(Xd,false);

% Set fracture material properties
stiffness = 1.0e3*rock.Young/h;      % Finite contact penalty (Pa/m)
fracture.cohesiveLaw     = 'mohrCoulomb';
fracture.initialAperture = 0.0;      % Initial aperture (m)
fracture.normalStiffness = stiffness;
fracture.shearStiffness  = stiffness;
fracture.contactPenalization = [];  % Use the stiffness specified above
fracture.frictionAngle  = pi/6.0;    % Friction angle (rad)
fracture.dilationAngle  = 0.0;
fracture.cohesion       = 0.0;
fracture.tensionCutOff  = 0.0;

% Add fracture to model
mdl.addPreExistingDiscontinuities(fracture);

%% PROCESS

anl = Anl_NonlinearQuasiStatic('LoadControl');
anl.increment     = 0.25;
anl.max_increment = 0.25;
anl.max_lratio    = 1.0;
anl.max_step      = 4;
anl.max_iter      = 40;
anl.tol           = 1.0e-8;

% Begin time measurement
caseTimer = tic;

% Run analysis
anl.run(mdl);

% Time measurement
elapsedTime = toc(caseTimer);
fprintf('Elapsed time: %.3f seconds\n', elapsedTime);
fprintf('Total iterations: %d\n', anl.totalIterations);

%% POST-PROCESS

% Plot mesh and fracture
mdl.plotField('Model'); hold on;
fracture.plotIntersectedGeometry();

% Analytical tangential jump (Eq. 45), measured from the first tip.
amplitude = 4.0 * (1.0-rock.nu^2) * sigma / rock.Young;
amplitude = amplitude * sin(alpha)* ...
    (cos(alpha) - tan(fracture.frictionAngle) * sin(alpha));
s_anal = linspace(0.0,Lf,401);
jump_anal = amplitude*sqrt(max(0.0,(Lf/2.0)^2-(s_anal-Lf/2.0).^2));

% Plot tangential displacement jump
mdl.plotFieldAlongDiscontinuiy('Dt',1,'x');
hold on;
plot(s_anal,jump_anal,'--k','LineWidth',1.2);
% xlabel('Distance from first tip (m)');
ylabel('Tangential displacement jump (m)');
legend('PorousLab','Analytical: 10 m fracture','Location','south');
ylim('auto');

% RMSE at fixed, uniformly spaced interior points along the fracture.
% Keep the same sampling points for different meshes and inclinations.
nSamples = 10000;
s_sample = ((1:nSamples)'-0.5)*Lf/nSamples;
jump_num = NaN(nSamples,1);
direction = (Xd(2,:)-Xd(1,:))/Lf;
for i = 1:numel(fracture.segment)
    segment = fracture.segment(i);
    s_segment = (segment.node-Xd(1,:))*direction';
    ids = s_sample >= min(s_segment) & s_sample < max(s_segment);
    jump_num(ids) = segment.intPoint(1).strain(1);
end
jump_ref = amplitude*sqrt(max(0.0,(Lf/2.0)^2-(s_sample-Lf/2.0).^2));
slipRMSE = 1000.0*sqrt(mean((jump_num-jump_ref).^2));  % mm

% Discretization counts for independent convergence runs
fprintf('\nFracture inclination: %g degrees\n',alpha*180.0/pi);
fprintf('Continuum elements: %d; total DOFs: %d; free DOFs: %d\n', ...
    mdl.nelem,mdl.ndof,numel(mdl.doffree));
fprintf('Tangential displacement jump RMSE: %.6g mm\n',slipRMSE);

%% AUXILIARY FUNCTIONS

function [node,elem] = singleFractureMesh(Lx,Lf,nSegments,alpha,growthFactor)
    % Graded Q4 mesh with fracture-tip edges and boundary midpoint nodes.
    h = Lf/nSegments;
    halfDomain = Lx/2.0;

    % Start with a uniform mesh around the fracture.
    xcoord = linspace(-Lf/2.0, Lf/2.0, nSegments+1);
    ycoord = (-Lf/2.0-0.5*h):h:(Lf/2.0+0.5*h);

    % Put the inclined fracture tips on vertical element edges. Adding only
    % these x-coordinates retains the horizontal grid and avoids tip vertices.
    xcoord = unique([xcoord, 0.5*Lf*cos(alpha)*[-1.0,1.0]]);

    % Increase the element size towards each outer boundary.
    coordinates = {xcoord, ycoord};
    for i = 1:2
        coord = coordinates{i};
        positive = coord(end);
        negative = -coord(1);
        spacing = h;
        while positive(end) < halfDomain
            spacing = growthFactor*spacing;
            if halfDomain-positive(end) <= 1.6*spacing
                positive(end+1) = halfDomain;
            else
                positive(end+1) = positive(end)+spacing;
            end
        end
        spacing = h;
        while negative(end) < halfDomain
            spacing = growthFactor*spacing;
            if halfDomain-negative(end) <= 1.6*spacing
                negative(end+1) = halfDomain;
            else
                negative(end+1) = negative(end)+spacing;
            end
        end
        coordinates{i} = [-fliplr(negative(2:end)), coord, positive(2:end)];
    end
    xcoord = coordinates{1};
    ycoord = coordinates{2};

    % Create nodes and counterclockwise Q4 connectivities.
    [X, Y] = meshgrid(xcoord, ycoord);
    node = [X(:), Y(:)];
    nx = numel(xcoord);
    ny = numel(ycoord);
    elem = cell((nx-1)*(ny-1),1);
    k = 0;
    for ix = 1:nx-1
        for iy = 1:ny-1
            k = k+1;
            n1 = (ix-1)*ny+iy;
            n2 = ix*ny+iy;
            elem{k} = [n1, n2, n2+1, n1+1];
        end
    end

    % Put nodes exactly at the four boundary midpoints for the roller supports.
    % The smooth adjustment acts only outside 60 m from the domain centre.
    [~, ixMid] = min(abs(xcoord));
    [~, iyMid] = min(abs(ycoord));
    xMid = xcoord(ixMid);
    yMid = ycoord(iyMid);
    r = hypot(node(:,1),node(:,2));
    blend = min(max((r-60.0)/(95.0-60.0),0.0),1.0);
    blend = blend.^2.*(3.0-2.0*blend);
    node(:,1) = node(:,1)-xMid*blend.* ...
        (halfDomain^2-node(:,1).^2)/(halfDomain^2-xMid^2);
    node(:,2) = node(:,2)-yMid*blend.* ...
        (halfDomain^2-node(:,2).^2)/(halfDomain^2-yMid^2);
end
