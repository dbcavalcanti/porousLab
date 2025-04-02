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
%% INITIALIZATION
close all; clear; clc;

% Path to source directory
src_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src');
addpath(genpath(src_dir));

print_header;

%% MODEL

% Create model
mdl = Model_M();

%% MESH

% Load the mesh
load('MeshStripFooting');
[node, elem] = convertToQuadraticMesh(node, elem);
mdl.re

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical    = 'druckerPrager';  % Mechanical constitutive law
rock.Young         = 1.0e+7;           % Young modulus (kPa)
rock.nu            = 0.48;             % Poisson ratio
rock.cohesion      = 490.0;            % Cohesion (kPa)
rock.frictionAngle = 20*pi/180;        % Friction angle (rad)
rock.dilationAngle = 20*pi/180;        % Dilation angle (rad)
rock.MCmatch       = "planestrain";

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

for i = 1:mdl.nnodes
    if ((mdl.NODE(i,1) <= 0.2) && (mdl.NODE(i,2) == 5.0))
        mdl.addLoadAtNode(i, [0.0 , -1.0e2])
    end
end

%% PROCESS

% Configure analysis
anl = Anl_NonlinearQuasiStatic();
anl.method     = 'GeneralizedDisplacement';
anl.adjustStep = true;
anl.increment  = 0.01;
anl.max_increment  = 1.0e3;
anl.max_lratio = 10.0;
anl.max_step   = 59;
anl.max_iter   = 100;
anl.trg_iter   = 4;

% Node and DOF used to plot Load Factor vs Displacement
ndId = mdl.closestNodeToPoint([0.5, 5.0]);
anl.setPlotDof(ndId, 2);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

anl.plotCurves();
% Plot contours
mdl.plotField('PEMAG');
