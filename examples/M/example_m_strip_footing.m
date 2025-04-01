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

% Set model options
mdl.gravityOn = true;

%% MESH

% Load the mesh
load('MeshStripFooting');

% Set mesh to model
mdl.setMesh(node, elem);

%% MATERIALS

% Create porous media
rock = PorousMedia('rock');
rock.mechanical    = 'druckerPrager';  % Mechanical constitutive law
rock.rho           = 2.0e+3;           % Density (kg/m3)
rock.Young         = 2.0e+7;           % Young modulus (Pa)
rock.nu            = 0.49;             % Poisson ratio
rock.cohesion      = 5.0e+4;           % Cohesion (Pa)
rock.frictionAngle = 20*pi/180;        % Friction angle (rad)
rock.dilationAngle = 20*pi/180;        % Dilation angle (rad)

% Set materials to model
mdl.setMaterial(rock);

%% BOUNDARY CONDITIONS

% Displacements
mdl.setDisplacementDirichletBCAtBorder('left',   [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom', [NaN, 0.0]);

for i = 1:mdl.nnodes
    if ((mdl.NODE(i,1) <= 0.2) && (mdl.NODE(i,2) == 5.0))
        mdl.addLoadAtNode(i, [0.0 , -1.0e3])
    end
end

%% PROCESS

% Configure analysis
anl = Anl_Nonlinear();
anl.method     = 'ArcLengthCylControl';
anl.adjustStep = true;
anl.increment  = 0.01;
anl.max_lratio = 10.0;
anl.max_step   = 50;
anl.max_iter   = 100;
anl.trg_iter   = 4;

% Node and DOF used to plot Load Factor vs Displacement
ndId = mdl.closestNodeToPoint([0.5, 5.0]);
anl.setPlotDof(ndId, 2);

% Run analysis
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('PEMAG');
