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

% Load nodes and elements
load('NODE_2');
load('ELEM_2');

% Set mesh to model
mdl.setMesh(NODE, ELEM);

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

% Analysis parameters
adapt_incr = true;    % Increment size adjustment
increment  = 0.01;     % Initial increment of load ratio
max_lratio = 10.0;     % Limit value of load ratio
max_step   = 15;      % Maximum number of steps
max_iter   = 100;     % Maximum number of iterations in each step
trg_iter   = 4;       % Desired number of iterations in each step
tol        = 1.0e-5;  % Numerical tolerance for convergence

% Run analysis
anl = Anl_Nonlinear('LoadControl', adapt_incr, increment, max_lratio, max_step, max_iter, trg_iter, tol);
anl.run(mdl);

%% POST-PROCESS

% Plot contours
mdl.plotField('Sxy');
