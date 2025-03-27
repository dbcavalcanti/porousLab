%% DESCRIPTION
%
% Elastoplastic pressurized cylinder example.
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

mdl = Model_M();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
ri = 0.1;      % Internal radius of the cylinder
re = 0.2;      % External radius of the cylinder
Nx = 10;       % Number of elements in the x-direction
Ny = 10;       % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(1.0, 1.0, Nx, Ny);

% Transform to cylindrical coordinates
r = ri + node(:, 1) * (re - ri);
theta = node(:, 2) * (pi / 2);

% Set the nodes with the transformed coordinates
node = [r .* cos(theta), r .* sin(theta)];
mdl.setMesh(node,elem);

% --- Material properties of the domain -----------------------------------

% Create the porous media
rock = PorousMedia('rock');
rock.mechanical = 'vonMises'; % Elastoplastic with von Mises criteria 
rock.Young = 2.1e11;          % Young modulus (Pa)
rock.nu    = 0.3;             % Poisson ratio
rock.sy0   = 2.40e8;          % Initial yield stress (Pa)
rock.Kp    = 0.0;             % Plastic modulus (Pa)

% Material parameters vector
mdl.mat  = struct('porousMedia',rock);

% --- Boundary conditions -------------------------------------------------

mdl.setDisplacementDirichletBCAtBorder('left',[0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('bottom',[NaN, 0.0]);

% Internal pressure (Pa)
pint = 192.0905814164710e06;

% Compute the radius of each node wrt to the center at (0,0)
r  = sqrt((mdl.NODE(:,1)).^2 + (mdl.NODE(:,2)).^2);
sn = mdl.NODE(:,2) ./ r;
cs = mdl.NODE(:,1) ./ r;

% Loaded nodes
internalNodes = (r-ri)<eps;

% Number of loaded nodes
nInternalNodes = sum(internalNodes);

% Force magnitude
F0 = (pint*mdl.t*ri*pi/2.0)/(nInternalNodes-1)/2.0;

% Count occurrences of each node
nodeCount = histcounts(mdl.ELEM(:), 1:(size(mdl.NODE,1)+1))';

% Apply the internal pressure to the nodes located at the internal face
mdl.LOAD(internalNodes,:) = F0 * nodeCount(internalNodes) .* [cs(internalNodes) , sn(internalNodes)];

%% RUN ANALYSIS

anl = Anl_Nonlinear('ArcLengthCylControl',true,0.01,2.0,100,100,4,1.0e-5);
ndId = mdl.closestNodeToPoint([ri,0.0]);
anl.setPlotDof(ndId,1)
anl.process(mdl);

%% POS-PROCESSING

mdl.plotField('S1');
mdl.plotField('Sr');
mdl.plotField('Sx');