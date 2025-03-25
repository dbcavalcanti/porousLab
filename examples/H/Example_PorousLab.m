%% DESCRIPTION
%
% Block crossed by a strong discontinuity.
%
% Physics:
% * Single-phase hydraulic (H)
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

% Create model
mdl = Model_H();

%% MESH GENERATION

% Mesh generation
Lx = 10.0;  % Horizontal dimension (m)
Ly = 3.0;   % Vertical dimension (m)
Nx = 63;    % Number of elements in the x-direction
Ny = 33;    % Number of elements in the y-direction

[mdl.NODE,mdl.ELEM] = regularMeshY(Lx,Ly,Nx,Ny);

% Obtain coordinates of discontinuty segments
letters = {
    'P', [0 0; 0 1; 0.5 1; 0.5 0.5];
    'O', [0 0; 0 1; 1 1; 1 0; 0 0];
    'R', [0 0; 0 1; 0.5 1; 0.5 0.5; 0 0.5; 0.5 0];
    'O', [0 0; 0 1; 1 1; 1 0; 0 0];
    'U', [0 1; 0 0; 1 0; 1 1];
    'S', [1 1; 0 1; 0 0.5; 1 0.5; 1 0; 0 0];
    'L', [0 1; 0 0; 1 0];
    'A', [0 0; 0.5 1; 1 0];
    'B', [0 0; 0 1; 0.5 1; 0.75 0.75; 0.5 0.5; 0 0.5; 0.5 0.5; 0.75 0.25; 0.5 0; 0 0]
};

% Compute total width for normalization
num_letters = length(letters);
max_letter_width = 1.5;  % Max relative width per letter (including spacing)
total_width = num_letters * max_letter_width;

% Compute scaling factors
scale_x = Lx / 1.1 /total_width; % Fit within Lx
scale_y = Ly * 0.9; % Scale to 0.9 of Ly to ensure margin

% Transform and store final coordinates
final_coords = cell(num_letters, 1);
x_offset = 0.7;  % Initial offset to center the text
y_offset = 0.1;  % Initial offset to center the text

for i = 1:num_letters
    data = letters{i,2};
    data(:,1) = data(:,1) * scale_x + x_offset;
    data(:,2) = data(:,2) * scale_y + y_offset;  % Scale to domain height
    final_coords{i} = data;  % Store in array
    x_offset = x_offset + max_letter_width * scale_x;  % Move to next letter position
end

% Discontinuity generation
fractures = [];

for i = 1:num_letters
    letter_coords = final_coords{i};
    for j = 1:size(letter_coords, 1) - 1
        % Extract segment endpoints
        segment = letter_coords(j:j+1, :);
        
        % Create a fracture for this segment
        fractures = [fractures, Discontinuity(segment, true)];
    end
    fractures = [fractures, Discontinuity([0.0, 0.5*Ly; Lx, 0.5*Ly], true)];
end

%% MATERIAL CREATION

% Create fluids
water = Fluid('water',1000.0,1.0e-3,2.2e9);

% Create the porous media
rock = PorousMedia('rock');
rock.K     = 1.0e-16;        % Intrinsic permeability (m2)
rock.phi   = 0.25;          % Porosity

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluid',water);

% Set the fracture material properties
for i = 1:length(fractures)
    fractures(i).initialAperture = 1.0e-1;
    fractures(i).fluid = water;
end
fractures(end).initialAperture = 1.0e-1;

%% ======================= BOUNDARY CONDITIONS ============================

% Pore pressure boundary conditions
CoordSupp  = [1 0.0 -1;1 Lx -1];         
CoordLoad  = [];            
CoordPresc = [1000000.0 0.0 -1;0.0 Lx -1];            
CoordInit  = [];                   
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

%% ===================== MODEL CONFIGURATION ==============================

% Using Gauss quadrature
mdl.intOrder = 2;

%% ========================= INITIALIZATION ===============================

% Create the discontinuity elements
% Intersect all fractures with the mesh
for i = 1:length(fractures)
    fractures(i).intersectMesh(mdl);
end

% Add the fracture to the model
mdl.addPreExistingDiscontinuities(fractures);

% Perform the basic pre-computations associated to the model
mdl.preComputations();

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Transient analysis parameters
tinit = 0.1;   % Initial time
dt    = 0.1;   % Time step
tf    = 5.5;  % Final time

% Solve the problem
anl = Anl_Transient(result,"Newton");
anl.setUpTransientSolver(tinit,dt,tf,1.0,0.001,true);
anl.process(mdl);

% anl = Anl_Linear(result);
% anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
% mdl.printResults();

mdl.plotField('Pressure'); hold on
for i = 1:length(fractures)
    fractures(i).plotIntersectedGeometry();
end