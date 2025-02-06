%% ========================================================================
%
% Block crossed by a strong discontinuity
%
% Author: Danilo Cavalcanti
%
%% ========================================================================
%
% Initialize workspace
clear
initWorkspace; 
%
%% ============================== MESH  ===================================

mdl = Model_H();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 10.0;     % Horizontal dimension (m)
Ly = 3.0;     % Vertical dimension (m)
Nx = 51;       % Number of elements in the x-direction
Ny = 31;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

letters = {
    'P', [0 0; 0 1; 0.5 1; 0.5 0.5]
    'O', [0 0; 0 1; 1 1; 1 0; 0 0];
    'R', [0 0; 0 1; 0.5 1; 0.5 0.5; 0 0.5; 0.5 0];
    'O', [0 0; 0 1; 1 1; 1 0; 0 0];
    'U', [0 1; 0 0; 1 0; 1 1];
    'S', [1 1; 0 1; 0 0.5; 1 0.5; 1 0; 0 0];
    'L', [0 1; 0 0; 1 0];
    'A', [0 0; 0.5 1; 1 0; 0.25 0.5; 0.75 0.5];
    'B', [0 0; 0 1; 0.5 1; 0.75 0.75; 0.5 0.5; 0 0.5; 0.5 0.5; 0.75 0.25; 0.5 0; 0 0]
};

% Compute total width for normalization
num_letters = length(letters);
max_letter_width = 1.2;  % Max relative width per letter (including spacing)
total_width = num_letters * max_letter_width;

% Compute scaling factors
scale_x = Lx / 1.1 /total_width; % Fit within Lx
scale_y = Ly * 0.9;  % Scale to 90% of Ly to ensure margin

% Initialize storage for final coordinates
final_coords = cell(num_letters, 1);

% Transform and store final coordinates
x_offset = 0.5; % Initial offset to center text
y_offset = 0.1; % Initial offset to center text
for i = 1:num_letters
    data = letters{i,2};
    data(:,1) = data(:,1) * scale_x + x_offset;
    data(:,2) = data(:,2) * scale_y + y_offset; % Scale to domain height
    final_coords{i} = data; % Store in array
    x_offset = x_offset + max_letter_width * scale_x; % Move to next letter position
end

% Create the discontinuity 
% Create discontinuities for each individual line segment
fractures = [];

for i = 1:num_letters
    letter_coords = final_coords{i};
    for j = 1:size(letter_coords, 1) - 1
        % Extract segment endpoints
        segment = letter_coords(j:j+1, :);
        % Create a fracture for this segment
        fractures = [fractures, Discontinuity(segment, true)];
    end
    fractures = [fractures, Discontinuity([0.0,0.2;Lx,0.5], true)];
end

%% ============================= MATERIAL =================================

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,2.0e9);

% Create the porous media
rock = PorousMedia('rock');
rock.K     = 9.8e-16;        % Intrinsic permeability (m2)
rock.phi   = 0.25;          % Porosity

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluid',water);

% Set the fracture material properties
for i = 1:length(fractures)
    fractures(i).initialAperture = 1.0;
    fractures(i).fluid = water;
end

%% ======================= BOUNDARY CONDITIONS ============================

% Pore pressure boundary conditions
CoordSupp  = [1 0.0 -1;1 Lx -1];         
CoordLoad  = [];            
CoordPresc = [20.0 0.0 -1;0.0 Lx -1];            
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
tinit = 1.0;   % Initial time
dt    = 1.0;   % Time step
tf    = 400;  % Final time

% Solve the problem
% anl = Anl_Transient(result,"Newton");
% anl.setUpTransientSolver(tinit,dt,tf,100.0,0.001,true);
% anl.process(mdl);

anl = Anl_Linear(result);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
mdl.printResults();

mdl.plotField('Pressure'); hold on
for i = 1:length(fractures)
    fractures(i).plotIntersectedGeometry();
end