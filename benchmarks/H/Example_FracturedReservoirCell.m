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
Lx = 200.0;     % Horizontal dimension (m)
Ly = 200.0;     % Vertical dimension (m)
Nx = 53;       % Number of elements in the x-direction
Ny = 53;       % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

% Create the discontinuity 
fractures = [];
nFrac = 10; rng(123);
count = 0;
% while count < nFrac
%     x1 = rand() * Lx;
%     y1 = rand() * Ly;
%     x2 = rand() * Lx;
%     y2 = rand() * Ly;
%     l = sqrt((x2-x1)^2+(y1-y2)^2);
%     if l > Lx/3.0
%         fractures = [fractures, Discontinuity([x1,y1;x2,y2], true)];
%         count = count + 1;
%     end
% end

fractures = [fractures, Discontinuity([93.0402, 150.0; 59.6560, 99.6233], true)];
fractures = [fractures, Discontinuity([92.6980, 190.4503; 67.1685, 92.4745], true)];
fractures = [fractures, Discontinuity([22.4000, 128.8991; 0.0000, 87.0630], true)];
fractures = [fractures, Discontinuity([26.6108, 178.4629; 0.0000, 129.0134], true)];
fractures = [fractures, Discontinuity([75.6925, 187.3910; 54.7212, 176.7641], true)];
fractures = [fractures, Discontinuity([92.2234, 176.9258; 72.2743, 153.8060], true)];
fractures = [fractures, Discontinuity([141.0454, 164.0137; 66.8328, 79.8375], true)];
fractures = [fractures, Discontinuity([128.4986, 184.7592; 91.2519, 96.4418], true)];
fractures = [fractures, Discontinuity([51.3924, 166.1533; 77.8526, 23.2657], true)];
fractures = [fractures, Discontinuity([31.8595, 114.4861; 48.3436, 32.3254], true)];
fractures = [fractures, Discontinuity([17.5208, 34.1698; 18.7207, 37.5684], true)];
fractures = [fractures, Discontinuity([31.8171, 22.0021; 28.4015, 67.5715], true)];
fractures = [fractures, Discontinuity([34.0408, 183.4709; 66.6623, 87.8132], true)];
fractures = [fractures, Discontinuity([20.7393, 107.4604; 3.4102, 154.1464], true)];
fractures = [fractures, Discontinuity([134.4233, 171.7966; 185.1331, 72.5896], true)];
fractures = [fractures, Discontinuity([75.6423, 195.1802; 118.4192, 97.5307], true)];
fractures = [fractures, Discontinuity([106.8604, 188.9048; 103.0021, 200.0000], true)];
fractures = [fractures, Discontinuity([151.5427, 118.7506; 78.9161, 100.1690], true)];
fractures = [fractures, Discontinuity([131.6637, 173.4711; 160.5112, 152.6563], true)];
fractures = [fractures, Discontinuity([159.8756, 161.0135; 169.4091, 147.8234], true)];
fractures = [fractures, Discontinuity([85.2698, 200.0000; 23.1479, 129.8266], true)];
fractures = [fractures, Discontinuity([25.9630, 115.9241; 0.0000, 81.7746], true)];
fractures = [fractures, Discontinuity([84.3391, 197.1414; 144.4823, 133.8286], true)];


%% ============================= MATERIAL =================================

% Create the fluids
water = Fluid('water',1000.0,1.0e-3,2.2e9);

% Create the porous media
rock = PorousMedia('rock');
rock.K     = 1.0e-14;        % Intrinsic permeability (m2)
rock.phi   = 0.25;          % Porosity

% Material parameters vector
mdl.mat  = struct( ...
    'porousMedia',rock, ...
    'fluid',water);

% Set the fracture material properties
% Set the fracture material properties
for i = 1:length(fractures)
    fractures(i).initialAperture = 1.0e-3;
    fractures(i).fluid = water;
end

%% ======================= BOUNDARY CONDITIONS ============================

% Pore pressure boundary conditions
CoordSupp  = [1 -1 0.0;1 -1 Ly];         
CoordLoad  = [];            
CoordPresc = [0.0 -1 0.0;1000000.0 -1 Ly];            
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
tf    = 127.0;  % Final time

% Solve the problem
anl = Anl_Transient(result,"Newton");
anl.setUpTransientSolver(tinit,dt,tf,500.0,0.001,true);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Plot pressure along a segment
mdl.printResults();

mdl.plotField('Pressure'); 
mdl.plotField('Model'); colorbar off, hold on
for i = 1:length(fractures)
    fractures(i).plotIntersectedGeometry();
end