%% ================ Terzaghi consolidation problem ========================
%
% Hydromechanical with single-phase flow validation problem
%
% Author: Danilo Cavalcanti
%
%% ========================================================================
%
% Initialize workspace
clear
initWorkspace; 
%
%% ========================================================================
% Run the case twice to compare the solution between:
%   Case 1: Linear elements for both pressure and displacement
%   Case 2: Quadratic displacement and linear pressure
% The goal is to see how the oscillations are reduced when analysing the
% case 2.

quadratic = [0, 1];
% quadratic = 0;

for i=1:size(quadratic,2)

    %% ============================== MESH  ===================================

    mdl = Model_HM();

    % --- Mesh of continuum elements ------------------------------------------

    % Mesh properties
    Lx = 1;     % Horizontal dimension (m)
    Ly = 30;     % Vertical dimension (m)
    Nx = 1;       % Number of elements in the x-direction
    Ny = 20;      % Number of elements in the y-direction

    % Generate the mesh
    [mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny);

    % Type of elements
    mdl.type = 'ISOQ4';

    % Different interpolation order?
    if quadratic(i) == true
        differentInterpConversion(mdl);
    end

    % Thickness (m)
    mdl.t = 1.0;

    %% ============================= MATERIAL =================================

    % Create the fluids
    water = Fluid('water',1000.0,1.0e-3,3.0e14);

    % Create the porous media
    rock = PorousMedia('rock');
    rock.Young = 2.5e7;         % Young modulus (Pa)
    rock.nu    = 0.2;           % Poisson ratio
    rock.rho   = 2000.0;        % Density (kg/m3)
    rock.phi   = 0.3;           % Porosity
    rock.Ks    = 1.5e17;        % Solid bulk modulus (Pa)
    rock.K     = 1e-14;         % Intrinsic permeability (m2)

    % Material parameters vector
    mdl.mat  = struct( ...
        'porousMedia',rock, ...
        'fluid',water);

    %% ======================= BOUNDARY CONDITIONS ============================
    % In case it is prescribed a pressure value different than zero, don't
    % forget also that you need to constraint these degrees of freedom.

    % Displacement boundary conditions
    CoordSupp  = [1 0 0 -1;
        1 0 Lx -1
        1 1 -1 0.0;];
    CoordLoad  = [];
    CoordPresc = [];

    % Define supports and loads
    [mdl.SUPP_u, mdl.LOAD_u, mdl.PRESCDISPL_u] = boundaryConditionsDisplacement(mdl.NODE, ...
        CoordSupp, CoordLoad, CoordPresc, Lx, Ly, Nx, Ny);

    % Apply pressure at the top (Pa)
    [mdl.LOAD_u] = pressureLoad(-1.0e4,[Lx, Ly],2,mdl.NODE,mdl.ELEM,mdl.LOAD_u);

    % Liquid pressure boundary conditions
    CoordSupp  = [1 -1 Ly]; % TODO. Check the condition
    CoordLoad  = [];
    CoordPresc = [];
    CoordInit  = [];

    % Define supports and loads
    [mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
        CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

    %% ===================== MODEL CONFIGURATION ==============================

    % Using Gauss quadrature
    mdl.intOrder = 2;

    %% ========================= INITIALIZATION ===============================

    % Perform the basic pre-computations associated to the model (dof
    % definition, etc.)
    mdl.preComputations();

    % Create the result object for the analysis
    ndPlot  = 3;
    dofPlot = 1; % 1 for X and 2 for Y
    result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

    %% ========================== RUN ANALYSIS ================================

    % Transient analysis parameters
    tinit = 0.0;          % Initial time
    dt    = 0.02;         % Time step
    tf    = 2;            % Final time
    dtmax = 0.1;          % Maximum time step
    dtmin = 0.001;        % Minimum time step

    % Solve the problem
    anl = Anl_Transient(result,"Newton");
    anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
    anl.process(mdl);

    %% ========================= CHECK THE RESULTS ============================

    % Plot pressure along a segment
    Xi  = [0.0 , 0.0];
    Xf  = [0.0 , Ly];
    npts = 500;
    mdl.plotPressureAlongSegment(Xi, Xf, npts,'y')
    % mdl.plotField('Pressure');

end

