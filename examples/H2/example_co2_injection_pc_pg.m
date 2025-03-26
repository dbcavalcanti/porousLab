%% ================ Two-Phase flow in porous media ====================
%
%% Reference:
%
% Benisch, K., Graupner, B., & Bauer, S. (2013). 
% The coupled OpenGeoSys-eclipse simulator for simulation of CO2
% storage–code comparison for fluid flow and geomechanical processes.
% Energy Procedia, 37, 3663-3671.
%
%% References to similar problems
%
% Graupner, B. J., Li, D., & Bauer, S. (2011). The coupled simulator 
% ECLIPSE–OpenGeoSys for the simulation of CO2 storage in saline 
% formations. Energy Procedia, 4, 3794-3800.
%
% Class, H., Ebigbo, A., Helmig, R. et al.
% A benchmark study on problems related to CO2 storage in geologic
% formations. Comput Geosci 13, 409–434 (2009).
% https://doi.org/10.1007/s10596-009-9146-x
%
% Vilarrasa, Víctor. "Thermo-hydro-mechanical impacts of carbon dioxide 
% (CO2) injection in deep saline aquifers." (2012).
%
% Zhou, Q., J. T. Birkholzer, C.-F. Tsang, and J. Rutqvist (2008), 
% A method for quick assessment of CO 2 storage capacity in closed
% and semi-closed saline formations. Int. J. Greenh. Gas Control, 2(4),
% 626–639, doi:10.1016/j.ijggc.2008.02.004
% 
% Author: Danilo Cavalcanti
%
%% ========================================================================
%
% Initialize workspace
clear
initWorkspace; 
%
%% ========================== MODEL CREATION ==============================

mdl = Model_H2_PcPg();

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 200.0;      % Horizontal dimension (m)
Ly = 6.0;       % Vertical dimension (m)
Nx = 30;         % Number of elements in the x-direction
Ny = 60;          % Number of elements in the y-direction

% Generate the mesh
[mdl.NODE,mdl.ELEM] = regularMeshY(Lx, Ly, Nx, Ny,[],[],'ISOQ4',true,false);

% Type of elements
mdl.type = 'ISOQ4';
mdl.isAxisSymmetric = true;

% Thickness (m)
mdl.t = 1.0;

% --- Material properties of the domain -----------------------------------

% Create the fluids
brine = Fluid('brine',1173.0,1.252e-3,1.0e25);
co2   = Fluid('co2'  ,848.0,8.1e-5,1.0e25);

% Porous media properties
% --------------------------   |   K(m2) | phi | biot |  Ks   | Slr | Sgr |   Pb  | lambda | LiqRelPerm |  GasRelPerm  |  capPressure
aquifer = PorousMedia('aquifer', 3.0e-12 , 0.26 , 1.0 , 1.0e25 , 0.35 , 0.0 , 1.0e4 , 2.0 , 'BrooksCorey', 'BrooksCorey','BrooksCorey');
aquifer.setMinLiquidRelPermeability(1.0e-9);
aquifer.setMinGasRelPermeability(1.0e-9);

% Activate gravity
aquifer.gravityOn = true;

% Material parameters vector
% Same material for all elements
mdl.mat  = struct( ...
    'porousMedia',aquifer, ...
    'liquidFluid',brine,...
    'gasFluid',co2);

% --- Boundary conditions -------------------------------------------------
% In case it is prescribed a pressure value different than zero, don't 
% forget also that you need to constraint these degrees of freedom.

% Capillary pressure boundary conditions
CoordSupp  = [1 Lx -1];    % Fix the pressure at the right border         
CoordLoad  = []; 
CoordPresc = [1.0e4 Lx -1];
CoordInit  = [];
           
% Define supports and loads
[mdl.SUPP_p, mdl.LOAD_p, mdl.PRESCDISPL_p, mdl.INITCOND_p] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

% The aquifer is initially saturated with brine. Then, the capillary
% pressure is set to be equal to the gas entry pressure
mdl.INITCOND_p = 1.0e4*ones(size(mdl.INITCOND_p,1),1);

% Gas pressure boundary conditions
day = 60 * 60 * 24;
qinj = (1.0 / day) * Ly / (Ny + 1);
CoordSupp  = [1 Lx -1];      % Fix the pressure at the right border        
CoordLoad  = [qinj 0 -1];                     
CoordPresc = [];         
CoordInit  = [];
           
% Define supports and loads
[mdl.SUPP_pg, mdl.LOAD_pg, mdl.PRESCDISPL_pg, mdl.INITCOND_pg] = boundaryConditionsPressure(mdl.NODE, ...
    CoordSupp, CoordLoad, CoordPresc, CoordInit, Lx, Ly, Nx, Ny);

% Set hydrostatic pressure 
depth = 1500.0;          % Depth of the aquifer (top border)
depth0 = depth + Ly;     % Depth of the aquifer (bottom border)
grav   = 9.81;           % Gravity acceleration (m/s2)
for i = 1:size(mdl.NODE,1)
    h = depth0 - mdl.NODE(i,2);
    pl = grav * brine.rho * h;
    pc = 1.0e4;
    mdl.INITCOND_pg(i) = pc + pl;
    if (mdl.SUPP_pg(i) == 1)
        mdl.PRESCDISPL_pg(i) = pc + pl;
    end
end

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
mdl.intOrder = 2;

% Diagonalize compressibility matrix
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% ========================= INITIALIZATION ===============================

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot),[],[],[]);

%% ========================== RUN ANALYSIS ================================

% Transient analysis parameters
tinit = 0.05*day;   % Initial time
dt    = 0.05*day;   % Time step
tf    = 41*day;  % Final time
dtmax = 0.1*day;
dtmin = 1.0e-3*day;

% Solve the problem
anl = Anl_Transient(result,"Picard");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.setRelativeConvergenceCriteria(true);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [Lx , 0.0];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'x')
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');
