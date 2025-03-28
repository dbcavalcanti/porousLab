%% DESCRIPTION
%
% Liakopoulos problem using the Pc-Pg two-phase flow formulation
%
% References:
%
% Schrefler, B. A., and Z. Xiaoyong (1993).
% A fully coupled model for water flow and airflow in deformable
% porous media, Water Resour. Res., 29(1), 155–167,
% doi:10.1029/92WR01737.
%
% Schrefler, B. A., & Scotta, R. (2001).
% A fully coupled dynamic model for two-phase fluid flow in deformable
% porous media. Computer methods in applied mechanics and engineering,
% 190(24-25), 3223-3246.
%
% Physics:
% * Hydromechanical with two-phase flow (H2M)
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

mdl = Model_H2M();

%% MODEL CREATION

% --- Mesh of continuum elements ------------------------------------------

% Mesh properties
Lx = 0.1;     % Horizontal dimension (m)
Ly = 1.0;     % Vertical dimension (m)
Nx = 4;       % Number of elements in the x-direction
Ny = 40;      % Number of elements in the y-direction

% Generate the mesh
[node,elem] = regularMesh(Lx, Ly, Nx, Ny);
mdl.setMesh(node,elem);

% --- Material properties of the domain -----------------------------------

% Create the fluids
water   = Fluid('water');
water.K = 2.0e9;
gas     = Fluid('gas');
gas.rho = 1.20;
gas.mu  = 1.8e-5;
gas.K   = 1.0e5;

% Create the porous media
rock = PorousMedia('rock');
rock.K     = 4.5e-13;   % Intrinsic permeability (m2)
rock.phi   = 0.2975;           % Porosity
rock.Ks    = 1.0e12;
rock.Sgr   = 0.2;
rock.lambda = 3.0;
rock.Young = 1.3e6;         % Young modulus (Pa)
rock.nu    = 0.4;           % Poisson ratio
rock.rho = 2000.0;
rock.liqRelPermeability = 'Liakopoulos';
rock.gasRelPermeability = 'BrooksCorey';
rock.capillaryPressure  = 'Liakopoulos';
rock.setMinLiquidRelPermeability(1.0e-4);
rock.setMinGasRelPermeability(1.0e-4);

% Set the material to the model
mdl.setMaterial(rock, water, gas);

% Activate the gravity
mdl.gravityOn = true;

% --- Boundary and initial conditions -------------------------------------

% Displacement boundary conditions
mdl.setDisplacementDirichletBCAtBorder('bottom',[0.0, 0.0]);
mdl.setDisplacementDirichletBCAtBorder('left',  [0.0, NaN]);
mdl.setDisplacementDirichletBCAtBorder('right', [0.0, NaN]);

% Liquid pressure
mdl.setPressureDirichletBCAtBorder('bottom',101025.0);
mdl.setInitialPressureAtDomain(101025.0);

% Gas pressure
mdl.setGasPressureDirichletBCAtBorder('top',101325.0);
mdl.setGasPressureDirichletBCAtBorder('bottom',101325.0);
mdl.setGasPressureDirichletBCAtBorder('left',101325.0);
mdl.setGasPressureDirichletBCAtBorder('right',101325.0);
mdl.setInitialGasPressureAtDomain(101325.0);

% --- Numerical model configuration ---------------------------------------

% Diagonalize compressibility matrix (mass lumping)
mdl.massLumping = true;
mdl.lumpStrategy = 2;

%% PRE-PROCESS

mdl.preComputations();

% Initialize the stress tensor
for el = 1:mdl.nelem
    for i = 1:mdl.element(el).type.nIntPoints
        mdl.element(el).type.intPoint(i).stressOld = [101325.0;101325.0;101325.0;0.0];
    end
end

%% PROCESS

% Transient analysis parameters
tinit = 0.1;          % Initial time
dt    = 0.1;          % Time step
tf    = 1.0;         % Final time
dtmax = 1.0;          % Time step
dtmin = 0.0001;          % Time step

anl = Anl_Transient("Picard");
anl.setUpTransientSolver(tinit,dt,tf,dtmax,dtmin,true);
anl.setRelativeConvergenceCriteria(true);
anl.maxIter = 15;
anl.run(mdl);

%% POST-PROCESS

% Plot pressure along a segment
Xi  = [0.0 , 0.0];
Xf  = [0.0 , Ly];
npts = 500;
mdl.plotPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotGasPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotCapillaryPressureAlongSegment(Xi, Xf, npts,'y')
mdl.plotField('CapillaryPressure');
mdl.plotField('GasPressure');
