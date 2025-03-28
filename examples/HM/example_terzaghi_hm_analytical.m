%% ****************** TERZAGHI 1D CONSOLIDATION PROBLEM *******************

% Author: Danilo Cavalcanti
% July 2023

% References:
% Nguyen, Vinh Phu, Haojie Lian, Timon Rabczuk, and Stéphane Bordas.
% "Modelling Hydraulic Fractures in Porous Media Using Flow Cohesive
% Interface Elements.” Engineering Geology, Special Issue: Characterisation
% of Fractures in Rock: from Theory to Practice (ROCKFRAC), 
% 225 (July 20, 2017): 68–82. https://doi.org/10.1016/j.enggeo.2017.04.010.
%
% Segura, J. M., and I. Carol. 
% “Coupled HM Analysis Using Zero-Thickness Interface Elements with 
% Double Nodes—Part II: Verification and Application.” International 
% Journal for Numerical and Analytical Methods in Geomechanics 32, 
% no. 18 (2008): 2103–23. https://doi.org/10.1002/nag.730.
%% ========================================================================

clear, clc,% close all

%% = Physical parameters ==================================================

% Solid parameters
E     = 1.0e6;              % Young's modulus (in Pa)
nu    = 0.3;                % Poisson's ratio (-)
n     = 0.3;                % Porosity (-)
k     = 1.15741e-12;        % Intrinsic permeability (in m^2)
Ks    = 1.0e25;             % Solid bulk modulus (in Pa)

% Fluid parameters
muf   = 1e-3;               % Fluid dynamic viscosity (Pa*s)

% Porous medium parameters
alpha = 1.0;                % Biot's coefficient
Sm    = 0.0;                % Storage coefficient
K     = E/3/(1-2*nu);       % Bulk modulus (Pa)
G     = E/2/(1+nu);         % Shear modulus (Pa)
mv    = 1 / ((4/3)*G + K);  % Confined compressibility (1/Pa)

% Consolidation coefficient
Cv    = k / (muf * (Sm + alpha*alpha * mv));    

%% = Problem parameters ===================================================

% Initial excess pore pressure (Pa)
p0    = 1.0e4;             

% Total thickness of the soil layer (m)
H     = 1;            

% Consolidation time (in seconds)        
% t     = [43, 100, 262, 302.5, 500];
t     = 100.0;

% Number of terms in the expansion
np    = 1000;

% Define depth intervals
dz    = 0.001;             % Depth increment (in meters)
z     = 0:dz:H;            % Depths within the soil layer
Z     = z/H;               % Normalized depth

%% = Terzaghi solution =================================================

% Initialize the figure for plotting
figure;
hold on

% Calculate the pore pressure profile
for i = 1:length(t)

    % Initialize the pressure vector
    p = zeros(size(z));

    % Calculate the time factor
    T = Cv * t(i) / (H^2);

    % Terzaghi's solution
    for j = 1:np
        M = pi/2*(2*j-1);
        p = p + 4/pi * ((-1)^(j-1)/(2*j-1)) * cos(M*Z) * exp(-M*M*T);
    end

    % Plot the pressure vector
    plot(p*p0, z/H, 'LineWidth', 2, ...
        'DisplayName', ['Analytical: t = ',num2str(t(i)), ' s']);

    % Print pressure at the bottom
    disp(['Pressure at the bottom edge: ', num2str(p0*p(1)) , ' Pa'])

end

% Figure formatting
xlabel('Pore Pressure (Pa)');
ylabel('y (m)');
set(gca,'FontSize',18)
grid on,box on; legend('Location','best');
