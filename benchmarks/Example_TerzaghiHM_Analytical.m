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

%% = Solution from Kratos Poromechanics ===================================

% Pore-pressure profile for t = 43s
p43 = [0; 2158.9905; 4165.3042; 5897.2466; 7285.0317; 8315.9697; 9024.5928;
        9473.5127; 9732.501; 9861.5225; 9899.9883];

% Pore-pressure profile for t = 100s
p100 = [0; 1419.2424; 2792.1113; 4075.8826; 5234.4844; 6240.3516; 
        7074.8774; 7727.5127; 8193.8379; 8473.1396; 8566.0986 ];

% Pore-pressure profile for t = 262s
p262 = [0; 729.05255; 1440.1156; 2115.6511; 2739.01; 3294.8467; 3769.4941; 
        4151.2988; 4430.8975; 4601.4414; 4658.7573];

p3025 = [0
622.80933
1230.2755
1807.4351
2340.0732
2815.0754
3220.7495
3547.113
3786.1375
3931.9441
3980.948
];

% Pore-pressure profile for t = 500s
p500 = [0; 291.41986; 575.66394; 845.73334; 1094.9779; 1317.2605; 
        1507.1077; 1659.8448; 1771.7112; 1839.9521; 1862.8873 ];

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

% Plot the pore-pressure profile obtained with Kratos Poromechanics
dz    = 0.1;               % Depth increment (in meters)
z     = 0:dz:H;            % Depths within the soil layer
Z     = z/H;               % Normalized depth

% plot(p43/p0 ,1-z/H,'ok', 'LineWidth', 1.5, 'DisplayName', 'Kratos     : t =  43 s');
plot(p100,1-z/H,'xk', 'LineWidth', 1.5, 'DisplayName', 'Numerical     : t = 100 s');
% plot(p3025/p0,1-z/H,'sk', 'LineWidth', 1.5, 'DisplayName', 'Kratos    : t = 302.5 s');
% plot(p262/p0,1-z/H,'sk', 'LineWidth', 1.5, 'DisplayName', 'Kratos     : t = 262 s');
% plot(p500/p0,1-z/H,'dk', 'LineWidth', 1.5, 'DisplayName', 'Kratos     : t = 500 s');

% Figure formatting
xlabel('Pore Pressure (Pa)');
ylabel('y (m)');
set(gca,'FontSize',18)
grid on,box on; legend('Location','best');

%% = Error calculation =================================================

err = zeros(5,1);

Pnum = [p43 , p100 , p262, p3025, p500] / p0;
Panl = zeros(length(Z),length(t));

z     = H:-dz:0;            % Depths within the soil layer
Z     = z/H;               % Normalized depth

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

    % Save solution for the current time
    Panl(:,i) = p;

end

err = vecnorm(Pnum - Panl,1);

disp('Error of each time:')
disp(err)
