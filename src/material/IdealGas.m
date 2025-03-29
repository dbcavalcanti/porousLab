%% Fluid class
%
% This class defines a fluid object
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef IdealGas < Fluid   
    %% Public attributes
    properties (SetAccess = private, GetAccess = public)
        R = 8.3144621;      % Universal gas constant (J/(mol*K)
        T = 293.15;         % Temperature (K)
        M = 0.02897;        % Molar mass (kg/mol)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = IdealGas(id, viscosity, compressibility)
            this = this@Fluid(id);
            this.rho = 0.0;
            if nargin > 1
                this.mu = viscosity;
                this.K  = compressibility;
            end
        end
    end
    %%
    methods
        % Get the fluid density based on the ideal Gas law
        function rho = getDensity(this,pg)
            rho = pg * this.M / (this.R * this.T);
        end
        % Set the value of the universal gas constant
        function setUniversalGasConstant(this,R)
            this.R = R;
        end
        % Set the value of the gas molar mass
        function setMolarMass(this,M)
            this.M = M;
        end
        % Set the value of the temperature
        function setTemperature(this,T)
            this.T = T;
        end
    end
end