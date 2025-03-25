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
classdef Fluid < handle  
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id   = '';
        rho  = 1000.0;      % Density (kg/m3)
        mu   = 1.0e-3;      % Viscosity (Pa*s)
        K    = 1.0e25;      % Compressibility/Bulk modulus (1/Pa)
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Fluid(id, density, viscosity, compressibility)
            if nargin == 1, this.id = id; end
            if nargin > 0
                this.id   = id;
                this.rho  = density;
                this.mu   = viscosity;
                this.K    = compressibility;
            end
        end
    end
    %%
    methods
        function rho = getDensity(this,~)
            rho = this.rho;
        end
    end
end