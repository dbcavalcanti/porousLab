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
classdef PorousMedia < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id                   = '';
        K                    = 0.0;             % Intrinsic permeability (m2)  
        phi                  = 0.0;             % Porosity
        biot                 = 0.0;             % Biot's coefficient
        Ks                   = 1.0e25;          % Solid bulk modulus (kPa)
        Slr                  = 0.0;             % Residual liquid saturation       
        Pb                   = 0.0;             % Gas-entry pressure (kPa)
        lambda               = 0.0;             % Curve-fitting parameter
        klrmin               = 1.0e-5;          % Minimun liquid relative permeability
        kgrmin               = 1.0e-5;          % Minimun gas relative permeability
        relativePermeability = 'BrooksCorey';
        capillaryPressure    = 'BrooksCorey';
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = PorousMedia(id, permeability, porosity, ...
                biotCoefficient, solidBulkModulus, ...
                residualLiquidSaturationDegree, ...
                gasEntryPressure,curveFittingParameter, ...
                relativePermeability,capillaryPressure)
            if nargin > 0
                this.id                   = id;
                this.K                    = permeability;
                this.phi                  = porosity;
                this.biot                 = biotCoefficient;
                this.Ks                   = solidBulkModulus;
                this.Slr                  = residualLiquidSaturationDegree;
                this.Pb                   = gasEntryPressure;
                this.lambda               = curveFittingParameter;
                this.relativePermeability = relativePermeability;
                this.capillaryPressure    = capillaryPressure;
            end
        end
    end
    %% Public methods
    methods
        function Se = effectiveSaturationDegree(this,Sl)
            Se = (Sl - this.Slr)/(1.0 - this.Slr);
        end

        function Km = intrinsicPermeabilityMatrix(this)
            Km = this.K * eye(2);
        end
    end
end