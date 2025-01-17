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
        Sgr                  = 0.0;             % Residual gas saturation 
        Pb                   = 0.0;             % Gas-entry pressure (kPa)
        lambda               = 0.0;             % Curve-fitting parameter
        liqRelPermeability   = 'BrooksCorey';   % Liquid relative permeability
        gasRelPermeability   = 'BrooksCorey';   % Gas relative permeability
        capillaryPressure    = 'BrooksCorey';   % Saturation degree function
        gravityOn            = false;           % Flag to consider gravity forces
        g                    = 9.81;            % Gravity accelaration (m/s2)
        b                    = [0.0;-1.0];      % Gravity force direction vector       
        SlPc_umat           = [];
        klr_umat            = [];
        kgr_umat            = [];
        m                   = 1;                % Expoent for the polynomial relationships
    end
    properties (SetAccess = protected, GetAccess = public)
        klrmin               = 1.0e-9;          % Minimum liquid relative permeability
        kgrmin               = 1.0e-9;          % Minimum gas relative permeability
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = PorousMedia(id, permeability, porosity, ...
                biotCoefficient, solidBulkModulus, ...
                residualLiquidSaturationDegree, ...
                residualGasSaturationDegree,...
                gasEntryPressure,curveFittingParameter, ...
                liqRelPermeability, gasRelPermeability, ...
                capillaryPressure)
            if nargin > 0
                this.id                   = id;
                this.K                    = permeability;
                this.phi                  = porosity;
                this.biot                 = biotCoefficient;
                this.Ks                   = solidBulkModulus;
                this.Slr                  = residualLiquidSaturationDegree;
                this.Sgr                  = residualGasSaturationDegree;
                this.Pb                   = gasEntryPressure;
                this.lambda               = curveFittingParameter;
                this.liqRelPermeability   = liqRelPermeability;
                this.gasRelPermeability   = gasRelPermeability;
                this.capillaryPressure    = capillaryPressure;
            end
        end
    end
    %% Public methods
    methods
        function Se = effectiveSaturationDegree(this,Sl)
            Se = (Sl - this.Slr)/(1.0 - this.Slr - this.Sgr);
        end

        function Km = intrinsicPermeabilityMatrix(this)
            Km = this.K * eye(2);
        end

        function setMinLiquidRelPermeability(this,klrmin)
            this.klrmin = klrmin;
        end

        function setMinGasRelPermeability(this,kgrmin)
            this.kgrmin = kgrmin;
        end

        function setUMATCapillaryPressureCurve(this,curve)
            curve = sortrows(curve, 2);
            this.SlPc_umat = curve;
        end

        function setUMATLiquidRelPermCurve(this,curve)
            curve = sortrows(curve, 1);
            this.klr_umat = curve;
        end

        function setUMATGasRelPermCurve(this,curve)
            curve = sortrows(curve, 1);
            this.kgr_umat = curve;
        end
    end
end