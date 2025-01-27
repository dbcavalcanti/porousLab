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
        Young                = [];              % Young modulus (Pa)
        nu                   = [];              % Poisson ratio
        sy0                  = [];              % Initial yield stress (Pa)
        Kp                   = [];              % Plastic modulus (Pa)
        kappa                = [];              % Ratio between the uniaxial compressive strength and the uniaxial tensile strength
        DamageThreshold      = [];              % Damage threshold
        FractureEnergyMode1  = [];              % Fracture energy associated with mode 1 (N/m)
        rho                  = [];              % Density (kg/m3)
        K                    = 0.0;             % Intrinsic permeability (m2)  
        phi                  = 0.0;             % Porosity
        biot                 = 1.0;             % Biot's coefficient
        Ks                   = 1.0e25;          % Solid bulk modulus (Pa)
        Slr                  = 0.0;             % Residual liquid saturation
        Sgr                  = 0.0;             % Residual gas saturation 
        Pb                   = 0.0;             % Gas-entry pressure (Pa)
        lambda               = 0.0;             % Curve-fitting parameter
        liqRelPermeability   = 'BrooksCorey';   % Liquid relative permeability
        gasRelPermeability   = 'BrooksCorey';   % Gas relative permeability
        capillaryPressure    = 'BrooksCorey';   % Saturation degree function
        mechanical           = 'elastic';       % Mechanical constitutive law
        gravityOn            = false;           % Flag to consider gravity forces
        g                    = 9.806;           % Gravity accelaration (m/s2)
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
            if nargin == 1, this.id = id; end
            if nargin > 1
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

        function setMechanicalConstitutiveLaw(this,law)
            this.mechanical = law;
        end

        function setMechanicalProperties(this,E,nu)
            this.Young = E;
            this.nu = nu;
        end

        function setDensity(this,rho)
            this.rho = rho;
        end

        function rho = getDensity(this)
            rho = this.rho;
        end
    end
end