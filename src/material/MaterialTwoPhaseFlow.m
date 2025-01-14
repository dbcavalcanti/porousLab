%% Material class
%
% This class defines an abstract stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef MaterialTwoPhaseFlow < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        liqRelativePermeability = [];
        gasRelativePermeability = [];
        capillaryPressure       = [];
        liquidFluid             = Fluid();
        gasFluid                = Fluid();
        porousMedia             = PorousMedia();
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialTwoPhaseFlow(matData)
            this.liquidFluid = matData.liquidFluid;
            this.gasFluid    = matData.gasFluid;
            this.porousMedia = matData.porousMedia;
            % Liquid phase relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.liqRelPermeability)
                this.liqRelativePermeability = RelativePermeabilityBrooksCoreyLiquid();
            end
            % Gas phase relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.gasRelPermeability)
                this.gasRelativePermeability = RelativePermeabilityBrooksCoreyGas();
            end
            % Saturation degree relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.capillaryPressure)
                this.capillaryPressure = CapillaryPressureBrooksCorey();
            elseif strcmp('UMAT',matData.porousMedia.capillaryPressure)
                curve = matData.porousMedia.SlPc_umat;
                this.capillaryPressure = CapillaryPressureUMAT(curve(:,1),curve(:,2));
            end
        end
    end
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Get the liquid saturation degree
        function Sl = saturationDegree(this,pc)
            Sl = this.capillaryPressure.saturationDegree(pc, this.porousMedia);
        end

        %------------------------------------------------------------------
        % Compute the permeability matrices
        function [Kl,Kg] = permeabilityMtrcs(this,Sl)
            K = this.porousMedia.intrinsicPermeabilityMatrix();
            % Get fluids dynamic viscosity
            mul = this.liquidFluid.mu;
            mug = this.gasFluid.mu;
            % Compute relative permeability coefficients
            klr = this.liqRelativePermeability.calculate(Sl, this.porousMedia);
            kgr = this.gasRelativePermeability.calculate(Sl, this.porousMedia);
            % Permeability matrices
            Kl = K * klr / mul;
            Kg = K * kgr / mug;
        end

        %------------------------------------------------------------------
        % Compute compressibility coefficients
        function [cll,cgg,clg,cgl] = compressibilityCoeffs(this,pc,Sl)
            % Get porous media parameters
            biot = this.porousMedia.biot;
            phi  = this.porousMedia.phi;
            Ks   = this.porousMedia.Ks;
            % Get the fluids bulk modulus
            Klb  = this.liquidFluid.K;
            Kgb  = this.gasFluid.K;
            % Gas saturation degree
            Sg   = 1.0 - Sl;
            % Derivative of the liquid saturation degree wrt pc
            dSldpc = this.capillaryPressure.derivativeSaturationDegree(pc, this.porousMedia);
            % Compressibility coefficients
            cll = Sl * ((biot - phi) / Ks) * (Sl + dSldpc * pc) - phi * dSldpc + phi * Sl / Klb;
            clg = Sl * ((biot - phi) / Ks) * (Sg - dSldpc * pc) + phi * dSldpc;
            cgl = Sg * ((biot - phi) / Ks) * (Sl + dSldpc * pc) + phi * dSldpc;
            cgg = Sg * ((biot - phi) / Ks) * (Sg - dSldpc * pc) - phi * dSldpc + phi * Sg / Kgb;
        end

        %------------------------------------------------------------------
        % Compute compressibility coefficients
        function [cgc,ccc] = compressibilityCoeffsPgPc(this,pc)
            % Get porous media parameters
            phi  = this.porousMedia.phi;
            % Derivative of the liquid saturation degree wrt pc
            dSldpc = this.capillaryPressure.derivativeSaturationDegree(pc, this.porousMedia);
            % Compressibility coefficients
            cgc = -phi * dSldpc;
            ccc =  phi * dSldpc;
        end
    end
end