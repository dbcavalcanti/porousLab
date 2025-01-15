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
            elseif strcmp('Liakopoulos',matData.porousMedia.liqRelPermeability)
                this.liqRelativePermeability = RelativePermeabilityLiakopoulosLiquid();
            elseif strcmp('UMAT',matData.porousMedia.liqRelPermeability)
                curve = matData.porousMedia.klr_umat;
                this.liqRelativePermeability = RelativePermeabilityUMAT(curve(:,1),curve(:,2));
            end
            % Gas phase relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.gasRelPermeability)
                this.gasRelativePermeability = RelativePermeabilityBrooksCoreyGas();
            elseif strcmp('UMAT',matData.porousMedia.gasRelPermeability)
                curve = matData.porousMedia.kgr_umat;
                this.gasRelativePermeability = RelativePermeabilityUMAT(curve(:,1),curve(:,2));
            end
            % Saturation degree relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.capillaryPressure)
                this.capillaryPressure = CapillaryPressureBrooksCorey();
            elseif strcmp('UMAT',matData.porousMedia.capillaryPressure)
                curve = matData.porousMedia.SlPc_umat;
                this.capillaryPressure = CapillaryPressureUMAT(curve(:,1),curve(:,2));
            elseif strcmp('Liakopoulos',matData.porousMedia.capillaryPressure)
                this.capillaryPressure = CapillaryPressureLiakopoulos();
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
        % Compute the permeability matrices
        function [Kl,Kg] = permeabilityMtrcsPgPc(this,Sl,pg,pc)
            K = this.porousMedia.intrinsicPermeabilityMatrix();
            % Get fluids dynamic viscosity
            mul = this.liquidFluid.mu;
            mug = this.gasFluid.mu;
            % Get fluid densities
            rhol = this.liquidFluid.getDensity(pg-pc);
            rhog = this.gasFluid.getDensity(pg);
            % Compute relative permeability coefficients
            klr = this.liqRelativePermeability.calculate(Sl, this.porousMedia);
            kgr = this.gasRelativePermeability.calculate(Sl, this.porousMedia);
            % Permeability matrices
            Kl = K * klr / mul;
            Kg = K * kgr / mug * (rhog / rhol);
        end

        %------------------------------------------------------------------
        % Compute compressibility coefficients
        function [cll, clg, cgl, cgg] = compressibilityCoeffs(this,pc,Sl)
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
        function [ccc, ccg, cgc,cgg] = compressibilityCoeffsPgPc(this,Sl,pc,pg)
            % Get porous media parameters
            phi  = this.porousMedia.phi;
            % Get fluid densities
            rhol = this.liquidFluid.getDensity(pg-pc);
            rhog = this.gasFluid.getDensity(pg);
            % Derivative of the liquid saturation degree wrt pc
            dSldpc = this.capillaryPressure.derivativeSaturationDegree(pc, this.porousMedia);
            % Derivative of the gas density wrt to pg
            drhogdpg = this.derivativeGasDensityWrtGasPressure(rhog,pg);
            % Compressibility coefficients
            ccc =  phi * dSldpc;
            ccg = 0.0;
            cgc = -phi * dSldpc * (rhog / rhol);
            cgg =  phi * (1.0 - Sl) * drhogdpg / rhol;
        end
        %------------------------------------------------------------------
        % Compute the derivative of the gas density wrt to the gas pressure
        function drhogdpg = derivativeGasDensityWrtGasPressure(this,rhog,pg)
            % Pertubation value
            pert = sqrt(eps);
            % Compute the gas density given a pertubation at the gas pressure
            rhogPert = this.gasFluid.getDensity(pg + pert);
            % Compute the derivative
            drhogdpg = (rhogPert - rhog)/(pert);
        end

    end
end