%% Material_H2 class
% This class defines the material properties and behavior for a porous 
% medium involving liquid and gas phases. It includes methods for 
% computing relative permeabilities, saturation degrees, permeability 
% matrices, and compressibility coefficients for the liquid and gas phases.
%
%% Methods
% * *saturationDegree*: Computes the liquid saturation degree based on 
%                       the capillary pressure pc.
% * *permeabilityMtrcs*: Computes the permeability matrices for the liquid
%                        and gas phases based on the liquid saturation 
%                        degree Sl, liquid pressure pl, and gas 
%                        pressure pg.
% * *permeabilityMtrcsPgPc*: Computes the permeability matrices for the
%                            liquid and gas phases considering the 
%                            gas-to-liquid density ratio.
% * *compressibilityCoeffs*: Computes the compressibility coefficients for
%                            the liquid and gas phases based on the liquid 
%                            saturation degree Sl, liquid pressure pl, 
%                            and gas pressure pg.
% * *compressibilityCoeffsPgPc*: Computes the compressibility coefficients
%                                considering the gas-to-liquid 
%                                density ratio.
% * *derivativeGasDensityWrtGasPressure*: Computes the derivative of the 
%                                         gas density with respect to the 
%                                         gas pressure Pg.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef Material_H2 < handle    
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
        function this = Material_H2(matData)
            this.liquidFluid = matData.liquidFluid;
            this.gasFluid    = matData.gasFluid;
            this.porousMedia = matData.porousMedia;
            % Liquid phase relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.liqRelPermeability)
                this.liqRelativePermeability = RelativePermeabilityBrooksCoreyLiquid();
            elseif strcmp('Liakopoulos',matData.porousMedia.liqRelPermeability)
                this.liqRelativePermeability = RelativePermeabilityLiakopoulosLiquid();
            elseif strcmp('PolynomialLiquid',matData.porousMedia.liqRelPermeability)
                this.liqRelativePermeability = RelativePermeabilityPolynomialLiquid();
            elseif strcmp('UMAT',matData.porousMedia.liqRelPermeability)
                curve = matData.porousMedia.klr_umat;
                this.liqRelativePermeability = RelativePermeabilityUMAT(curve(:,1),curve(:,2));
            end
            % Gas phase relative permeability function
            if strcmp('BrooksCorey',matData.porousMedia.gasRelPermeability)
                this.gasRelativePermeability = RelativePermeabilityBrooksCoreyGas();
            elseif strcmp('PolynomialGas',matData.porousMedia.gasRelPermeability)
                this.gasRelativePermeability = RelativePermeabilityPolynomialGas();
            elseif strcmp('UMAT',matData.porousMedia.gasRelPermeability)
                curve = matData.porousMedia.kgr_umat;
                this.gasRelativePermeability = RelativePermeabilityUMAT(curve(:,1),curve(:,2));
            end
            % Saturation degree function
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
        % Get the liquid saturation degree derivative wrt the capillary
        % pressure
        function dSldpc = derivativeSaturationDegree(this,pc)
            dSldpc = this.capillaryPressure.derivativeSaturationDegree(pc, this.porousMedia);
        end

        %------------------------------------------------------------------
        % Compute the relative permeabilities
        function [klr, kgr] = relativePermeabilities(this,Sl)
            klr = this.liqRelativePermeability.calculate(Sl, this.porousMedia);
            kgr = this.gasRelativePermeability.calculate(Sl, this.porousMedia);
        end

        %------------------------------------------------------------------
        % Compute the relative permeabilities
        function [dklrdSl, dkgrdSl] = derivativeRelPerm(this,Sl)
            dklrdSl = this.liqRelativePermeability.derivative(Sl, this.porousMedia);
            dkgrdSl = this.gasRelativePermeability.derivative(Sl, this.porousMedia);
        end

        %------------------------------------------------------------------
        % Compute the permeability matrices
        function [Kll, Klg, Kgl, Kgg] = permeabilityMtrcs(this,Sl,pl,pg)
            K = this.porousMedia.intrinsicPermeabilityMatrix();
            % Get fluids dynamic viscosity
            mul = this.liquidFluid.mu;
            mug = this.gasFluid.mu;
            % Compute relative permeability coefficients
            klr = this.liqRelativePermeability.calculate(Sl, this.porousMedia);
            kgr = this.gasRelativePermeability.calculate(Sl, this.porousMedia);
            % Permeability matrices
            Kll = K * klr / mul;
            Klg = zeros(2);
            Kgl = zeros(2);
            Kgg = K * kgr / mug;
        end

        %------------------------------------------------------------------
        % Compute the permeability matrices
        function [Kll, Klg, Kgl, Kgg] = permeabilityMtrcsPgPc(this,Sl,pl,pg)
            K = this.porousMedia.intrinsicPermeabilityMatrix();
            % Get fluids dynamic viscosity
            mul = this.liquidFluid.mu;
            mug = this.gasFluid.mu;
            % Get fluid densities
            rhol = this.liquidFluid.getDensity(pl);
            rhog = this.gasFluid.getDensity(pg);
            % Compute relative permeability coefficients
            klr = this.liqRelativePermeability.calculate(Sl, this.porousMedia);
            kgr = this.gasRelativePermeability.calculate(Sl, this.porousMedia);
            % Permeability matrices
            Kll = - K * klr / mul;
            Klg = K * klr / mul;
            Kgl = zeros(2);
            Kgg = K * kgr / mug * (rhog / rhol);
        end

        %------------------------------------------------------------------
        % Compute the compressibility coefficients
        function [cll, clg, cgl, cgg] = compressibilityCoeffs(this,Sl,pl,pg)
            % Get porous media parameters
            biot = this.porousMedia.biot;
            phi  = this.porousMedia.phi;
            Ks   = this.porousMedia.Ks;
            % Get the fluids bulk modulus
            Klb  = this.liquidFluid.K;
            Kgb  = this.gasFluid.K;
            % Gas saturation degree
            Sg   = 1.0 - Sl;
            % Capillary pressure
            pc = pg - pl;
            % Derivative of the liquid saturation degree wrt pc
            dSldpc = this.capillaryPressure.derivativeSaturationDegree(pc, this.porousMedia);
            % Compressibility coefficients
            cll = Sl * ((biot - phi) / Ks) * (Sl + dSldpc * pc) - phi * dSldpc + phi * Sl / Klb;
            clg = Sl * ((biot - phi) / Ks) * (Sg - dSldpc * pc) + phi * dSldpc;
            cgl = Sg * ((biot - phi) / Ks) * (Sl + dSldpc * pc) + phi * dSldpc;
            cgg = Sg * ((biot - phi) / Ks) * (Sg - dSldpc * pc) - phi * dSldpc + phi * Sg / Kgb;
        end

        %------------------------------------------------------------------
        % Compute the compressibility coefficients
        function [ccc, ccg, cgc,cgg] = compressibilityCoeffsPgPc(this,Sl,pl,pg)
            % Get porous media parameters
            phi  = this.porousMedia.phi;
            % Get fluid densities
            rhol = this.liquidFluid.getDensity(pl);
            rhog = this.gasFluid.getDensity(pg);
            % Capillary pressure
            pc = pg - pl;
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