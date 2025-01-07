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
        relativePermeability = [];   
        capillaryPressure    = [];
        fluids(2,1)          = Fluid();
        porousMedia          = PorousMedia();
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialTwoPhaseFlow(matData)
            this.fluids(1) = matData.fluids(1);
            this.fluids(2) = matData.fluids(2);
            this.porousMedia = matData.porousMedia;
            if strcmp('BrooksCorey',matData.porousMedia.relativePermeability)
                this.relativePermeability = RelativePermeabilityBrooksCorey();
                this.capillaryPressure    = CapillaryPressureBrooksCorey();
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
            mul = this.fluids(1).mu;
            mug = this.fluids(2).mu;
            % Compute relative permeability coefficients
            klr = this.relativePermeability.liquidRelativePermeability(Sl, this.porousMedia);
            kgr = this.relativePermeability.gasRelativePermeability(Sl, this.porousMedia);
            % Permeability matrices
            Kl = K * klr / mul;
            Kg = K * kgr / mug;
        end

        %------------------------------------------------------------------
        % Compute the permeability matrices
        function [Kl,Kg] = permeabilityMtrcsPgPc(this,Sl)
            K = this.porousMedia.intrinsicPermeabilityMatrix();
            % Get fluids dynamic viscosity
            mul = this.fluids(1).mu;
            mug = this.fluids(2).mu;
            % Get fluids density
            rhoL = this.fluids(1).rho/1000.0;
            rhoG = this.fluids(2).rho/1000.0;
            % Compute relative permeability coefficients
            klr = this.relativePermeability.liquidRelativePermeability(Sl, this.porousMedia);
            kgr = this.relativePermeability.gasRelativePermeability(Sl, this.porousMedia);
            % Permeability matrices
            Kl = K * klr * rhoL / mul;
            Kg = K * kgr * rhoG / mug;
        end

        %------------------------------------------------------------------
        % Compute compressibility coefficients
        function [cll,cgg,clg,cgl] = compressibilityCoeffs(this,pc,Sl)
            % Get porous media parameters
            biot = this.porousMedia.biot;
            phi  = this.porousMedia.phi;
            Ks   = this.porousMedia.Ks;
            % Get the fluids bulk modulus
            Klb  = this.fluids(1).K;
            Kgb  = this.fluids(2).K;
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
            % Get fluids density
            rhoL = this.fluids(1).rho/1000.0;
            rhoG = this.fluids(2).rho/1000.0;
            % Derivative of the liquid saturation degree wrt pc
            dSldpc = this.capillaryPressure.derivativeSaturationDegree(pc, this.porousMedia);
            % Compressibility coefficients
            cgc = -phi * rhoG * dSldpc;
            ccc =  phi * rhoL * dSldpc;
        end
    end
end