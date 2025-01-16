%% RegularElementPcPg class
%
% This class defines a finite element of a two-phase flow formulation using
% the capillary pressure (Pc) and the gas pressure (Pg) as primary
% variables.
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef RegularElementPcPg < RegularElement    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElementPcPg(type, node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy)
            this = this@RegularElement(type, node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy);
        end
    end
    
    %% Public methods
    methods

        % -----------------------------------------------------------------
        % Compute the permeability tensors
        function [kll, klg, kgl, kgg] = permeabilityTensors(~,ip,pg,pc,Sl)
             [kll, klg, kgl, kgg] = ip.constitutiveMdl.permeabilityMtrcsPgPc(Sl,pg-pc,pg);
        end

        % -----------------------------------------------------------------
        % Compute the compressibility coefficients
        function [cll, clg, cgl, cgg] = compressibilityCoeffs(~,ip,pg,pc,Sl)
             [cll, clg, cgl, cgg] =  ip.constitutiveMdl.compressibilityCoeffsPgPc(Sl,pg-pc,pg);
        end

        %------------------------------------------------------------------
        % Add contribution of the gravity forces to the external force vct
        function [fec,feg] = addGravityForces(this,fec,feg,Bp,kl,kg,pl,pg,c)

            % Get gravity vector
            grav = this.mat.porousMedia.g * this.mat.porousMedia.b;

            % Get fluid densities
            rhol = this.mat.liquidFluid.getDensity(pl);
            rhog = this.mat.gasFluid.getDensity(pg);

            kl = -kl;

            % Compute the contribution of the gravitational forces
            fec = fec + Bp' * kl * rhol * grav * c;
            feg = feg + Bp' * kg * rhog * grav * c * (rhog/rhol);

        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the capillary pressure
        function pc = getNodalCapillaryPressure(this)
            pc = this.ue(1:this.nglp);
        end

        %------------------------------------------------------------------
        % Function to get the nodal values of the liquid pressure
        function pl = getNodalLiquidPressure(this)
            pc = this.getNodalCapillaryPressure();
            pg = this.getNodalGasPressure();
            pl = pg - pc;
        end

    end
end