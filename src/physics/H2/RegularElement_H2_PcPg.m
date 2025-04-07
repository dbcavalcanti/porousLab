%% RegularElement_H2_PcPg Class
% This class represents a finite element for a two-phase flow formulation 
% using capillary pressure (Pc) and gas pressure (Pg) as primary variables. 
% It extends the _RegularElement_H2_ class and provides methods for 
% computing permeability tensors, compressibility coefficients, and 
% gravitational force contributions, as well as retrieving nodal values 
% for capillary and liquid pressures.
% 
%% Methods
% * *permeabilityTensors*: Computes the permeability tensors for the 
%                          element.
% * *compressibilityCoeffs*: Computes the compressibility coefficients 
%                            for the element.
% * *addGravityForces*: Adds the contribution of gravity forces to the 
%                       external force vector.
% * *getNodalCapillaryPressure*: Retrieves the nodal capillary pressure 
%                                values.
% * *getNodalLiquidPressure*: Retrieves the nodal liquid pressure values.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef RegularElement_H2_PcPg < RegularElement_H2  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = RegularElement_H2_PcPg(type, node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy, ...
                isAxisSymmetric)
            this = this@RegularElement_H2(type, node, elem, t, ...
                mat, intOrder, glp, glpg, massLumping, lumpStrategy, ...
                isAxisSymmetric);
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
            grav = this.g * this.mat.porousMedia.b;

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