%% MaterialDiscontinuity_H class
% This class represents a material discontinuity in a porous medium, 
% characterized by its initial aperture and the fluid properties. 
% It provides methods to compute the longitudinal permeability 
% coefficient and compressibility based on the material's properties.
%
%% Methods
% * *longitudinalPermeability*: Computes the longitudinal permeability 
%                               coefficient based on the cubic law.
% * *compressibility*: Computes  the compressibility of the material 
%                      discontinuity based on its aperture and fluid 
%                      properties.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef MaterialDiscontinuity_HM < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        initialAperture = 0.0;
        leakoff         = 1.0;
        fluid           = Fluid();
        parameters      = [];
        mechanical      = [];
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_HM(matData)
            % Create the material struct
            this.parameters = struct( ...
                'initialAperture',    matData.initialAperture, ...
                'normalStiffness',    matData.normalStiffness, ...
                'shearStiffness',     matData.shearStiffness,...
                'contactPenalization',matData.contactPenalization);
            this.initialAperture = matData.initialAperture;
            this.fluid = matData.fluid;
            this.leakoff = matData.leakoff;
            % Mechanical constitutive behavior
            if strcmp('elastic',matData.cohesiveLaw)
                this.mechanical = MechanicalCohesiveLinearElastic();
            end
        end
    end
    %% Public methods
    methods

        % -----------------------------------------------------------------
        % Evaluate the mechanical constitutive law
        function [stress,D] = mechanicalLaw(this,ip)
            [stress,D] = this.mechanical.eval(this.parameters,ip);
        end

        % -----------------------------------------------------------------
        % Get the number of state variables associated with the mechanical
        % constitutive law
        function nstVar = getNumberStateVar(this)
            nstVar = this.mechanical.nstVar;
        end

        % -----------------------------------------------------------------
        % Check if the material is elasto-plastic or not
        function flag = hasPlasticStrain(this)
            flag = this.mechanical.isElastoPlastic();
        end

        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficient based on the
        % cubic's law
        function kl = longitudinalPermeability(this)
            w = this.initialAperture();
            kl = w*w*w/12.0/this.fluid.mu;
        end

        %------------------------------------------------------------------
        % Computes  the compressibility of the material discontinuity 
        % based on its aperture and fluid properties.
        function c = compressibility(this)
            w = this.initialAperture();
            c = w / this.fluid.K;
        end

    end
end