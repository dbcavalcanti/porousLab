%% Material_M class
%
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef MaterialDiscontinuity_H < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        initialAperture = 0.0;
        fluid           = Fluid();
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_H(matData)
            this.initialAperture = matData.initialAperture;
            this.fluid = matData.fluid;
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Computes the longitudinal permeability coefficient based on the
        % cubic's law
        function kl = longitudinalPermeability(this)
            w = this.initialAperture();
            kl = w*w*w/12.0/this.fluid.mu;
        end

        %------------------------------------------------------------------
        function c = compressibility(this)
            w = this.initialAperture();
            c = w / this.fluid.K;
        end


    end
end