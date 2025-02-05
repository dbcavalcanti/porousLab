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
classdef MaterialDiscontinuity_M < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        parameters  = [];
        mechanical  = [];
    end  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialDiscontinuity_M(matData)
            % Create the material struct
            this.parameters = struct( ...
                'initialAperture',    matData.initialAperture, ...
                'normalStiffness',    matData.normalStiffness, ...
                'shearStiffness',     matData.shearStiffness,...
                'contactPenalization',matData.contactPenalization);
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

        function flag = hasPlasticStrain(this)
            flag = this.mechanical.isElastoPlastic();
        end

    end
end