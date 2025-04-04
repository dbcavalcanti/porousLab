%% MechanicalCohesiveLinearElastic Class
% This class implements a linear elastic cohesive law for mechanical 
% materials. It provides methods to compute the stress vector and the 
% constitutive matrix based on the material properties and the strain 
% state at integration points.
%
%% Methods
% * *eval*: Computes the stress vector and the constitutive matrix for the 
%           given material and integration point.
% * *isElastoPlastic*: Static method that indicates that the material is 
%                      not elasto-plastic.
% * *elasticConstitutiveMatrix*: Static method that computes the elastic 
%                                constitutive matrix based on the material 
%                                properties and the strain state at the 
%                                integration point.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalCohesiveLinearElastic < handle  
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        nstVar = 0;   % Number of state variables
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalCohesiveLinearElastic()
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,De] = eval(this,material,ip)

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);

            % Stress vector
            stress = De * (ip.strain - ip.strainOld) + ip.stressOld;

        end    
    end
    %% Public methods
    methods (Static)
        %------------------------------------------------------------------
        % Flag to return that the material is not elasto-plastic
        function flag = isElastoPlastic()
            flag = false;
        end
        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix
        function De = elasticConstitutiveMatrix(material,ip)

            % Elastic material properties
            kn = material.normalStiffness;
            ks = material.shearStiffness;

            % Normal jump
            dn = ip.strain(2);

            % Closure model
            if dn < -1.0e-6
                kn = kn * 1.0e3;
            end

            % Assemble constitutive matrix
            De = [ ks   0.0;
                   0.0   kn ];
        end 
    end
end