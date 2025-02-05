%% Material_Elastic class
%
% This class defines an linear elastic stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
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