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
classdef MechanicalLinearElastic < handle  
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        nstVar = 0;   % Number of state variables
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalLinearElastic()
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
        % Compute the elastic constitutive matrix
        function De = elasticConstitutiveMatrix(material,ip)

            % Elastic material properties
            E  = material.Young;
            nu = material.nu;

            if strcmp(ip.anm,'PlaneStress')

                De = [ 1.0    nu    0.0;
                       nu    1.0    0.0;
                       0.0   0.0  (1-nu)/2.0 ] * E/(1.0 - (nu*nu));

            elseif strcmp(ip.anm,'PlaneStrain')

                De = [ 1.0-nu    nu       0.0;
                         nu    1.0-nu     0.0;
                        0.0     0.0    (1-2.0*nu)/2.0 ] * E/(1.0 + nu)/(1.0 - 2.0*nu);

            else
                De = [];
            end
        end

    end
end