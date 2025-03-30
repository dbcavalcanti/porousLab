%% MechanicalLinearElastic class
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
classdef MechanicalNonlinearElasticDruckerPrager < MechanicalLinearElastic  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalNonlinearElasticDruckerPrager()
            this = this@MechanicalLinearElastic();
        end
    end

    %% Public methods
    methods
        % Compute the stress vector and the constitutive matrix
        function [stress,De] = eval(this,material,ip)

            % Compute out-of-plane strain (plane stress problems)
            this.elasticOutOfPlaneStrain(material,ip)

            % Decompose the strain tensor
            ev = this.volumetricStrain(ip.strain);
            ed = this.deviatoricStrain(ip.strain);

            % Norm of the strain tensor
            normE = this.normTensor(ip.strain);

            % Get material parameters
            K = this.bulkModulus(material);

            % Compute the stress vector
            stress(1:3) = K * ev + 2 * mu * ed(1:3);

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);

        end
        %% Elastic constants
        function G = shearModulus(~,material)
            G = material.Young / (2.0 * (1.0 + material.nu));
        end

        function K = bulkModulus(~,material)
            K = material.Young / (3.0 * (1.0 - 2.0*material.nu));
        end 
    end

    %% Public methods
    methods (Static)
        function flag = isElastoPlastic()
            flag = false;
        end
        %% Elastic tensors
        % Compute the elastic constitutive matrix
        function De = elasticConstitutiveMatrix(material,ip)

            % Elastic material properties
            E  = material.Young;
            nu = material.nu;

            if strcmp(ip.anm,'PlaneStress')

                c = E/(1.0 - (nu*nu));
                De = [  c   ,   c*nu , 0.0  ,  0.0;
                       c*nu ,    c   , 0.0  ,  0.0;
                       0.0  ,   0.0  , 1.0  ,  0.0;
                       0.0  ,   0.0  , 0.0  , c*(1-nu)/2.0 ];
                

            elseif strcmp(ip.anm,'PlaneStrain')

                De = [ 1.0-nu ,   nu   ,   nu   ,    0.0;
                         nu   , 1.0-nu ,   nu   ,    0.0;
                         nu   ,  nu    , 1.0-nu ,    0.0;
                        0.0   ,  0.0   ,   0.0  , (1-2.0*nu)/2.0 ];

                De = De * E/(1.0 + nu)/(1.0 - 2.0*nu);

            else
                De = [];
            end
        end

        %------------------------------------------------------------------
        % Compute the elastic flexibility matrix
        function Ce = elasticFlexibilityMatrix(material,ip)

            % Elastic material properties
            E  = material.Young;
            nu = material.nu;

            if strcmp(ip.anm,'PlaneStress')

                Ce = [  1.0/E ,  -nu/E  ,  0.0  ,  0.0;
                        -nu/E ,  1.0/E  ,  0.0  ,  0.0;
                         0.0  ,  0.0    ,  1.0  ,  0.0;
                         0.0  ,  0.0    ,  0.0  , 2.0*(1+nu)/E ];

            elseif strcmp(ip.anm,'PlaneStrain')

                Ce = [  1.0/E ,  -nu/E  ,  -nu/E ,  0.0;
                        -nu/E ,  1.0/E  ,  -nu/E ,  0.0;
                        -nu/E ,  -nu/E  ,  1.0/E ,  0.0;
                         0.0  ,  0.0    ,  0.0   , 2.0*(1+nu)/E ];
            else
                Ce = [];
            end
        end   
    end
end