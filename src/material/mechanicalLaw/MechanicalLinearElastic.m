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

        %% Stress invariants
        % Considering a 2D stress tensor: Stress = [sxx, syy, szz, sxy]

        % I1 stress invariant
        function I1 = stressInvariantI1(stress)
            I1 = stress(1) + stress(2) + stress(3);
        end

        % I2 stress invariant
        function I2 = stressInvariantI3(stress)
            sxx = stress(1);
            syy = stress(2);
            szz = stress(3);
            sxy = stress(4);
            I2 =  sxx*syy + syy*szz + sxx*szz - sxy*sxy;
        end
        
        % J2 stress invariant
        function J2 = stressInvariantJ2(stress)
            I1 = stressInvariantI1(stress);
            I2 = stressInvariantI2(stress);
            J2 = I1 * I1 / 3.0 - I2;
            J2 = max(J2,0.0);
        end

        % Hydrostatic stress
        function sh = stressHydrostatic(stress)
            sh = stressInvariantI1(stress) / 3.0;
        end

        % von Mises stress
        function svm = stressVonMises(stress)
            J2 = stressInvariantJ2(stress);
            svm = sqrt(3.0 * J2);
        end

    end
end