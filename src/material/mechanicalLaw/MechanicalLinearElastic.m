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

        %% Elastic constants
 
        function G = shearModulus(~,material)
            G = material.Young / (2.0 * (1.0 + material.nu));
        end

        function K = bulkModulus(~,material)
            K = material.Young / (3.0 * (1.0 - 2.0*material.nu));
        end
        
        %% Stress utilities
        % Considering a 2D stress tensor: Stress = [sxx, syy, szz, sxy]

        % Hydrostatic stress
        function sh = hydrostaticStress(this,stress)
            sh = this.stressInvariantI1(stress) / 3.0;
        end

        % Deviatoric stress
        function sd = deviatoricStress(this,stress)
            sh = this.hydrostaticStress(stress);
            Im = [1.0;1.0;1.0;0.0];
            sd = stress - sh * Im;
        end

        % I1 stress invariant
        function I1 = stressInvariantI1(~,stress)
            I1 = stress(1) + stress(2) + stress(3);
        end

        % I2 stress invariant
        function I2 = stressInvariantI2(~,stress)
            sxx = stress(1);
            syy = stress(2);
            szz = stress(3);
            sxy = stress(4);
            I2 =  sxx*syy + syy*szz + sxx*szz - sxy*sxy;
        end

        % J2 stress invariant
        function J2 = stressInvariantJ2(this,stress)
            I1 = this.stressInvariantI1(stress);
            I2 = this.stressInvariantI2(stress);
            J2  = I1*I1/3.0 - I2;
            J2  = max(J2,0.0);
        end

        % Gradient of the J2 stress invariant
        function dJ2 = gradientJ2(~,stress)
            sxx = stress(1);
            syy = stress(2);
            szz = stress(3);
            sxy = stress(4);
            dJ2 = zeros(4,1);
            dJ2(1) = (2.0 * sxx - syy - szz)/3.0;
            dJ2(2) = (2.0 * syy - sxx - szz)/3.0;
            dJ2(3) = (2.0 * szz - syy - sxx)/3.0;
            dJ2(4) = 2.0 * sxy;
        end

        % Hessian of the J2 stress invariant
        function d2J2 = hessianJ2(~)
            d2J2 = [  2.0 , -1.0 , -1.0 , 0.0;
                     -1.0 ,  2.0 , -1.0 , 0.0;
                     -1.0 , -1.0 ,  2.0 , 0.0;
                      0.0 ,  0.0  , 0.0 , 6.0 ]/3.0;
        end

        % von Mises stress
        function sVM = vonMisesStress(this,stress)
            J2 = this.stressInvariantJ2(stress);
            sVM = sqrt(3.0 * J2);
        end

        % von Mises stress gradient
        function dsVM = vonMisesStressGradient(this,stress)
            J2      = this.stressInvariantJ2(stress);
            dsVMdJ2 = 0.5*sqrt(3.0/J2);
            dJ2     = this.gradientJ2(stress);
            dsVM    = dsVMdJ2 * dJ2;
        end

        % von Mises stress hessian matrix
        function d2sVM = vonMisesStressHessian(this,stress)
            J2    = this.stressInvariantJ2(stress);
            dJ2   = this.gradientJ2(stress);
            d2J2  = this.hessianJ2();
            dsVMdJ2 = 0.5*sqrt(3.0/J2);
            d2sVM = dsVMdJ2 * d2J2 - (0.25 * sqrt(3.0/J2/J2/J2))*(dJ2 * dJ2');
        end

        % Get plane stress constitutive matrix
        function D = planeStressConstitutiveMatrix(~,D)
            % Get the sub-matrices
            D11 = [ D(1,1) , D(1,2) , D(1,4) ;
                    D(2,1) , D(2,2) , D(2,4) ;
                    D(4,1) , D(4,2) , D(4,4) ];
            D21 = [ D(3,1) ; D(3,2) ; D(3,4)];
            D12 = [ D(1,3) ; D(2,3) ; D(4,3)];
            D22 = D(3,3);
            % Compute the plane stress matrix
            Dt  = D11 - D12 * D21' / D22;
            % Assemble
            D = [ Dt(1,1) , Dt(1,2) , 0.0 , Dt(1,3) ;
                  Dt(2,1) , Dt(2,2) , 0.0 , Dt(2,3) ;
                   0.0    ,   0.0   , 1.0 ,   0.0   ;
                  Dt(3,1) , Dt(3,2) , 0.0 , Dt(3,3) ];
        end

        %% Strain utilities
        % Considering a 2D strain tensor: strain = [exx, eyy, ezz, 2exy]

        % Compute the elastic out-of-plane strain
        function elasticOutOfPlaneStrain(~,material,ip)
            if strcmp(ip.anm,'PlaneStress')
                nu = material.nu;
                ip.strain(3) = -nu/(1.0-nu)*(ip.strain(1)+ip.strain(2));
            end
        end

        % I1 strain invariant
        function I1 = strainInvariantI1(~,strain)
            I1 = strain(1) + strain(2) + strain(3);
        end

        % I2 strain invariant
        function I2 = strainInvariantI2(~,strain)
            exx = strain(1);
            eyy = strain(2);
            ezz = strain(3);
            exy = strain(4) / 2.0;
            I2 =  exx*eyy + eyy*ezz + exx*ezz - exy*exy;
        end

        % J2 strain invariant
        function J2 = strainInvariantJ2(this,strain)
            I1 = this.strainInvariantI1(strain);
            I2 = this.strainInvariantI2(strain);
            J2  = I1*I1/3.0 - I2;
            J2  = max(J2,0.0);
        end

        % Gradient of the J2 stress invariant
        function dI1 = gradientStrainInvariantI1(~)
            dI1 = [1.0;1.0;1.0;0.0];
        end

        % Gradient of the J2 strain invariant
        function dJ2 = gradientStrainInvariantJ2(~,strain)
            exx = strain(1);
            eyy = strain(2);
            ezz = strain(3);
            exy = strain(4) / 2.0;
            dJ2 = zeros(4,1);
            dJ2(1) = (2.0 * exx - eyy - ezz)/3.0;
            dJ2(2) = (2.0 * eyy - exx - ezz)/3.0;
            dJ2(3) = (2.0 * ezz - eyy - exx)/3.0;
            dJ2(4) = 2.0 * exy;
        end

        function I4 = fourthOrderSymTensor(~)
            I4 = [ 1.0 , 0.0 , 0.0 , 0.0;
                   0.0 , 1.0 , 0.0 , 0.0;
                   0.0 , 0.0 , 1.0 , 0.0;
                   0.0 , 0.0 , 0.0 , 0.5 ];
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