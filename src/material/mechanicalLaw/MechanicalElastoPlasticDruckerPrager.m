%% MechanicalElastoPlasticDruckerPrager Class
% This class implements the Drucker-Prager criteria for the elasto-plastic
% material law. It provides methods for evaluating stress, constitutive
% matrices, yield conditions, flow vectors, and their gradients, as well
% as handling plastic strain updates.
%
%% Methods
% * *eval*: Computes the stress vector and the constitutive matrix for the 
%           material at a given integration point. Handles both elastic 
%           and plastic steps.
% * *yieldCondition*: Defines the yield function based on the Drucker-
%                     Prager criteria.
% * *yieldStressGradient*: Computes the gradient of the yield function 
%                          with respect to the stress vector.
% * *flowVector*: Computes the flow vector for the plastic potential.
% * *flowVectorGradient*: Computes the gradient of the flow vector with 
%                         respect to the stress vector.
% * *hardening*: Returns the hardening value.
% * *hardeningStressGradient*: Returns the gradient of the hardening law 
%                              with respect to the stress vector  
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalElastoPlasticDruckerPrager < MechanicalElastoPlastic  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalElastoPlasticDruckerPrager()
            this = this@MechanicalElastoPlastic();
            this.nstVar = 5;   % Hardening + Kinematic hardening
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = eval(this,material,ip)

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);
            Ce = this.elasticFlexibilityMatrix(material,ip);
            Dt = De;

            % Trial stress vector
            stress = De * (ip.strain - ip.plasticstrain);

            % Evaluate the yield condition
            f = this.yieldCondition(material,ip,stress);

            % Elastic step
            if f < 0.0, return, end

            % Material parameters
            tanPhi = tan(material.frictionAngle);
            tanPsi = tan(material.dilationAngle);
            coh    = material.cohesion;
            Id     = this.gradientI1(stress);

            % Elastic properties
            K = this.bulkModulus(material);
            G = this.shearModulus(material);

            % Stress invariants
            p = this.hydrostaticStress(stress);
            q = this.vonMisesStress(stress);

            % Deviatoric stresses
            s = this.deviatoricStress(stress);

            % Plastic multiplier
            lambda = (q + p * tanPhi - coh) / (3.0 * G + K * tanPhi * tanPsi);

            % Stress update
            factor = 1.0 - 3.0 * G * lambda / q;
            if factor >= 0.0
                s = factor * s;
                p = p - lambda * K * tanPsi;
                stress = s + p * Id;
                df = this.yieldStressGradient(material,ip,stress);
                n  = this.flowVector(material,ip,stress);
                dn = this.flowVectorGradient(material,ip,stress);
                Psi = inv(Ce + lambda * dn);
                Dt  = Psi - (Psi * (n * df') * Psi)/(df' * (Psi * n));
            else 
                % Return to the apex of the surface
                stress = (coh / tanPhi) * Id;
                Dt = 1.0e-8*eye(4,4);
            end

            % Update the plastic strain
            ip.plasticstrain = ip.strain - Ce * stress;
        end

        %------------------------------------------------------------------
        % Yield function definition
        function f = yieldCondition(this,material,~,stress)
            % Material parameters
            coh = material.cohesion;
            phi = material.frictionAngle;
            % Stress invariants
            I1 = this.stressInvariantI1(stress);
            J2 = max(this.stressInvariantJ2(stress),1.0e-8);
            % Yield surface
            f  = tan(phi) * I1 / 3.0 + sqrt(3.0 * J2) - coh;
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the stress vector
        function df = yieldStressGradient(this,material,ip,stress)
            % Material parameters
            phi = material.frictionAngle;
            % Stress invariants gradients
            dI1 = this.gradientI1(stress);
            dJ2 = this.gradientJ2(stress);
            J2 = max(this.stressInvariantJ2(stress),1.0e-8);
            % Derivatives of the yield surface wrt to the invariants
            dfdI1 = tan(phi) / 3.0;
            dfdJ2 = 0.5 * sqrt(3.0 / J2);
            % Yield surface gradient
            df = dfdI1 * dI1 + dfdJ2 * dJ2;
            if strcmp(ip.anm,'PlaneStress')
                df(3) = 0.0;
            end
        end

        %------------------------------------------------------------------
        % Flow vector
        function n = flowVector(this,material,ip,stress)
            % Material parameters
            psi = material.dilationAngle;
            % Deviatoric stress invariant
            J2 = max(this.stressInvariantJ2(stress),1.0e-8);
            % Stress invariants gradients
            dI1 = this.gradientI1(stress);
            dJ2 = this.gradientJ2(stress);
            % Derivatives of the yield surface wrt to the invariants
            dfdI1 = tan(psi) / 3.0;
            dfdJ2 = 0.5 * sqrt(3.0 / J2);
            % Yield surface gradient
            n = dfdI1 * dI1 + dfdJ2 * dJ2;
            if strcmp(ip.anm,'PlaneStress')
                n(3) = 0.0;
            end
        end

        %------------------------------------------------------------------
        % Flow vector gradient
        function dn = flowVectorGradient(this,~,ip,stress)
            % Deviatoric stress invariant 
            J2   = max(this.stressInvariantJ2(stress),1.0e-8);
            dJ2  = this.gradientJ2(stress);
            d2J2 = this.hessianJ2();
            % Derivatives of the yield surface wrt to the invariants
            dfdJ2 = 0.5 * sqrt(3.0 / J2);
            d2fdJ2 = -0.25 * sqrt(3.0 / J2 / J2 / J2);
            % Yield surface gradient
            dn = d2fdJ2 * (dJ2 * dJ2') + dfdJ2 * d2J2;
            if strcmp(ip.anm,'PlaneStress')
                dn(3,:) = 0.0;
                dn(:,3) = 0.0;
                dn(3,3) = 1.0;
            end
        end

        %------------------------------------------------------------------
        % Hardening law
        function h = hardening(~,~,~,~)
            h = 0.0;
        end

        %------------------------------------------------------------------
        % Gradient of the hardening law wrt to the stress vector
        function dh = hardeningStressGradient(~,~,~,~)
            dh = 0.0;
        end

    end
end