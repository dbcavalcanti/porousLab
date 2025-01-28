%% MechanicalIsotropicDamage class
%
% This class defines a isotropic damage material model
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef MechanicalIsotropicDamage < MechanicalLinearElastic  
    properties (SetAccess = public, GetAccess = public)
        lc = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalIsotropicDamage(lc)
            this = this@MechanicalLinearElastic();
            this.nstVar = 2;  % Scalar damage + Damage threshold
            this.lc = lc;     % Element characteristic length
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = eval(this,material,ip)

            % Constitutive matrix
            De = this.elasticConstitutiveMatrix(material,ip);

            % Effective stress vector
            effstress = De * (ip.strain - ip.strainOld) + ip.stressOld;

            % Update the damage threshold
            DrDstrain = this.damageCriteria(material,ip);

            % Damage law
            DdDr = this.damageLaw(material,ip);

            % Tangent constitutive matrix
            Dt = this.tangentConstitutiveMatrix(ip,effstress,De,DrDstrain,DdDr);

            % Update the stress vector
            stress = (1.0 - ip.statevar(1)) * effstress;

        end

        %------------------------------------------------------------------
        % von Mises equivalent strain and its derivative wrt to the strain
        % tensor
        function [eqstrain,Deqstrain] = equivalentStrain(this,material,ip)
            % Get the material parameters
            k = material.kappa;
            nu = material.nu;
            % Compute out-of-plane strain (plane stress problems)
            this.elasticOutOfPlaneStrain(material,ip)
            % Compute the strain invariants
            I1 = this.strainInvariantI1(ip.strain);
            J2 = this.strainInvariantJ2(ip.strain);
            % Compute the equivalent strain
            eqstrain = (k - 1.0) / (2.0 * k * (1.0 - 2.0 * nu)) * I1;
            eqstrain = eqstrain + sqrt(((k - 1.0) / (1.0 - 2.0 * nu) * I1) ^ 2.0 + 12.0 * k * J2 / ((1.0 + nu) * (1.0 + nu))) / (2.0 * k);
            % Get the gradient of the invariants
            dI1 = this.gradientStrainInvariantI1();
            dJ2 = this.gradientStrainInvariantJ2(ip.strain);
            I4  = this.fourthOrderSymTensor();
            % Derivative of equivalent strain wrt to strain invariants
            DqstrainDI1 = I1 * (k - 1.0) * (k - 1.0) / (2.0 * k * (2.0 * nu - 1.0) * (2.0 * nu - 1.0) * sqrt(12.0 * J2 * k / ((nu + 1.0) * (nu + 1.0)) + I1 * I1 * (k - 1.0) * (k - 1.0) / ((2.0 * nu - 1) * (2.0 * nu - 1)))) - (k - 1.0) / (2.0 * k * (2.0 * nu - 1.0));
            DqstrainDJ2 = 3.0 / (sqrt(12.0 * J2*k / ((nu + 1.0) * (nu + 1.0)) + I1 * I1 * (k - 1.0) * (k - 1.0) / ((2.0 * nu - 1.0) * (2.0 * nu - 1.0))) * (nu + 1) * (nu + 1));
            % Derivate of the equivalent strain wrt to the strain tensor
            Deqstrain = I4 * (DqstrainDI1 * dI1 + DqstrainDJ2 * dJ2);
        end

        %------------------------------------------------------------------
        % Damage criteria
        % The state variable r is the maximum equivalente strain value in
        % the entire load history
        function DrDstrain = damageCriteria(this,material,ip)
            [eqstrain,Deqstrain] = this.equivalentStrain(material,ip);
            ip.statevarOld(2) = max(ip.statevarOld(2),material.DamageThreshold);
            if eqstrain > ip.statevarOld(2)
                ip.statevar(2) = eqstrain;
                DrDstrain = Deqstrain;
            else
                ip.statevar(2) = ip.statevarOld(2);
                DrDstrain = zeros(4,1);
            end
        end

        %------------------------------------------------------------------
        % Damage exponential law
        function DdDr = damageLaw(this,material,ip)
            % Get the material parameters
            E    = material.Young;
            r0   = material.DamageThreshold;
            Gf   = material.FractureEnergyMode1;
            beta = E * r0 * this.lc / Gf;
            % Scalar damage
            r = ip.statevar(2);
            d = 1.0 - r0 / r * exp(-beta * (r - r0));
            if d < eps
                d = eps;
            elseif (d > (1.0 - eps))
                d = 1.0 - eps;
            end
            % Update value at the integration point
            ip.statevar(1) = d;
            % Compute the derivative of the damage wrt damage threshold
            DdDr = r0 / r * (1.0 / r + beta) * exp(-beta * (r - r0));
        end

        %------------------------------------------------------------------
        % Damage exponential law
        function Dt = tangentConstitutiveMatrix(this,ip,effstress,De,DrDstrain,DdDr)
            d = ip.statevar(1);
            Dt = (1.0 - d) * De - DdDr * effstress * DrDstrain';
            if strcmp(ip.anm,'PlaneStress')
                Dt = this.planeStressConstitutiveMatrix(Dt);
            end
        end

    end
    %%
    methods (Static)
        function flag = isElastoPlastic()
            flag = false;
        end
    end
end