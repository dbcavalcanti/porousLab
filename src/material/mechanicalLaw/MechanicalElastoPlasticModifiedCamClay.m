%% MechanicalElastoPlasticModifiedCamClay Class
% This class implements the Modified Cam Clay criteria for the 
% elasto-plastic material law. It provides methods for evaluating stress, 
% constitutive matrices, yield conditions, flow vectors, and their 
% gradients, as well as handling plastic strain updates.
%
%% Methods
% * *eval*: Computes the stress vector and the constitutive matrix for the 
%           material at a given integration point. Handles both elastic 
%           and plastic steps.
% * *alternativeStressIntegration*: Implements an alternative stress 
%                                   integration algorithm for the material.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef MechanicalElastoPlasticModifiedCamClay < MechanicalElastoPlastic  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MechanicalElastoPlasticModifiedCamClay()
            this = this@MechanicalElastoPlastic();
            % State variables: [pc ; ]
            % where:
            %   pc: pre-consolidation pressure
            this.nstVar = 2;
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = eval(this,material,ip)
            if strcmp(material.stressIntAlgorithm,'implicit')
                 [stress,Dt] = eval@MechanicalElastoPlastic(this,material,ip);
            elseif strcmp(material.stressIntAlgorithm,'alternative')
                [stress,Dt] = this.alternativeStressIntegration(material,ip);
            else
                disp('Error: the given stress integration algorithm is not available');
                disp('Tags of the methods available: ''implicit'', ''alternative''');
                error('Error: stressIntAlgorithm is not available');
            end
        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive matrix
        function [stress,Dt] = alternativeStressIntegration(this,material,ip)

        end

        %------------------------------------------------------------------
        % Yield function definition
        function f = yieldCondition(this,material,~,stress,state)
            M = material.mccCriticalStateSlope;
            p = this.pressureStress(stress);
            q = this.vonMisesStress(stress);
            pc = state(1);
            f = q * q - M * M * p * (p - pc);
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the stress vector
        function df = yieldStressGradient(this,material,~,stress,state)
            % Get material parameter
            M = material.mccCriticalStateSlope;
            % Get the pre-consolidation pressure
            pc = state(1);
            % Compute the pressure and the von Mises stresses
            p = this.pressureStress(stress);
            q = this.vonMisesStress(stress);
            % Derivatives of the yield function wrt p and q
            dfdp = M * M * (2.0 * p - pc);
            dfdq = 2.0 * q;
            % Derivatives of p and q wrt to the stress tensor (Voigt)
            dpdstress = this.gradientPressureStress();
            dqdstress = this.vonMisesStressGradient(stress);
            % Derivative of the yield function
            df = dfdp * dpdstress + dfdq * dqdstress;
        end

        %------------------------------------------------------------------
        % Gradient of the yield function wrt to the state variables vector
        % The state variable here is the pre-consolidation pressure.
        function dfda = yieldStateGradient(this,material,~,stress,~)
            % Get material parameter
            M = material.mccCriticalStateSlope;
            % Compute the pressure stresses
            p = this.pressureStress(stress);
            % Derivative of the yield function wrt pre-consolidation
            % pressure
            dfda = -M*M*p;
        end

        %------------------------------------------------------------------
        % Flow vector (associative plasticity)
        function n = flowVector(this,material,ip,stress,state)
            n = this.yieldStressGradient(material,ip,stress,state);
        end

        %------------------------------------------------------------------
        % Flow vector gradient
        % The flow vector can be defined in a more simplified way as:
        %           n = 3 * s + M^2 * (2*p - pc) * dp
        % where: 
        %   s: deviatoric stress tensor in Voigt notation as a vector,
        %       which is equal to the gradient of J2.
        %   dp: gradient of the pressure stress
        function dn = flowStressGradient(this,material,~,~,~)
            % Get material parameter
            M = material.mccCriticalStateSlope;
            % Get auxiliar derivatives
            dp = this.gradientPressureStress();
            ds = this.hessianJ2();
            % Gradient of the flow vector
            dn = 3.0 * ds + 2.0 * M * M * (dp * dp');
        end

        %------------------------------------------------------------------
        % Flow vector gradient wrt to the state variables vector
        function dnda = flowStateGradient(this,material,~,~,~)
            % Get material parameter
            M = material.mccCriticalStateSlope;
            % Get auxiliar derivatives
            dp = this.gradientPressureStress();
            % Derivative wrt the pre-consolidation pressure
            dnda = - M * M * dp;
        end

        %------------------------------------------------------------------
        % Internal/state variables evolution
        function h = stateEvolution(~,~,~,~,~)

        end

        %------------------------------------------------------------------
        % Gradient of the internal/state variables law wrt to the stress vector
        function dhds = stateStressGradient(~,~,ip,~,~)

        end

        %------------------------------------------------------------------
        % Gradient of the internal/state variables law wrt to the state variables
        function dhda = stateStateGradient(~,~,~,~,~)

        end

    end
end
