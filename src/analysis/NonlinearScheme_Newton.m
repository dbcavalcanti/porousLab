%% NonlinearScheme_Newton Class
%
% Considers a fully implicit time integration scheme.
%
classdef NonlinearScheme_Newton < NonlinearScheme 
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = NonlinearScheme_Newton()
            this = this@NonlinearScheme();
        end
    end
    %% Public methods
    % Implementation of the abstract methods declared in super-class
    methods (Static)

        %------------------------------------------------------------------
        function [J, r] = assembleLinearSystem(C, K, fi, fe, dfidx, x, xOld, dt)
            % Compute residual vector
            r = fi + K * x + C * (x - xOld) / dt - fe;
            % Jacobian matrix
            J = K + dfidx + C / dt;
        end
        
        %------------------------------------------------------------------
        function [J, r] = assemblePerturbedLinearSystem(C, K, fi, fe, dfidx, x, xOld, dt, ngle)
            
            % Compute residual vector (non perturbed terms)
            r = fi(:,2) + K(:,ngle+1:2*ngle) * x + C(:,ngle+1:2*ngle) * (x - xOld) / dt - fe(:,2);
            
            % Non perturbed Jacobian
            nonPerturbedJ = K(:,ngle+1:2*ngle) + dfidx(:,ngle+1:2*ngle) + (C(:,ngle+1:2*ngle) / dt);
            
            % Perturbed Jacobian
            % Notice that, as the assemly is done element by element (and
            % not for DOF), the perturbation is taken as middle value
            % between the displacement perturbation (1e-8) and the pressure
            % perturbation (1e-3)
            
            % Needed factors
            eps = 1e-5;
            x_dot = (x - xOld)/dt;
            
            % Perturbed terms
            perturbedC = (C(:,2*ngle+1:end) - C(:,1:ngle))/(2*eps) * x_dot;
            perturbedK = (K(:,2*ngle+1:end) - K(:,1:ngle))/(2*eps) * x;
            perturbedfi = (fi(:,3) - fi(:,1))/(2*eps);

            % Final Jacobian assembly
            J = nonPerturbedJ + perturbedC + perturbedK - perturbedfi;

            % TODO. Check how OGS fills the Jacobian as they loop for every
            % DOF and the resultant column vector is added at the
            % corresponding jacobian column.
            % For example, horizontal displacement at 2n node correspond to
            % column 5, thus, the perturbed terms will be added directly to
            % column 5

        end 
        
        %------------------------------------------------------------------
        function bf = applyBCtoRHS(~, b, ~, doffree, ~)
            bf = b(doffree);
        end

        %------------------------------------------------------------------
        function b = addNodalForces(b,fe)
            b = b - fe;
        end

    end

    methods

        %------------------------------------------------------------------
        function [X, dx] = eval(~,J,r,X,~,freedof,~)
            % Compute the increment of the variables
            dx = -J\r;
            % Update the variables
            X(freedof) = X(freedof) + dx;
        end

        %------------------------------------------------------------------
        function convFlg = convergence(this,~,~,~,r,~,iter)
            fprintf("\t\t iter.: %3d , ||R|| = %7.3e \n",iter,norm(r));
            if (norm(r) < this.tol) && (iter > 1)
                convFlg = true;
            else
                convFlg = false;
            end
        end

    end
    
end