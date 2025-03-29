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
    methods
        %------------------------------------------------------------------
        function [J,r] = assembleLinearSystem(~,C,K,fi,fe,dfidx,x,xOld,dt)
            % Compute residual vector
            r = fi + K * x + C * (x - xOld) / dt - fe;

            % Jacobian matrix
            J = K + dfidx + C / dt;
        end

        %------------------------------------------------------------------
        function bf = applyBCtoRHS(~,~,b,~,doffree,~)
            bf = b(doffree);
        end

        %------------------------------------------------------------------
        function b = addNodalForces(~,b,fe)
            b = b - fe;
        end

        %------------------------------------------------------------------
        function [X,dx] = eval(~,J,r,X,~,freedof,~)
            % Compute increment of variables
            dx = -J\r;

            % Update variables
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
