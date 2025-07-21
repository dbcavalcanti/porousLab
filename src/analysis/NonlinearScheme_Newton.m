%% NonlinearScheme_Newton Class
% This class inherits from the base class 'NonlinearScheme' to implement
% a fully implicit time integration scheme using the Newton-Raphson method for solving nonlinear systems.
%
%% Authors
% * Danilo Cavalcanti
% 
%% Class definition
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
        % Assemble the Jacobian matrix and the residual vector.
        function [J,r] = assembleLinearSystem(~,C,K,fi,fe,dfidx,x,xOld,dt)
            % Residual vector
            r = fi + K * x + C * (x - xOld) / dt - fe;

            % Jacobian matrix
            J = K + dfidx + C / dt;
        end

        %------------------------------------------------------------------
        % Apply boundary conditions to the right-hand side of the system.
        function bf = applyBCtoRHS(~,~,b,~,doffree,~)
            bf = b(doffree);
        end

        %------------------------------------------------------------------
        % Add nodal forces to the right-hand side vector.
        function b = addNodalForces(~,b,fe)
            b = b - fe;
        end

        %------------------------------------------------------------------
        % Evaluate the solution increment and updates the solution vector.
        function [X,dx] = eval(~,J,r,X,~,freedof,~)
            % Compute increment of variables
            dx = -J\r;

            % Update variables
            X(freedof) = X(freedof) + dx;
        end

        %------------------------------------------------------------------
        % Check for convergence of the nonlinear scheme.
        function convFlg = convergence(this,~,XOld,dx,r,~,iter,echo)
            if echo
                fprintf("\t\t iter.: %3d , ||R|| = %7.3e  , ||dx||/||X0|| = %7.3e \n",iter,norm(r),norm(dx)/norm(XOld));
            end
            if ((norm(r) < this.tol) || (norm(dx)/norm(XOld)) < this.tol) && (iter > 1)
                convFlg = true;
            else
                convFlg = false;
            end
        end
    end
end
