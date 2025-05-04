%% NonlinearScheme_Newton Class
% This class implements a fully implicit time integration scheme using the
% Newton-Raphson method for solving nonlinear systems. It inherits from the
% _NonlinearScheme_ base class and provides methods for assembling the 
% linear system, applying boundary conditions, adding nodal forces, 
% evaluating the solution, and checking for convergence.
%
%% Methods
% * *assembleLinearSystem*: Assembles the Jacobian matrix and the residual 
%                           vector for the nonlinear system.
% * *applyBCtoRHS*: Applies boundary conditions to the right-hand side 
%                   vector.
% * *addNodalForces*: Adds nodal forces to the right-hand side vector. 
% * *eval*: Evaluates the solution increment and updates the solution 
%           vector.
% * *convergence*: Checks for convergence of the nonlinear solver.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
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
        % Assembles the linear system
        function [J,r] = assembleLinearSystem(~,C,K,fi,fe,dfidx,x,xOld,dt)
            % Compute residual vector
            r = fi + K * x + C * (x - xOld) / dt - fe;

            % Jacobian matrix
            J = K + dfidx + C / dt;
        end

        %------------------------------------------------------------------
        % Applies the boundary conditions to the right-hand side vector
        function bf = applyBCtoRHS(~,~,b,~,doffree,~)
            bf = b(doffree);
        end

        %------------------------------------------------------------------
        % Adds nodal forces to the right-hand side vector
        function b = addNodalForces(~,b,fe)
            b = b - fe;
        end

        %------------------------------------------------------------------
        % Evaluates the solution increment and updates the solution vector
        function [X,dx] = eval(~,J,r,X,~,freedof,~)
            % Compute increment of variables
            dx = -J\r;

            % Update variables
            X(freedof) = X(freedof) + dx;
        end
        
        %------------------------------------------------------------------
        % Checks the convergence of the nonlinear solver
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
