%% Anl_Transient class
%
classdef Anl_Transient < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        theta       = 1.0;    % Implicit time integration scheme parameter
        ti          = 0.01;   % Initial time
        tf          = 1.0;    % Final time
        dt          = 0.001;  % Time step
        dtMax       = 0.001;  % Maximum time step
        dtMin       = 0.001;  % Minimum time step
        adaptStep   = false;  % Adaptive step size
        maxIter     = 250;    % Maximum number of iterations
        maxAttempts = 10;     % Maximum attempts to converge
        nlscheme    = [];     % Nonlinear solution scheme
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Transient(nlscheme)
            this = this@Anl('Transient');

            if strcmp(nlscheme,'Picard')
                this.nlscheme = NonlinearScheme_Picard();
            elseif strcmp(nlscheme,'Newton')
                this.nlscheme = NonlinearScheme_Newton();
            else
                disp("Error creating the Analysis object.");
                disp("Nonlinear solution scheme was not provided.");
                disp("Available options:");
                disp("   Picard");
                disp("   Newton");
                error("Error: Nonlinear solution scheme was not provided");
            end
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        function run(this,mdl)
            % Initialize analysis parameters
            t    = this.ti;
            t0   = this.ti;
            step = 1;

            % Initialize model object
            mdl.preComputations();
            
            disp("*** Performing transient nonlinear analysis...")

            % Initialize solution vector
            X  = mdl.U;
            dx = zeros(mdl.ndof);

            % Start transient process
            attempt = 1;
            brokenStep = false;
            while (t0 < this.tf)
                fprintf("\t Time: %12.5f s \n", t);

                % Update transient solution
                XOld = X;

                % Start iterative process
                convFlg = false;
                attemptOld = attempt;
                attempt = 1;

                while attempt < this.maxAttempts
                    iter = 1;
                    
                    while true
                        % Compute model global matrices
                        [A,b] = mdl.getLinearSystem(X,XOld,this.nlscheme,this.dt);

                        % Apply Dirichlet boundary conditions
                        [A,b] = mdl.applyDirichletBC(A,b,X,this.nlscheme);

                        % Update variables
                        [X,dx] = this.nlscheme.eval(A,b,X,dx,mdl.doffree,iter);

                        % Check convergence
                        convFlg = this.nlscheme.convergence(X,XOld,dx,b,mdl.doffree,iter);
                        if convFlg == true
                            break;
                        end

                        % Check maximum number of iterations
                        iter = iter + 1;
                        if (iter > this.maxIter)
                            break
                        end
                    end

                    % Check convergence
                    if convFlg == true
                        break;
                    end

                    % Reduce time step 
                    this.dt = max(this.dt/4.0, this.dtMin);

                    % Clean previous attempt
                    X = XOld;

                    % Update attempt counter
                    attempt = attempt + 1;
                end

                if convFlg == false
                    disp("Solution did not converge!");
                    break;
                end

                % Update state variables
                mdl.updateStateVar();

                % Update time step
                if (this.adaptStep == true) && (attempt == 1) && (brokenStep == false) && (attemptOld == 1)
                    this.dt = min(2 * this.dt, this.dtMax);
                end
                
                % Update time
                t0 = t;
                if (t + this.dt) > this.tf
                    this.dt = this.tf - t;
                end
                t = t + this.dt;
                step = step + 1;
            end

            % Update state variables
            mdl.updateStateVar();

            % Save final result
            mdl.U = X;
            for i = 1:mdl.nelem
                gle = mdl.element(i).type.gle;
                mdl.element(i).type.ue = mdl.U(gle);
            end

            disp("*** Analysis completed!");
        end

        %------------------------------------------------------------------
        function setPicardRelaxation(this,flag)
            this.nlscheme.applyRelaxation = flag;
        end

        %------------------------------------------------------------------
        function setRelativeConvergenceCriteria(this,flag)
            this.nlscheme.normalizeError = flag;
        end

        %------------------------------------------------------------------
        function printStep(~,X,mdl)
            for i = 1:mdl.nnodes
                fprintf("  %4d: \t", i);
                for j = 1:mdl.ndof_nd
                    fprintf("  %8.4f ", X(mdl.ID(i,j)));
                end
                fprintf("\n");
            end
        end

        %------------------------------------------------------------------
        function [normRp,normRpg] = decomposeResidualVct(~,r,Fext,mdl)
            rp = r(mdl.pFreeDof);
            rpg = r(mdl.pgFreeDof);

            nfep = norm(Fext(mdl.pFreeDof));
            nfepg = norm(Fext(mdl.pgFreeDof));

            if nfep < 1.0e-15
                nfep = 1.0;
            end
            if nfepg < 1.0e-15
                nfep = 1.0;
            end

            normRp = norm(rp)/nfep;
            normRpg = norm(rpg)/nfepg; 
        end

        %------------------------------------------------------------------
        function setUpTransientSolver(this,ti,dt,tf,dtMax,dtMin,adaptStep)
            if nargin == 4
                dtMax = dt;
                dtMin = dt;
                adaptStep = false;
            end
            this.ti = ti;
            this.dt = dt;
            this.tf = tf;
            this.adaptStep = adaptStep;
            this.dtMax = dtMax;
            this.dtMin = dtMin;
        end
    end
end
