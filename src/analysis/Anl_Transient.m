%% Anl_Transient Class
%
classdef Anl_Transient < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        theta = 1.0;    % Implicit time integration scheme
        dt    = 0.001;  % Time increment;
        tf    = 1.0;    % Final time
        tinit = 0.01;
        dtMax = 0.001;
        dtMin = 0.001;
        adaptStep = false;
        maxIter = 250;
        maxAttempts = 10;
        nlscheme = [];
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Transient(result,nlscheme)
            this = this@Anl('Transient',result);
            if strcmp(nlscheme,'Picard')
                this.nlscheme = NonlinearScheme_Picard();
            elseif strcmp(nlscheme,'Newton')
                this.nlscheme = NonlinearScheme_Newton();
            else
                disp("Error creating the Analysis object...");
                disp("Non available nonlinear solution scheme...");
                disp("Nonlinear solution schemes availables:");
                disp("   Picard");
                disp("   Newton");
                error("Error: Non available nonlinear solution scheme");
            end
            
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function process(this,mdl)

            disp("*** Initialize nonlinear transient analysis...")

            % Initialize transient analysis parameters
            t        = this.tinit;
            t0       = this.tinit;
            maxSteps = length(this.tinit:this.dt:this.tf);
            step     = 1;

            % Initialize the model object
            mdl.preComputations();

            % Initialize solution vector
            X  = mdl.U;
            dx = zeros(mdl.ndof);

            % Initialize the output vectors
            this.result.time = zeros(maxSteps,1);
            this.result.p    = zeros(maxSteps,1);
            attempt = 1;
            brokenStep = false;

            % Transient analysis
            while (t0 < this.tf)
                fprintf("\t Time: %12.5f s \n",t);

                % Update transient solution
                XOld = X;

                % Iterative solution
                convFlg = false;
                attemptOld = attempt;
                attempt = 1;
                
                while attempt < this.maxAttempts

                    iter = 1;
                    % Iterative process
                    while true
    
                        % Compute model global matrices
                        [A,b] = mdl.getLinearSystem(X,XOld,this.nlscheme,this.dt);

                        % Apply Dirichlet boundary conditions
                        [A,b] = mdl.applyDirichletBC(A, b, X, this.nlscheme);

                        % Update variables
                        [X, dx] = this.nlscheme.eval(A,b,X,dx,mdl.doffree,iter);
    
                        % Check convergence
                        convFlg = this.nlscheme.convergence(X,XOld,dx,b,mdl.doffree,iter);
                        if convFlg == true, break;  end
    
                        % Check maximum number of iterations
                        iter = iter + 1;
                        if (iter > this.maxIter)
                            break
                        end
                    end

                    % Check convergence
                    if convFlg == true, break, end

                    % Reduce the time step 
                    this.dt = max(this.dt / 4.0, this.dtMin);
                    
                    % Clean the previous attempt
                    X = XOld;

                    % Update counter variable
                    attempt = attempt + 1;
                end

                if convFlg == false
                    disp("Transient analysis: solution did not converge!");
                    break
                end

                % Update the state variables
                mdl.updateStateVar();

                % Update time step
                if (this.adaptStep == true) && (attempt == 1) && (brokenStep == false) && (attemptOld == 1)
                    this.dt = min(this.dt * 2,this.dtMax);
                end

                % Save results
                this.result.time(step) = t;
                if isempty(this.result.coordP) == false
                    elemPlot = findElementInMesh(mdl.NODE, mdl.ELEM,this.result.coordP);
                    gle = mdl.element(elemPlot).type.gle;
                    this.result.p(step) = mdl.element(elemPlot).type.pressureField(this.result.coordP,X(gle));
                else
                    this.result.p(step) = 0.0;
                end  
                
                % Update time
                t0 = t;
                if (t + this.dt) > this.tf
                    this.dt = this.tf - t;
                end
                t = t + this.dt;
                step = step + 1;

            end
            
            mdl.updateStateVar();

            % Save final result
            mdl.U = X;
            for el = 1:mdl.nelem
                gle = mdl.element(el).type.gle;
                mdl.element(el).type.ue = mdl.U(gle);
            end

            % Remove unused steps in the output vectors
            this.result.p = this.result.p(1:(step-1));
            % this.result.ST(:) = this.result.ST(:,1:(step-1));
            this.result.time = this.result.time(1:(step-1));

            disp("*** Analysis completed.")

        end

        function setPicardRelaxation(this,flag)
            this.nlscheme.applyRelaxation = flag;
        end

        function setRelativeConvergenceCriteria(this,flag)
            this.nlscheme.normalizeError = flag;
        end

        function printStep(~,X,mdl)

            for i = 1:mdl.nnodes
                fprintf("  %4d: \t",i);
                for j = 1:mdl.ndof_nd
                    fprintf("  %8.4f ",X(mdl.ID(i,j)))
                end
                fprintf("\n");
            end

        end


        function [normRp , normRpg] = decomposeResidualVct(~,r, Fext,mdl)

            rp = r(mdl.pFreeDof);
            rpg = r(mdl.pgFreeDof);

            nfep = norm(Fext(mdl.pFreeDof));
            nfepg = norm(Fext(mdl.pgFreeDof));

            if nfep < 1.0e-15, nfep = 1.0; end
            if nfepg < 1.0e-15, nfep = 1.0; end

            normRp = norm(rp)/nfep;
            normRpg = norm(rpg)/nfepg; 

        end

        function setUpTransientSolver(this,tinit,dt,tf,dtMax,dtMin,adaptStep)
            if nargin == 4
                dtMax = dt;
                dtMin = dt;
                adaptStep = false;
            end
            this.tinit = tinit;
            this.dt = dt;
            this.tf = tf;
            this.adaptStep = adaptStep;
            this.dtMax = dtMax;
            this.dtMin = dtMin;
        end
    end
    
end