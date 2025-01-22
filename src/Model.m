%% Model class
%
% This class defines a finite element model that has a strong discontinuity
%
%% Author
% Danilo Cavalcanti
%
%
%% Class definition
classdef Model < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        physics             = [];            % Physics of the problem
        NODE                = [];            % Nodes of the fem mesh
        ELEM                = [];            % Nodes connectivity
        t                   = 1.0;           % Thickness
        mat                 = [];            % Vector with material properties
        type                = 'ISOQ4';       % Typf of element used
        intOrder            = 2;             % Number of integration points
        nnodes              = 1;             % Number of nodes
        nelem               = 1;             % Number of elements
        nnd_el              = 4;             % Number of nodes per element
        doffree             = [];
        doffixed            = [];
        ndof_nd             = 2;             % Number of dof per node
        ndof                = 1;             % Number of regular degrees of freedom
        ndoffree            = 0;             % Number of free degrees of freedom
        ndoffixed           = 0;             % Number of fixed degrees of freedom
        Dof                 = [];            % Vector with all regular dofs
        ID                  = [];            % Each line of the ID matrix contains the global numbers for the node DOFs (DX, DY)
        F                   = [];            % Global force vector
        U                   = [];            % Global displacement vector
        element             = [];            % Array with the element's objects
        nDofElemTot         = 0.0;
        sqrNDofElemTot      = 0.0;
        matID               = [];            % Vector with the material id of each element
        massLumping         = false;
        lumpStrategy        = 1;
        isAxisSymmetric     = false;
    end
    
    %% Constructor method
    methods
        function this = Model()
        end
    end

    %% Abstract methods
    methods(Abstract)
        createNodeDofIdMtrx(this);
        assembleElementDofs(this);
        initializeElements(this);
        initializeDisplacementVct(this);
        addNodalLoad(this,Fref);
        updateResultVertexData(this,type);
        printResultsHeader();
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function preComputations(this)
            
            % Initialize basic variables
            this.initializeBasicVariables();

            % Create nodes DOF ids matrix
            this.createNodeDofIdMtrx();

            % Assemble the regular dofs to each element
            this.assembleElementDofs();

            % Initialize elements
            this.initializeElements();
            
            % Compute auxiliar variables for assemblage of sparse matrices
            this.initializeSparseMtrxAssemblageVariables();

            % Initialize the displacement vector
            this.initializeDisplacementVct();
            
            % Initialize the external force vector
            this.initializeExtForceVct();

        end

        %------------------------------------------------------------------
        function initializeBasicVariables(this)
            this.nnodes   = size(this.NODE,1);
            this.nelem    = size(this.ELEM,1);     
            this.nnd_el   = size(this.ELEM,2);            
            this.ndof     = this.ndof_nd * this.nnodes; 
            if isempty(this.matID)
                this.matID  = ones(this.nelem,1);
            end
        end

        %------------------------------------------------------------------
        function Lce = getElementsCharacteristicLength(this)
            Lce=zeros(size(this.ELEM,1),1);
            for el = 1:size(this.ELEM,1)
            
                % Vertices of the element el coordinates
                vx = this.NODE(this.ELEM(el,:),1); 
                vy = this.NODE(this.ELEM(el,:),2);
            
                % Number of vertices 
                nv = length(this.ELEM(el,:)); 
            
                % Shifted vertices
                vxS = vx([2:nv 1]);
                vyS = vy([2:nv 1]); 
            
                % Compute the area of the element (trapezoidal rule)
                temp = vx.*vyS - vy.*vxS;
                Ae   = 0.5*sum(temp);
                
                % Characteristic lenght (quadrilateral elements)
                Lce(el) = sqrt(Ae);
            end
        end

        %------------------------------------------------------------------
        function initializeSparseMtrxAssemblageVariables(this)
            this.nDofElemTot = 0;
            this.sqrNDofElemTot = 0;
            for el = 1:this.nelem
                this.nDofElemTot = this.nDofElemTot + length(this.element(el).type.gle);
                this.sqrNDofElemTot = this.sqrNDofElemTot + length(this.element(el).type.gle)*length(this.element(el).type.gle);
            end
        end
        
        %------------------------------------------------------------------
        function initializeExtForceVct(this)
            this.F = zeros(this.ndof,1);
        end

        %------------------------------------------------------------------
        % Compute the jacobian matrix and the residual vector associated
        % with the coupled system discretized using an implicit
        % time-integration scheme.
        function [J, r] = setTransientSystem(this, dt, X, DX, Fext)   

            % Compute model global matrices
            [K, C, Fint] = this.globalMatrices(DX);

            % Get components associated with the free dofs
            freedof = 1:this.ndoffree;
            pfreedof = this.pFreeDof;

            % Compute the Jacobian matrix
            J =  K(freedof,freedof) + C(freedof,freedof) / dt;

            % Compute the residual vector: 
            % r = Fint + Q * P - Fext + C * DX / dt
            r = Fint(freedof) + K(freedof,pfreedof)*X(pfreedof) - Fext(freedof) + C(freedof,freedof) * DX(freedof) / dt;

        end

        %------------------------------------------------------------------
        % Compute the matrix and vector of the problem applying and
        % implicit time discretization
        function [A,b] = setLinearSystem(this, X, Fext)   

            % Compute model global matrices
            K = this.globalMatrices(X);

            % Get components associated with the free dofs
            freedof  = 1:this.ndoffree;
            fixeddof = [1:this.ndof, freedof];

            % Compute the right-handside matrix
            A =  K(freedof,freedof);

            % Compute the left-handside vector
            b =  Fext(freedof) - K(freedof,fixeddof)*X(fixeddof);

        end

        %------------------------------------------------------------------
        % Global system matrices
        function [K, C, Fi, Fe] = globalMatrices(this,U)   

            % Indices for the assembling the matrices and vector
            iDof = zeros(this.sqrNDofElemTot,1);
            jDof = zeros(this.sqrNDofElemTot,1);
            eDof = zeros(this.nDofElemTot,1);

            % Initialize the components of the matrices and vector
            K_ij  = zeros(this.sqrNDofElemTot,1);
            C_ij  = zeros(this.sqrNDofElemTot,1);
            fi_i  = zeros(this.nDofElemTot,1);
            fe_i  = zeros(this.nDofElemTot,1);

            % Initialize auxiliar variables
            counterK = 0;
            counterF = 0;

            % Update the element displacement vector of each element
            for el = 1:this.nelem
                this.element(el).type.ue = U(this.element(el).type.gle);
            end
            
            % Compute and assemble element data
            for el = 1:this.nelem

                % Get the vector of the element dof  
                gle_i  = this.element(el).type.gle;
                gle_j  = gle_i;
                nGlei  = length(gle_i);
                nGlej  = length(gle_j);
                nGleij = nGlei*nGlej;

                % Get the indices for assemblage
                iDofEl = repmat(gle_i',1,nGlej);
                jDofEl = repmat(gle_j,nGlei,1);
                iDof(counterK+1:counterK+nGleij) = iDofEl(:);
                jDof(counterK+1:counterK+nGleij) = jDofEl(:);
                eDof(counterF+1:counterF+nGlei)  = gle_i';
            
                % Get local matrices and vector
                [K_e,C_e,fi_e,fe_e] = this.element(el).type.elementData();

                % Store in the global vectors
                K_ij(counterK+1:counterK+nGleij) = K_e(:);
                C_ij(counterK+1:counterK+nGleij) = C_e(:);
                fi_i(counterF+1:counterF+nGlei)  = fi_e(:);
                fe_i(counterF+1:counterF+nGlei)  = fe_e(:);
            
                % Update auxiliar variables
                counterK = counterK + nGleij;
                counterF = counterF + nGlei;
                
            end

            % Assemble the matrices and vector
            K  = sparse(iDof,jDof,K_ij);
            C  = sparse(iDof,jDof,C_ij);
            Fi = sparse(eDof,ones(this.nDofElemTot,1),fi_i);
            Fe = sparse(eDof,ones(this.nDofElemTot,1),fe_i);

            % Add contribution of the nodal forces to the external force
            % vector
            Fe = Fe + this.addNodalLoad(Fe);

        end

        function A = applyDirichletBC(this, A)
            % Iterate over all fixed DOFs
            for i = 1:this.ndoffixed
                fixedDOF = this.doffixed(i); 
                
                % Set the row and column of the fixed DOF to zero
                A(fixedDOF, :) = 0; % Clear the row
                A(:, fixedDOF) = 0; % Clear the column
                
                % Set the diagonal entry to 1
                A(fixedDOF, fixedDOF) = 1.0;
            end
        end

        %------------------------------------------------------------------
        % Update the state variables from all integration points
        function updateStateVar(this)
            
            for el = 1:this.nelem
                % Update the element state variables
                this.element(el).type.updateStateVar();
            end

        end

      % -----------------------------------------------------------------
        % Print the nodal displacements
        function printResults(this)
            fprintf('\n******** NODAL RESULTS ********\n');
            this.printResultsHeader();
            for i = 1:this.nnodes
                fprintf("  %4d: \t",i);
                for j = 1:this.ndof_nd
                    fprintf("  %8.4f ",this.U(this.ID(i,j)))
                end
                fprintf("\n");
            end
        end


        % -----------------------------------------------------------------
        % Plot the mesh with the boundary conditions
        function plotMeshWithMatId(this)
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.mesh(true);
        end

        % -----------------------------------------------------------------
        % Plot the deformed mesh
        function plotDeformedMesh(this,amplFactor)

            this.updateResultVertices('Deformed',amplFactor);
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.mesh();

        end

        % -----------------------------------------------------------------
        % Plot the deformed mesh
        function plotField(this,fieldPlot)

            this.updateResultVertexData(fieldPlot)
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.mesh();
%             clim([0.0 0.8]); % set your range here
            c = colorbar;
%             set(c,'Limits',[0.0 0.8])

        end

    end
end