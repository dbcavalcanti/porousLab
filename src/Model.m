%% Model class
%
% Abstract class to create a Finite Element model.
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
        mat                 = [];            % Struct with material properties
        type                = 'ISOQ4';       % Type of element used
        intOrder            = 2;             % Order of the numerical integration quadrature
        nnodes              = 1;             % Number of nodes
        nelem               = 1;             % Number of elements
        nnd_el              = 4;             % Number of nodes per element
        doffree             = [];            % Vector with the free dofs
        doffixed            = [];            % Vector with the fixed dofs
        ndof_nd             = 2;             % Number of dof per node
        ndof                = 1;             % Number of degrees of freedom
        ndoffree            = 0;             % Number of free degrees of freedom
        ndoffixed           = 0;             % Number of fixed degrees of freedom
        Dof                 = [];            % Vector with all regular dofs
        ID                  = [];            % Each line of the ID matrix contains the global numbers for the node DOFs
        U                   = [];            % Global displacement vector
        element             = [];            % Array with the element's objects
        nDofElemTot         = 0.0;           % Aux value used to sparse matrix assemblage
        sqrNDofElemTot      = 0.0;           % Aux value used to sparse matrix assemblage
        matID               = [];            % Vector with the material id of each element
        massLumping         = false;         % Tag for applying a mass lumping process
        lumpStrategy        = 1;             % Id of the mass lumping strategy
        isAxisSymmetric     = false;         % Flag for axissymetric models
        enriched            = false;         % Flag to use embedded formulation
        discontinuitySet    = [];            % Array with the discontinuity objects
    end
    
    %% Constructor method
    methods
        function this = Model()
            disp("*** Initializing model...")
        end
    end

    %% Abstract methods
    methods(Abstract)
       
        % Assemble the nodes' Dirichlet conditions matrix
        SUPP = dirichletConditionMatrix(this);

        % Assemble the nodes' Neumann conditions matrix
        LOAD = neumannConditionMatrix(this);

        % Assemble the nodes' initial condition matrix
        INITCOND = initialConditionMatrix(this)

        % Assemble the nodes' prescribed Dirichlet condition values matrix
        PRESCDISPL = prescribedDirichletMatrix(this);

        % Assemble the degrees of freedom to the elements
        assembleElementDofs(this);

        % Initialize the elements objects
        initializeElements(this);

        % Set the fields available for visualization
        updateResultVertexData(this,type);

        % Configure the header to printed when printing results
        printResultsHeader();
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function createNodeDofIdMtrx(this)
            % Initialize the ID matrix and the number of fixed dof
            this.ID = zeros(this.nnodes,this.ndof_nd);
            this.ndoffixed = 0;

            % Get the Dirichlet conditions matrix
            SUPP = this.dirichletConditionMatrix();
            
            % Assemble the ID matrix
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.ID(i,j) = (i - 1) * this.ndof_nd + j;
                    if (SUPP(i,j) == 1)
                        this.ndoffixed = this.ndoffixed + 1;
                    end
                end
            end

            % Vector with all the dofs
            this.Dof = 1:this.ndof;
            
            % Number of free dof
            this.ndoffree = this.ndof - this.ndoffixed;
            
            % Initialize the counters
            this.doffixed = zeros(this.ndoffixed,1);
            this.doffree  = zeros(this.ndoffree,1);
            
            % Update the ID matrix with the free dof numbered first
            countFree = 1;
            countFixed = 1;
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    if SUPP(i,j) == 1
                        this.doffixed(countFixed) = this.ID(i,j);
                        countFixed = countFixed + 1;
                    else 
                        this.doffree(countFree) = this.ID(i,j);
                        countFree = countFree + 1;
                    end
                end
            end
        end

        %------------------------------------------------------------------
        function preComputations(this)

            disp("*** Pre-processing...");
            
            % Initialize basic variables
            this.initializeBasicVariables();

            % Create nodes DOF ids matrix
            this.createNodeDofIdMtrx();

            % Assemble the regular dofs to each element
            this.assembleElementDofs();

            % Initialize elements
            this.initializeElements();

            % Assemble discontinuity segments to the elements
            this.assembleDiscontinuitySegments();
            
            % Compute auxiliar variables for assemblage of sparse matrices
            this.initializeSparseMtrxAssemblageVariables();

            % Initialize the displacement vector
            this.initializeDisplacementVct();

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
        function initializeDisplacementVct(this)

            % Initialize the displacement vector 
            this.U = zeros(this.ndof,1);

            % Get initial conditions matrix
            INITCOND = this.initialConditionMatrix();

            % Set the initial values
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.U(this.ID(i,j)) = INITCOND(i,j);
                end
            end

            % Get the Dirichlet conditions matrix
            SUPP = this.dirichletConditionMatrix();
            PRESCDISPL = this.prescribedDirichletMatrix();

            % Set the prescribed values
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    if (SUPP(i,j) == 1.0)
                        this.U(this.ID(i,j)) = PRESCDISPL(i,j);
                    end
                end
            end

            % Save initial dofs to the elements
            for el = 1 : this.nelem
                this.element(el).type.ue = this.U(this.element(el).type.gle);
            end

        end

        %------------------------------------------------------------------
        function assembleDiscontinuitySegments(this)
            if (this.enriched == false)
                return
            elseif isempty(this.discontinuitySet)
                return
            else
                nDiscontinuities = size(this.discontinuitySet,1);
                for i = 1:nDiscontinuities
                    % Loop through the segments of this discontinuity
                    nDiscontinuitySeg = size(this.discontinuitySet(i).elemID,1);
                    for j = 1:nDiscontinuitySeg
                        el = this.discontinuitySet(i).elemID(j);
                        dseg = this.discontinuitySet(i).segment(j);
                        this.element(el).type.addDiscontinuitySegment(dseg);
                    end
                end
            end
        end

        %------------------------------------------------------------------
        % Add contribution of nodal loads to reference load vector.
        function Fref = addNodalLoad(this,Fref)
            LOAD = this.neumannConditionMatrix();
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    Fref(this.ID(i,j)) = Fref(this.ID(i,j)) + LOAD(i,j);
                end
            end
        end

        %------------------------------------------------------------------
        function Lce = getElementsCharacteristicLength(this)
            Lce=zeros(size(this.ELEM,1),1);
            for el = 1:size(this.ELEM,1)
                % Characteristic lenght (quadrilateral elements)
                Lce(el) = this.getElementCharacteristicLength(el);
            end
        end

        %------------------------------------------------------------------
        function Lce = getElementCharacteristicLength(this,el)

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
            Lce = sqrt(Ae);
            if strcmp(this.type,'CST')||strcmp(this.type,'LST')
                Lce = Lce * sqrt(2.0);
            end
        end
        
        %------------------------------------------------------------------
        % Compute the mean characteristic length of the elements associated with
        % each node
        function Lc = getNodeCharacteristicLength(this)
            Lce = this.getElementsCharacteristicLength();
            Lc = zeros(size(this.NODE,1),1);
            for i = 1:size(this.NODE,1)
                % Get the elements associated with this node
                idElem = any(this.ELEM == i, 2);
                % Compute the mean characteristic lenght of these nodes
                Lc(i) = mean(Lce(idElem));
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
        function [K, C, Fi, Fe, dfidx] = globalMatrices(this,U)   

            % Indices for the assembling the matrices and vector
            iDof = zeros(this.sqrNDofElemTot,1);
            jDof = zeros(this.sqrNDofElemTot,1);
            eDof = zeros(this.nDofElemTot,1);

            % Initialize the components of the matrices and vector
            K_ij      = zeros(this.sqrNDofElemTot,1);
            C_ij      = zeros(this.sqrNDofElemTot,1);
            dfidx_ij  = zeros(this.sqrNDofElemTot,1);
            fi_i      = zeros(this.nDofElemTot,1);
            fe_i      = zeros(this.nDofElemTot,1);

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
                [K_e,C_e,fi_e,fe_e,dfidx_e] = this.element(el).type.elementData();

                % Store in the global vectors
                K_ij(counterK+1:counterK+nGleij)     = K_e(:);
                C_ij(counterK+1:counterK+nGleij)     = C_e(:);
                dfidx_ij(counterK+1:counterK+nGleij) = dfidx_e(:);
                fi_i(counterF+1:counterF+nGlei)      = fi_e(:);
                fe_i(counterF+1:counterF+nGlei)      = fe_e(:);
            
                % Update auxiliar variables
                counterK = counterK + nGleij;
                counterF = counterF + nGlei;
                
            end

            % Assemble the matrices and vector
            K      = sparse(iDof,jDof,K_ij);
            C      = sparse(iDof,jDof,C_ij);
            dfidx  = sparse(iDof,jDof,dfidx_ij);
            Fi     = sparse(eDof,ones(this.nDofElemTot,1),fi_i);
            Fe     = sparse(eDof,ones(this.nDofElemTot,1),fe_i);

            % Add contribution of the nodal forces to the external force
            % vector
            Fe = this.addNodalLoad(Fe);

        end

        %------------------------------------------------------------------
        % Global system matrices
        function [A,b] = getLinearSystem(this,U,UOld,nonlinearScheme,dt)   

            % Indices for the assembling the matrices and vector
            iDof = zeros(this.sqrNDofElemTot,1);
            jDof = zeros(this.sqrNDofElemTot,1);
            eDof = zeros(this.nDofElemTot,1);

            % Initialize the components of the matrices and vector
            A_ij  = zeros(this.sqrNDofElemTot,1);
            b_i  = zeros(this.nDofElemTot,1);

            % Initialize auxiliar variables
            counterK = 0;
            counterF = 0;

            % Update the element displacement vector of each element
            for el = 1:this.nelem
                this.element(el).type.DTime = dt;
                this.element(el).type.ue    = U(this.element(el).type.gle);
                this.element(el).type.ueOld = UOld(this.element(el).type.gle);
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
                [A_e,b_e] = this.element(el).type.elementLinearSystem(nonlinearScheme);

                % Store in the global vectors
                A_ij(counterK+1:counterK+nGleij) = A_e(:);
                b_i(counterF+1:counterF+nGlei)   = b_e(:);
            
                % Update auxiliar variables
                counterK = counterK + nGleij;
                counterF = counterF + nGlei;
                
            end

            % Assemble the matrices and vector
            A = sparse(iDof,jDof,A_ij);
            b = sparse(eDof,ones(this.nDofElemTot,1),b_i);
            b = full(b);

            % Add contribution of the nodal forces
            Fe = sparse(this.ndof,1);
            Fe = this.addNodalLoad(Fe);
            b = nonlinearScheme.addNodalForces(b,Fe);

        end

        %------------------------------------------------------------------
        function [Aff,bf] = applyDirichletBC(this, A, b, X, nlscheme)
            Aff = A(this.doffree,this.doffree);
            bf  = nlscheme.applyBCtoRHS(A, b, X, this.doffree,this.doffixed);
        end

        %------------------------------------------------------------------
        % Update the state variables from all integration points
        function updateStateVar(this)
            
            for el = 1:this.nelem
                % Update the element state variables
                this.element(el).type.updateStateVar();
            end

        end

        %------------------------------------------------------------------
        function resequenceNodes(this)
            % Get auxiliar variables
            nNode   = size(this.NODE,1);
            nElem   = size(this.ELEM,1);
            nNdElem = size(this.ELEM,2);
            % Size of the connectivity matrix
            nn = nElem*(nNdElem * nNdElem);
            % Get connectivity matrix
            i=zeros(nn,1); j=zeros(nn,1); s=zeros(nn,1); index=0;
            for el = 1:nElem
              eNode=this.ELEM(el,:);
              ElemSet=index+1:index+nNdElem^2;
              i(ElemSet) = kron(eNode,ones(nNdElem,1))';
              j(ElemSet) = kron(eNode,ones(1,nNdElem))';
              s(ElemSet) = 1;
              index = index + nNdElem^2;
            end
            K = sparse(i,j,s,nNode, nNode);
            % Apply a Symmetric reverse Cuthill-McKee permutation
            p = symrcm(K);
            cNode(p(1:nNode))=1:nNode;
            % Rebuild the nodes and elements matrices
            this.rebuildConnectivity(cNode);
        end

        %------------------------------------------------------------------
        function rebuildConnectivity(this,cNode)
            [~,ix,jx] = unique(cNode);
            this.NODE = this.NODE(ix,:);
            for el=1:size(this.ELEM,1)
                for i = 1:size(this.ELEM,2)
                    this.ELEM(el,i) = jx(this.ELEM(el,i));
                end
            end
        end

        %------------------------------------------------------------------
        function addPreExistingDiscontinuities(this,dSet)
            this.discontinuitySet = dSet;
            this.useEnrichedFormulation(true);
            this.initializeDiscontinuitySegments();
        end

        %------------------------------------------------------------------
        function initializeDiscontinuitySegments(this)
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                nDiscontinuitySeg = this.discontinuitySet(i).getNumberOfDiscontinuitySegments();
                for j = 1:nDiscontinuitySeg
                    this.discontinuitySet(i).segment(j).t = this.t;
                end
            end
        end

        %------------------------------------------------------------------
        function n = getNumberOfDiscontinuities(this)
            n = size(this.discontinuitySet,1);
        end

        %------------------------------------------------------------------
        function useEnrichedFormulation(this,flag)
            this.enriched = flag;
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
        % Plot given field over the mesh
        function plotField(this,fieldPlot,range)
            if nargin < 3, range = []; end

            this.updateResultVertexData(fieldPlot)
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.mesh();
            if isempty(range)
                colorbar;
            else
                clim(range);
                c = colorbar;
                set(c,'Limits',range)
            end

        end

    end
end