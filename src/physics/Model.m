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
        DIRICHLET_TAG       = [];            % Matrix with tags indicating the Dirichlet BC
        DIRICHLET_VAL       = [];            % Matrix with values of the Dirichlet BC
        LOAD                = [];            % Matrix with the nodal Neumann BC
        INIT                = [];            % Matrix with the initial values of each dof
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
        nDofElemTot         = 0;             % Aux value used to sparse matrix assemblage
        sqrNDofElemTot      = 0;             % Aux value used to sparse matrix assemblage
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
        function setMesh(this,node,elem)
            % Set the mesh nodes coordinates and connectivity
            this.NODE = node;
            this.ELEM = elem;

            % Initialize basic variables
            this.initializeBasicVariables();

        end

        %------------------------------------------------------------------
        function initializeBasicVariables(this)
            this.nnodes        = size(this.NODE,1);
            this.nelem         = size(this.ELEM,1);     
            this.nnd_el        = size(this.ELEM,2);            
            this.ndof          = this.ndof_nd * this.nnodes; 
            this.DIRICHLET_TAG = zeros(this.nnodes,this.ndof_nd);
            this.DIRICHLET_VAL = zeros(this.nnodes,this.ndof_nd);
            this.LOAD          = zeros(this.nnodes,this.ndof_nd);
            this.INIT          = zeros(this.nnodes,this.ndof_nd);
            if (this.nnd_el == 3)
                this.type = 'CST';
            elseif (this.nnd_el == 4)
                this.type = 'ISOQ4';
            elseif (this.nnd_el == 6)
                this.type = 'LST';
            elseif (this.nnd_el == 8)
                this.type = 'ISOQ8';
            end
        end

        %------------------------------------------------------------------
        function checkMaterialId(this)
            if isempty(this.matID)
                this.matID  = ones(this.nelem,1);
            end
        end

        %------------------------------------------------------------------
        function createNodeDofIdMtrx(this)
            % Initialize the ID matrix and the number of fixed dof
            this.ID = zeros(this.nnodes,this.ndof_nd);
            this.ndoffixed = 0;
            
            % Assemble the ID matrix
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.ID(i,j) = (i - 1) * this.ndof_nd + j;
                    if (this.DIRICHLET_TAG(i,j) == 1)
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
                    if this.DIRICHLET_TAG(i,j) == 1
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
        function setDirichletBCAtNode(this, nodeId, dofId, value)
            this.DIRICHLET_TAG(nodeId,dofId) = repmat( ~isnan(value), sum(nodeId), 1);
            this.DIRICHLET_VAL(nodeId,dofId) = repmat(value, sum(nodeId), 1);
        end

        %------------------------------------------------------------------
        function setDirichletBCAtPoint(this, X, dofId, value)
            nodeId = this.closestNodeToPoint(X);
            this.DIRICHLET_TAG(nodeId,dofId) = repmat( ~isnan(value), sum(nodeId), 1);
            this.DIRICHLET_VAL(nodeId,dofId) = repmat(value, sum(nodeId), 1);
        end

        %------------------------------------------------------------------
        function setDirichletBCAtBorder(this, border, dofId, value)
            % Get the nodes at the given border
            if strcmp(border,'left')
                nodeId = abs(this.NODE(:,1)-min(this.NODE(:,1)))<1.0e-12;
            elseif strcmp(border,'right')
                nodeId = abs(this.NODE(:,1)-max(this.NODE(:,1)))<1.0e-12;
            elseif strcmp(border,'top')
                nodeId = abs(this.NODE(:,2)-max(this.NODE(:,2)))<1.0e-12;
            elseif strcmp(border,'bottom')
                nodeId = abs(this.NODE(:,2)-min(this.NODE(:,2)))<1.0e-12;
            else
                disp('Warning: non-supported border.');
                disp('Available borders tag: ''left'',''right'', ''top'',''bottom''');
            end
            this.DIRICHLET_TAG(nodeId,dofId) = repmat(~isnan(value), sum(nodeId), 1);
            this.DIRICHLET_VAL(nodeId,dofId) = repmat(value, sum(nodeId), 1);
        end
        
        %------------------------------------------------------------------
        function nd = closestNodeToPoint(this,X)
            if size(X,1) == 2, X = X'; end
            d = vecnorm((this.NODE - X)');
            [~,id] = sort(d);
            nd = id(1);
        end

        %------------------------------------------------------------------
        function preComputations(this)

            disp("*** Pre-processing...");
            
            % Check and initialize the material ID vector
            this.checkMaterialId();

            % Create nodes DOF ids matrix
            this.createNodeDofIdMtrx();

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
        function dof = getElementDofs(this,el,dofId)
            dof = reshape(this.ID(this.ELEM(el,:),dofId)', 1, this.nnd_el*length(dofId));
        end

        %------------------------------------------------------------------
        function initializeDisplacementVct(this)

            % Initialize the displacement vector 
            this.U = zeros(this.ndof,1);

            % Set the initial values
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.U(this.ID(i,j)) = this.INIT(i,j);
                end
            end

            % Set the prescribed values
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    if (this.DIRICHLET_TAG(i,j) == 1.0)
                        this.U(this.ID(i,j)) = this.DIRICHLET_VAL(i,j);
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
                nDiscontinuities = length(this.discontinuitySet);
                for i = 1:nDiscontinuities
                    % Loop through the segments of this discontinuity
                    k = 1;
                    for j = 1:size(this.discontinuitySet(i).Xlin,1)-1
                        el = this.discontinuitySet(i).elemID(j);
                        if (el > 0)
                            dseg = this.discontinuitySet(i).segment(k);
                            this.element(el).type.addDiscontinuitySegment(dseg);
                            k = k + 1;
                        end
                    end
                end
            end
        end

        %------------------------------------------------------------------
        % Add contribution of nodal loads to reference load vector.
        function Fref = addNodalLoad(this,Fref)
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    Fref(this.ID(i,j)) = Fref(this.ID(i,j)) + this.LOAD(i,j);
                end
            end
        end

        %------------------------------------------------------------------
        function Lce = getElementsCharacteristicLength(this)
            Lce=zeros(this.nelem,1);
            for el = 1:this.nelem
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
        % Compute the mean characteristic length of the elements associated 
        % with each node
        function Lc = getNodeCharacteristicLength(this)
            Lce = this.getElementsCharacteristicLength();
            Lc = zeros(this.nnodes,1);
            for i = 1:this.nnodes
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
        function addPreExistingDiscontinuities(this,dSet,additionalData)
            if nargin > 2
                this.addDiscontinuityData(additionalData);
            end
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
            n = length(this.discontinuitySet);
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
                    fprintf("  %8.4e ",this.U(this.ID(i,j)))
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