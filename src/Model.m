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
        physics             = 'hydroTwoPhase';
        NODE                = [];            % Nodes of the fem mesh
        ELEM                = [];            % Nodes connectivity
        t                   = 1.0;           % Thickness
        mat                 = [];            % Vector with material properties
        type                = 'ISOQ4';       % Typf of element used
        SUPP_p              = [];            % Matrix with support conditions
        LOAD_p              = [];            % Matrix with load conditions
        PRESCDISPL_p        = [];            % With with prescribed displacements
        INITCOND_p          = [];            % Initial displacement 
        SUPP_pg             = [];            % Matrix with support conditions
        LOAD_pg             = [];            % Matrix with load conditions
        PRESCDISPL_pg       = [];            % With with prescribed displacements
        INITCOND_pg         = [];            % Initial displacement 
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
        pDof                = [];            % Vector with the pressure regular dofs
        pgDof               = [];            % Vector with the pressure regular dofs
        pFreeDof            = [];            % Vector with the pressure regular dofs
        pgFreeDof           = [];            % Vector with the pressure regular dofs
        Dof                 = [];            % Vector with all regular dofs
        ID                  = [];            % Each line of the ID matrix contains the global numbers for the node DOFs (DX, DY)
        GLP                 = [];            % Matrix with the regular dof of each element
        GLPg                 = [];            % Matrix with the regular dof of each element
        F                   = [];            % Global force vector
        U                   = [];            % Global displacement vector
        element             = [];            % Array with the element's objects
        nDofElemTot         = 0.0;
        sqrNDofElemTot      = 0.0;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Model(physics,NODE, ELEM, t, mat, type, ...
                SUPP_p, LOAD_p, PRESCDISPL_p, INITCOND_p,...
                SUPP_pg, LOAD_pg, PRESCDISPL_pg, INITCOND_pg, intOrder)
            if (nargin > 0)
                this.physics            = physics;
                this.NODE               = NODE;
                this.ELEM               = ELEM;
                this.type               = type;
                this.t                  = t;
                this.mat                = mat;
                this.SUPP_p             = SUPP_p;
                this.SUPP_pg            = SUPP_pg;
                this.LOAD_p             = LOAD_p;
                this.LOAD_pg            = LOAD_pg;
                this.PRESCDISPL_p       = PRESCDISPL_p;
                this.PRESCDISPL_pg      = PRESCDISPL_pg;
                this.INITCOND_p         = INITCOND_p;
                this.INITCOND_pg        = INITCOND_pg;
                this.intOrder           = intOrder;
            end
        end
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
            this.assembleElementsRegularDofs();

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
            this.nnodes     = size(this.NODE,1);
            this.nelem      = size(this.ELEM,1);     
            this.nnd_el     = size(this.ELEM,2);    
            this.ndof_nd    = 2;               
            this.ndof       = this.ndof_nd * this.nnodes;         
        end

        %------------------------------------------------------------------
        %  Assemble nodes DOF ids matrix
        %   Each line of the ID matrix contains the global numbers for 
        %   the node DOFs (P). Free DOFs are numbered first.
        function createNodeDofIdMtrx(this)
            % Initialize the ID matrix and the number of fixed dof
            this.ID = zeros(this.nnodes,this.ndof_nd);
            this.ndoffixed = 0;

            SUPP = [this.SUPP_p , this.SUPP_pg];
            
            % Assemble the ID matrix
            for i = 1:this.nnodes
                for j = 1:this.ndof_nd
                    this.ID(i,j) = (i - 1) * this.ndof_nd + j;
                    if (SUPP(i,j) == 1)
                        this.ndoffixed = this.ndoffixed + 1;
                    end
                end
            end
            
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
        function assembleElementsRegularDofs(this)

            this.GLP = zeros(this.nelem, this.nnd_el);
            for el = 1:this.nelem
                this.GLP(el,:) = reshape(this.ID(this.ELEM(el,:),1)',1,...
                    this.nnd_el);
            end

            this.GLPg = zeros(this.nelem, this.nnd_el);
            for el = 1:this.nelem
                this.GLPg(el,:) = reshape(this.ID(this.ELEM(el,:),2)',1,...
                    this.nnd_el);
            end

            % Vector with all regular dofs
            this.pDof = unique(this.GLP);
            this.pgDof = unique(this.GLPg);
            this.Dof  = [this.pDof(:); this.pgDof(:)];

            % Vector will free regular dofs
            this.pFreeDof  = intersect(this.pDof,this.doffree);
            this.pgFreeDof = intersect(this.pgDof,this.doffree);
        end

        %------------------------------------------------------------------
        function initializeElements(this)
            % Initialize the vector with the Element's objects
            elements(this.nelem,1) = Element(); 

            % Assemble the properties to the elements' objects
            for el = 1 : this.nelem
                elements(el) = RegularElement(...
                        this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                        this.t, this.mat, this.intOrder,this.GLP(el,:),this.GLPg(el,:));

                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
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
        function initializeDisplacementVct(this)

            % Initialize the displacement vector 
            this.U = zeros(this.ndof,1);

            % Set the initial values
            for i = 1:this.nnodes
                this.U(this.ID(i,1)) = this.INITCOND_p(i,1);
                this.U(this.ID(i,2)) = this.INITCOND_pg(i,1);
            end

            % Set the prescribed values
            for i = 1:this.nnodes
                if (this.SUPP_p(i,1) == 1.0)
                    this.U(this.ID(i,1)) = this.PRESCDISPL_p(i,1);
                end
                if (this.SUPP_pg(i,1) == 1.0)
                    this.U(this.ID(i,2)) = this.PRESCDISPL_pg(i,1);
                end
            end

            % Save initial dofs to the elements
            for el = 1 : this.nelem
                this.element(el).type.ue = this.U(this.element(el).type.gle);
            end

        end
        
        %------------------------------------------------------------------
        function initializeExtForceVct(this)
            this.F = zeros(this.ndof,1);
        end
        
        %------------------------------------------------------------------
        % Add contribution of nodal loads to reference load vector.
        function Fref = addNodalLoad(this,Fref)
            for i = 1:this.nnodes
                Fref(this.ID(i,1)) = Fref(this.ID(i,1)) + this.LOAD_p(i,1);
                Fref(this.ID(i,2)) = Fref(this.ID(i,2)) + this.LOAD_pg(i,1);
            end
        end

        %------------------------------------------------------------------
        % Add contribution of nodal loads to reference load vector.
        function Fref = addBodyForceLoad(this,Fref)
            for el = 1:this.nelem

                % Get the vector of the element dof
                gle = this.element(el).type.glu;
                
                % Element body force vector
                Fb = this.element(el).type.bodyForce();
            
                % Assemble
                Fref(gle) = Fref(gle) + Fb;
                
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
        function [K, C, Fi] = globalMatrices(this,U)   

            % Indices for the assembling the matrices and vector
            iDof = zeros(this.sqrNDofElemTot,1);
            jDof = zeros(this.sqrNDofElemTot,1);
            eDof = zeros(this.nDofElemTot,1);

            % Initialize the components of the matrices and vector
            Kij  = zeros(this.sqrNDofElemTot,1);
            Cij  = zeros(this.sqrNDofElemTot,1);
            fi   = zeros(this.nDofElemTot,1);

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
                [Ke,Ce,fe] = this.element(el).type.elementData();

                % Store in the global vectors
                Kij(counterK+1:counterK+nGleij) = Ke(:);
                Cij(counterK+1:counterK+nGleij) = Ce(:);
                fi(counterF+1:counterF+nGlei)   = fe(:);
            
                % Update auxiliar variables
                counterK = counterK + nGleij;
                counterF = counterF + nGlei;
                
            end

            % Assemble the matrices and vector
            K  = sparse(iDof,jDof,Kij);
            C  = sparse(iDof,jDof,Cij);
            Fi = sparse(eDof,ones(this.nDofElemTot,1),fi);

            % Apply the boundary conditions
            K = this.applyDirichletBC(K);
            C = this.applyDirichletBC(C);

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
        % Global fluid-flow and compressibility matrices and discharge vector
        function Fint = globalInternalForceVector(this,dX,X)   
            
            % Initialize global matrices
            Fint = zeros(this.ndof, 1);
            
            for el = 1:this.nelem

                % Get the vector of the element dof
                gle = this.element(el).type.gle;
            
                % Get local stiffness matrix
                fie = this.element(el).type.elementInternalForceVector(dX(gle),X(gle));
            
                % Assemble
                Fint(gle) = Fint(gle) + fie;
                
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
            fprintf('\n  Node           Pl        Pg\n');
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
        function plotDisplacementAlongSegment(this, dir, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotDisplacementAlongSegment(dir, Xi, Xf, npts,axisPlot);
        end

        % -----------------------------------------------------------------
        % Plot the mesh with the boundary conditions
        function plotPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        % -----------------------------------------------------------------
        % Plot the mesh with the boundary conditions
        function plotMeshWithBC(this)
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.mesh();
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
            % clim([0.0 10]); % set your range here
            c = colorbar;
            % set(c,'Limits',[0.0 10])

        end

        % -----------------------------------------------------------------
        % Returns a cell with the elements that share the nodes of the
        % fracture mesh
        function idCell = findElementsShareNode(this)

            % Initialize variables
            nfracNodes = size(this.NODE_D,1);
            idMtrx     = zeros(nfracNodes,this.nelem);
            
            % Fill the idMtrx with ones if the node is in the element domain/boundary
            for i = 1:nfracNodes
                for el = 1:this.nelem
                    idMtrx(i,el) = findElementInMesh(this.NODE, this.ELEM(el,:), this.NODE_D(i,:));
                end
            end
            
            % Convert into a cell
            idCell = cell(nfracNodes,1);
            for i = 1:nfracNodes
                idCell{i} = find(idMtrx(i,:));
            end
            
        end

        %------------------------------------------------------------------
        % Update the result nodes coordinates of each element
        function updateResultVertices(this,configuration,factor)
            
            for el = 1:this.nelem
                
                % Initialize the vertices array
                vertices = this.element(el).type.result.vertices0;

                % Get the updated vertices:
                if strcmp(configuration,'Deformed')

                    % Update the nodal displacement vector associated to the
                    % element. This displacement can contain the enhancement
                    % degrees of freedom.
                    this.element(el).type.ue = this.U(this.element(el).type.gle); 

                    % Update the vertices based on the displacement vector
                    % associated to the element
                    for i = 1:length(this.element(el).type.result.faces)
                        X = vertices(i,:);
                        u = this.element(el).type.displacementField(X);
                        vertices(i,:) = X + factor*u';
                    end
                end
                this.element(el).type.result.setVertices(vertices);

            end
            
        end

        %------------------------------------------------------------------
        % Update the result nodes data of each element
        function updateResultVertexData(this,type)
            for el = 1:this.nelem
                % Update the nodal displacement vector associated to the
                % element. This displacement can contain the enhancement
                % degrees of freedom.
                this.element(el).type.ue = this.U(this.element(el).type.gle); 
                vertexData = zeros(length(this.element(el).type.result.faces),1);
                for i = 1:length(this.element(el).type.result.faces)
                    X = this.element(el).type.result.vertices(i,:);
                    if strcmp(type,'Model')
                        vertexData(i) = 0.0;
                    elseif strcmp(type,'Pressure')
                        p = this.element(el).type.pressureField(X);
                        vertexData(i) = p;
                    elseif strcmp(type,'Ux')
                        u = this.element(el).type.displacementField(X);
                        vertexData(i) = u(1);
                    elseif strcmp(type,'Uy')
                        u = this.element(el).type.displacementField(X);
                        vertexData(i) = u(2);
                    elseif strcmp(type,'Sx')
                        s = this.element(el).type.stressField(X);
                        vertexData(i) = s(1);
                    elseif strcmp(type,'Sy')
                        s = this.element(el).type.stressField(X);
                        vertexData(i) = s(2);
                    end
                end
                this.element(el).type.result.setVertexData(vertexData);
            end
        end

    end
end