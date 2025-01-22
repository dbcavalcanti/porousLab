%% Model_H2 class
%
% Two-phase flow finite element model
%
%% Author
% Danilo Cavalcanti
%
%
%% Class definition
classdef Model_H2 < Model    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        SUPP_p              = [];            % Matrix with support conditions
        LOAD_p              = [];            % Matrix with load conditions
        PRESCDISPL_p        = [];            % With with prescribed displacements
        INITCOND_p          = [];            % Initial displacement 
        SUPP_pg             = [];            % Matrix with support conditions
        LOAD_pg             = [];            % Matrix with load conditions
        PRESCDISPL_pg       = [];            % With with prescribed displacements
        INITCOND_pg         = [];            % Initial displacement 
        pDof                = [];            % Vector with the pressure dofs
        pgDof               = [];            % Vector with the pressure dofs
        pFreeDof            = [];            % Vector with the pressure dofs
        pgFreeDof           = [];            % Vector with the pressure dofs
        GLP                 = [];            % Matrix with the regular dof of each element
        GLPg                = [];            % Matrix with the regular dof of each element
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Model_H2()
            this = this@Model();
            this.ndof_nd = 2;
            this.physics = 'H2';
        end
    end
    
    %% Public methods
    methods

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
        function assembleElementDofs(this)

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
                % Create the material for the element
                emat =struct( ...
                        'porousMedia',this.mat.porousMedia(this.matID(el)), ...
                        'liquidFluid',this.mat.liquidFluid,...
                        'gasFluid',this.mat.gasFluid);
                if strcmp(this.physics,'H2')
                    elements(el) = RegularElement_H2(...
                            this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.t, emat, this.intOrder,this.GLP(el,:),this.GLPg(el,:), ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                elseif strcmp(this.physics,'H2_PcPg')
                    elements(el) = RegularElement_H2_PcPg(...
                            this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.t, emat, this.intOrder,this.GLP(el,:),this.GLPg(el,:), ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
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
        % Add contribution of nodal loads to reference load vector.
        function Fref = addNodalLoad(this,Fref)
            for i = 1:this.nnodes
                Fref(this.ID(i,1)) = Fref(this.ID(i,1)) + this.LOAD_p(i,1);
                Fref(this.ID(i,2)) = Fref(this.ID(i,2)) + this.LOAD_pg(i,1);
            end
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
        function plotGasPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotGasPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        % -----------------------------------------------------------------
        % Plot the mesh with the boundary conditions
        function plotCapillaryPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotCapillaryPressureAlongSegment(Xi, Xf, npts,axisPlot);
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
                    elseif strcmp(type,'LiquidPressure')
                        p = this.element(el).type.pressureField(X);
                        vertexData(i) = p;
                    elseif strcmp(type,'CapillaryPressure')
                        p = this.element(el).type.capillaryPressureField(X);
                        vertexData(i) = p;
                    elseif strcmp(type,'GasPressure')
                        p = this.element(el).type.gasPressureField(X);
                        vertexData(i) = p;
                    elseif strcmp(type,'LiquidSaturation')
                        Sl = this.element(el).type.liquidSaturationField(X);
                        vertexData(i) = Sl;
                    elseif strcmp(type,'GasSaturation')
                        Sg = this.element(el).type.gasSaturationField(X);
                        vertexData(i) = Sg;
                    end
                end
                this.element(el).type.result.setVertexData(vertexData);
            end
        end

    end

    %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           Pl        Pg\n');
        end

    end
end