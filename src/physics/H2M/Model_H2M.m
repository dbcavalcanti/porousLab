%% Model class
%
% This class defines a finite element model that has a strong discontinuity
%
%% Author
% Danilo Cavalcanti
%
%
%% Class definition
classdef Model_H2M < Model    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        %% Degrees of freedom vectors
        % Vector with all the dofs of each type
        uDof                = [];
        pDof                = [];
        pgDof               = [];
        % Vector with all the FREE dofs of each type
        uFreeDof            = [];
        pFreeDof            = [];
        pgFreeDof           = [];
        % Matrix with the dofs of each type of each element
        GLU                 = [];
        GLP                 = [];
        GLPg                = [];
        %% Matrix indicating the Dirichlet BCs 
        SUPP_u              = [];
        SUPP_p              = [];
        SUPP_pg             = [];
        %% Matrix with the prescribed BC values 
        PRESCDISPL_u        = [];
        PRESCDISPL_p        = [];
        PRESCDISPL_pg       = [];
        %% Matrix with the Neumann BCs
        LOAD_u              = [];
        LOAD_p              = [];
        LOAD_pg             = [];
        %% Matrix with the initial conditions
        INITCOND_p          = []; 
        INITCOND_pg         = [];
        %% Additional data
        isPlaneStress       = false;
    end
    
    %% Constructor method
    methods
        function this = Model_H2M()
            this = this@Model();
            this.ndof_nd = 4;        % Number of dofs per node
            this.physics = 'H2M';    % Tag with the physics name
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function SUPP = dirichletConditionMatrix(this)
            SUPP = [this.SUPP_u ,this.SUPP_p , this.SUPP_pg];  
        end

        %------------------------------------------------------------------
        function LOAD = neumannConditionMatrix(this)
            LOAD = [this.LOAD_u ,this.LOAD_p , this.LOAD_pg];  
        end

        %------------------------------------------------------------------
        function INITCOND = initialConditionMatrix(this)
            INITCOND_u = zeros(this.nnodes,2);
            INITCOND = [INITCOND_u ,this.INITCOND_p , this.INITCOND_pg];  
        end

        %------------------------------------------------------------------
        function PRESCDISPL = prescribedDirichletMatrix(this)
            PRESCDISPL = [this.PRESCDISPL_u ,this.PRESCDISPL_p , this.PRESCDISPL_pg];  
        end

        %------------------------------------------------------------------
        function assembleElementDofs(this)

            this.GLU = zeros(this.nelem, this.nnd_el*2);
            for el = 1:this.nelem
                this.GLU(el,:) = reshape(this.ID(this.ELEM(el,:),1:2)',1,...
                    this.nnd_el*2);
            end
            this.GLP = zeros(this.nelem, this.nnd_el);
            for el = 1:this.nelem
                this.GLP(el,:) = reshape(this.ID(this.ELEM(el,:),3)',1,...
                    this.nnd_el);
            end
            this.GLPg = zeros(this.nelem, this.nnd_el);
            for el = 1:this.nelem
                this.GLPg(el,:) = reshape(this.ID(this.ELEM(el,:),4)',1,...
                    this.nnd_el);
            end

            % Vector with all regular dofs
            this.uDof = unique(this.GLU);
            this.pDof = unique(this.GLP);
            this.pgDof = unique(this.GLPg);
            this.Dof  = [this.uDof(:); this.pDof(:); this.pgDof(:)];

            % Vector will free regular dofs
            this.uFreeDof  = intersect(this.uDof,this.doffree);
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
                elements(el) = RegularElement_H2M(...
                            this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.t, emat, this.intOrder,this.GLU(el,:),this.GLP(el,:),this.GLPg(el,:), ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                            this.isPlaneStress);
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
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

        % -----------------------------------------------------------------
        % Plot the deformed mesh
        function plotDeformedMesh(this,amplFactor)

            this.updateResultVertices('Deformed',amplFactor);
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.mesh();

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
                        vertexData(i) = this.matID(el);
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
        %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        % Print header of the results
        function printResultsHeader()
            fprintf('\n  Node           ux        uy        Pl        Pg\n');
        end

    end
end