%% Model_H2 class
%
% Two-phase flow finite element model.
%
% Each node has two degrees of freedom (dof). The liquid phase pressure (p)
% and the gas phase pressure (pg).
%
%% Author
% Danilo Cavalcanti
%
%
%% Class definition
classdef Model_H2 < Model    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        %% Degrees of freedom vectors
        % Vector with all the dofs of each type
        pDof                = [];
        pgDof               = [];
        % Vector with all the FREE dofs of each type
        pFreeDof            = [];
        pgFreeDof           = [];
        % Matrix with the dofs of each type of each element
        GLP                 = [];
        GLPg                = [];
        %% Matrix indicating the Dirichlet BCs 
        SUPP_p              = [];
        SUPP_pg             = [];
        %% Matrix with the prescribed BC values 
        PRESCDISPL_p        = [];
        PRESCDISPL_pg       = [];
        %% Matrix with the Neumann BCs
        LOAD_p              = [];
        LOAD_pg             = [];
        %% Matrix with the initial conditions
        INITCOND_p          = []; 
        INITCOND_pg         = [];
    end
    %% Constructor method
    methods
        function this = Model_H2()
            this = this@Model();
            this.ndof_nd = 2;       % Number of dofs per node
            this.physics = 'H2';    % Tag with the physics name
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function SUPP = dirichletConditionMatrix(this)
            SUPP = [this.SUPP_p , this.SUPP_pg];
        end

        %------------------------------------------------------------------
        function LOAD = neumannConditionMatrix(this)
            LOAD = [this.LOAD_p , this.LOAD_pg];  
        end

        %------------------------------------------------------------------
        function INITCOND = initialConditionMatrix(this)
            INITCOND = [this.INITCOND_p , this.INITCOND_pg];  
        end

        %------------------------------------------------------------------
        function PRESCDISPL = prescribedDirichletMatrix(this)
            PRESCDISPL = [this.PRESCDISPL_p , this.PRESCDISPL_pg];  
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
                elements(el) = RegularElement_H2(...
                            this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                            this.t, emat, this.intOrder,this.GLP(el,:),this.GLPg(el,:), ...
                            this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
        end

        % -----------------------------------------------------------------
        function plotPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        % -----------------------------------------------------------------
        function plotGasPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotGasPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        % -----------------------------------------------------------------
        function plotCapillaryPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotCapillaryPressureAlongSegment(Xi, Xf, npts,axisPlot);
        end

        %------------------------------------------------------------------
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