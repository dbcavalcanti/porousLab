%% Model_H2 class
%
% Single-phase fluid flow finite element model.
%
% Each node has one degree of freedom (dof) which is the pore-pressure (p)
%
%% Author
% Danilo Cavalcanti
%
%
%% Class definition
classdef Model_H < Model    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        %% Degrees of freedom vectors
        % Matrix with the dofs of each type of each element
        GLP                 = [];
        %% Matrix indicating the Dirichlet BCs 
        SUPP_p              = [];
        %% Matrix with the prescribed BC values 
        PRESCDISPL_p        = [];
        %% Matrix with the Neumann BCs
        LOAD_p              = [];
        %% Matrix with the initial conditions
        INITCOND_p          = []; 
    end
    %% Constructor method
    methods
        function this = Model_H()
            this = this@Model();
            this.ndof_nd = 1;       % Number of dofs per node
            this.physics = 'H';     % Tag with the physics name
            disp("*** Physics: Single-phase fluid flow");
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function SUPP = dirichletConditionMatrix(this)
            SUPP = this.SUPP_p;  
        end

        %------------------------------------------------------------------
        function LOAD = neumannConditionMatrix(this)
            LOAD = this.LOAD_p;  
        end

        %------------------------------------------------------------------
        function INITCOND = initialConditionMatrix(this)
            INITCOND = this.INITCOND_p;  
        end

        %------------------------------------------------------------------
        function PRESCDISPL = prescribedDirichletMatrix(this)
            PRESCDISPL = this.PRESCDISPL_p;  
        end
        %------------------------------------------------------------------
        function assembleElementDofs(this)

            this.GLP = zeros(this.nelem, this.nnd_el);
            for el = 1:this.nelem
                this.GLP(el,:) = reshape(this.ID(this.ELEM(el,:),1)',1,...
                    this.nnd_el);
            end
            
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
                        'fluid',this.mat.fluid);
                if (this.enriched == false)
                    elements(el) = RegularElement_H(...
                                this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                                this.t, emat, this.intOrder,this.GLP(el,:), ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                else
                    elements(el) = EnrichedElement_H(...
                                this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                                this.t, emat, this.intOrder,this.GLP(el,:), ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric);
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
        end


        % -----------------------------------------------------------------
        function seg = initializeDiscontinuitySegArray(~,n)
            seg(n,1) = DiscontinuityElement_H([],[]);
        end

        % -----------------------------------------------------------------
        function seg = initializeDiscontinuitySegment(~,nodeD,matD)
            seg = DiscontinuityElement_H(nodeD,matD);
        end

        %------------------------------------------------------------------
        function mat = createMaterialDataStructure(this)
            mat = struct( ...
                'fluid',this.fluid, ...
                'initialAperture',this.initialAperture);
        end

        % -----------------------------------------------------------------
        function plotPressureAlongSegment(this, Xi, Xf, npts,axisPlot)
            if nargin < 4, npts = 10; end
            EFEMdraw = EFEMDraw(this);
            EFEMdraw.plotPressureAlongSegment(Xi, Xf, npts,axisPlot);
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
                        vertexData(i) = this.matID(el);
                    elseif strcmp(type,'Pressure')
                        p = this.element(el).type.pressureField(X);
                        vertexData(i) = p;
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
            fprintf('\n  Node           Pl\n');
        end

    end
end