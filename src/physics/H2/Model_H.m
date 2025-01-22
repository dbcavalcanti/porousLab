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
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function createNodeDofIdMtrx(this)
            % Initialize the ID matrix and the number of fixed dof
            this.ID = zeros(this.nnodes,this.ndof_nd);
            this.ndoffixed = 0;

            SUPP = [this.SUPP_p];
            
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

            % Vector with all regular dofs
            this.Dof = unique(this.GLP);
            
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
        
        %------------------------------------------------------------------
        function initializeDisplacementVct(this)

            % Initialize the displacement vector 
            this.U = zeros(this.ndof,1);

            % Set the initial values
            for i = 1:this.nnodes
                this.U(this.ID(i,1)) = this.INITCOND_p(i,1);
            end

            % Set the prescribed values
            for i = 1:this.nnodes
                if (this.SUPP_p(i,1) == 1.0)
                    this.U(this.ID(i,1)) = this.PRESCDISPL_p(i,1);
                end
            end

            % Save initial dofs to the elements
            for el = 1 : this.nelem
                this.element(el).type.ue = this.U(this.element(el).type.gle);
            end

        end
 
        %------------------------------------------------------------------
        function Fref = addNodalLoad(this,Fref)
            for i = 1:this.nnodes
                Fref(this.ID(i,1)) = Fref(this.ID(i,1)) + this.LOAD_p(i,1);
            end
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
                        vertexData(i) = 0.0;
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