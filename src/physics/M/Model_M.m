%% Model_M class
%
% Mechanical finite element model.
%
% Each node has 2 degrees of freedom (dof):
% - 2 displacement components (ux,uy)
%
%% Author
% Danilo Cavalcanti
%
%
%% Class definition
classdef Model_M < Model   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        %% Degrees of freedom vectors
        % Matrix with the dofs of each type of each element
        GLU                 = [];
        %% Matrix indicating the Dirichlet BCs 
        SUPP_u              = [];
        %% Matrix with the prescribed BC values 
        PRESCDISPL_u        = [];
        %% Matrix with the Neumann BCs
        LOAD_u              = [];
        %% Additional data
        isPlaneStress       = false;
    end
    
    %% Constructor method
    methods
        function this = Model_M()
            this = this@Model();
            this.ndof_nd = 2;       % Number of dofs per node
            this.physics = 'M';     % Tag with the physics name
            disp("*** Physics: Mechanical");
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function SUPP = dirichletConditionMatrix(this)
            SUPP = this.SUPP_u;  
        end

        %------------------------------------------------------------------
        function LOAD = neumannConditionMatrix(this)
            LOAD = this.LOAD_u;  
        end

        %------------------------------------------------------------------
        function INITCOND = initialConditionMatrix(this)
            INITCOND_u = zeros(this.nnodes,2);
            INITCOND = INITCOND_u;  
        end

        %------------------------------------------------------------------
        function PRESCDISPL = prescribedDirichletMatrix(this)
            PRESCDISPL = this.PRESCDISPL_u;  
        end

        %------------------------------------------------------------------
        function assembleElementDofs(this)

            this.GLU = zeros(this.nelem, this.nnd_el*2);
            for el = 1:this.nelem
                this.GLU(el,:) = reshape(this.ID(this.ELEM(el,:),1:2)',1,...
                    this.nnd_el*2);
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
                        'lc',this.getElementCharacteristicLength(el));
                if (this.enriched == false)
                    elements(el) = RegularElement_M(...
                                this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                                this.t, emat, this.intOrder,this.GLU(el,:), ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress);
                else
                    elements(el) = EnrichedElement_M(...
                                this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                                this.t, emat, this.intOrder,this.GLU(el,:), ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress);
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;

            % Assemble the discontinuity segments to the elements
            if (this.enriched) && ~isempty(this.discontinuitySet)
                nDiscontinuities = size(this.discontinuitySet,1);
                for i = 1:nDiscontinuities
                    % Loop through the segments of this discontinuity
                    nDiscontinuitySeg = size(this.discontinuitySet(i).elemID,1);
                    for j = 1:nDiscontinuitySeg
                        el = this.discontinuitySet.elemID(el);
                        this.element(el).type.addDiscontinuitySegment();
                    end
                end
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
                    elseif strcmp(type,'Sxy')
                        s = this.element(el).type.stressField(X);
                        vertexData(i) = s(3);
                    elseif strcmp(type,'S1')
                        s = this.element(el).type.stressField(X);
                        sp = this.element(el).type.principalStress(s);
                        vertexData(i) = sp(1);
                    elseif strcmp(type,'S2')
                        s = this.element(el).type.stressField(X);
                        sp = this.element(el).type.principalStress(s);
                        vertexData(i) = sp(2);
                    elseif strcmp(type,'Sr')
                        s = this.element(el).type.stressField(X);
                        sp = this.element(el).type.stressCylindrical(s,X);
                        vertexData(i) = sp(1);
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
            fprintf('\n  Node           ux        uy\n');
        end

    end
end