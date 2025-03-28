%% Model_M class
%
% Mechanical finite element model.
%
% Each node has 2 degrees of freedom (dof):
% - 2 displacement components (ux,uy)
%
%% Author
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class definition
classdef Model_M < Model   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        %% Model data
        isPlaneStress       = false;
        %% Embedded related data
        addStretchingMode   = false;
        addRelRotationMode  = false;
    end
    
    %% Constructor method
    methods
        function this = Model_M(printFlag)
            if nargin == 0, printFlag = true; end
            this = this@Model();
            this.ndof_nd = 2;       % Number of dofs per node
            this.physics = 'M';     % Tag with the physics name
            if (printFlag)
                disp("*** Physics: Mechanical");
            end 
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        function setMaterial(this,porousMedia)
            if nargin < 2
                disp('Error in setMaterial: insuficient number of inputs.');
                disp('The HM physics requires two attributes: porousMedia.');
                error('Error in setMaterial.');
            end
            if ~isa(porousMedia,'PorousMedia')
                disp('Error in setMaterial: porousMedia is not a PorousMedia object.');
                error('Error in setMaterial.');
            end
            this.mat  = struct( ...
            'porousMedia',porousMedia);
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
                dof_e = this.getElementDofs(el,[1,2]);
                if (this.enriched == false)
                    elements(el) = RegularElement_M(...
                                this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                                this.t, emat, this.intOrder,dof_e, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress);
                else
                    elements(el) = EnrichedElement_M(...
                                this.type,this.NODE(this.ELEM(el,:),:), this.ELEM(el,:),...
                                this.t, emat, this.intOrder,dof_e, ...
                                this.massLumping, this.lumpStrategy, this.isAxisSymmetric, ...
                                this.isPlaneStress,this.addRelRotationMode,this.addStretchingMode);
                end
                elements(el).type.initializeIntPoints();
            end
            this.element = elements;
            
        end

        % -----------------------------------------------------------------
        function setDisplacementDirichletBCAtNode(this, nodeId, value)
            this.setDirichletBCAtNode(nodeId, [1,2], value);
        end

        % -----------------------------------------------------------------
        function setDisplacementDirichletBCAtPoint(this, X, value)
            this.setDirichletBCAtPoint(X, [1,2], value);
        end

        % -----------------------------------------------------------------
        function setDisplacementDirichletBCAtBorder(this, border, value)
            this.setDirichletBCAtBorder(border, [1,2], value);
        end

        % -----------------------------------------------------------------
        function addLoadAtNode(this, nodeId, value)
            this.setNeumannBCAtNode(nodeId, [1,2], value);
        end

        % -----------------------------------------------------------------
        function addLoadAtPoint(this, X, value)
            this.setNeumannBCAtPoint(X, [1,2], value);
        end

        % -----------------------------------------------------------------
        function addLoadAtBorder(this, border, dir, p)
            % Get the nodes at the given border
            if strcmp(border,'left')
                ref = min(this.NODE(:,1));
            elseif strcmp(border,'right')
                ref = max(this.NODE(:,1));
            elseif strcmp(border,'top')
                ref = max(this.NODE(:,2));
            elseif strcmp(border,'bottom')
                ref = min(this.NODE(:,2));
            else
                disp('Warning: non-supported border.');
                disp('Available borders tag: ''left'',''right'', ''top'',''bottom''');
            end
            % Get number of linear interpolation points
            nLinNodes = this.nnd_el;
            quadMesh  = false;
            if strcmp(this.type,'LST') || strcmp(this.type,'ISOQ8')
                nLinNodes = nLinNodes / 2;
                quadMesh  = true;
            end
            for el = 1:this.nelem 
                % Get the number of edges of the element
                nEdges = nLinNodes;
            
                % Get the coordinates of the element
                cX = [this.NODE(this.ELEM(el,1:nLinNodes),1); this.NODE(this.ELEM(el,1),1)];
                cY = [this.NODE(this.ELEM(el,1:nLinNodes),2); this.NODE(this.ELEM(el,1),2)];
            
                % Get the nodes of the borders
                NdBorders = [this.ELEM(el,1:nLinNodes), this.ELEM(el,1)];
            
                % Loop through the edges of the element
                for j = 1:nEdges
            
                    % coordinates of the edge
                    edgeX = [cX(j) , cX(j+1)];
                    edgeY = [cY(j) , cY(j+1)];
            
                    % select the edge
                    if dir == 1
                        edge = edgeX;
                    elseif dir == 2
                        edge = edgeY;
                    end
            
                    % check if the edge belong to the boundary
                    if norm(edge-ref) < 1.0e-12
                        
                        % Compute the length of the edge
                        dx = edgeX(2) - edgeX(1);
                        dy = edgeY(2) - edgeY(1);
                        l = sqrt(dx*dx + dy*dy);
            
                        % id of the nodes of the edge
                        idNds = [NdBorders(j); NdBorders(j+1)];
            
                        % Equivalent nodal load
                        if quadMesh == false
                            feq = [0.5*p*l;0.5*p*l];
                        else
                            feq = [p*l;p*l;4.0*p*l]/6.0;
                            idNds = [idNds; this.ELEM(el,j+nLinNodes)];
                        end
            
                        % Add contribution to the LOAD matrix
                        this.LOAD(idNds,dir) = this.LOAD(idNds,dir) + feq;
                    end
                end
            end
        end

        % -----------------------------------------------------------------
        function seg = initializeDiscontinuitySegArray(~,n)
            seg(n,1) = DiscontinuityElement_M([],[]);
        end

        % -----------------------------------------------------------------
        function seg = initializeDiscontinuitySegment(~,nodeD,matD)
            seg = DiscontinuityElement_M(nodeD,matD);
        end

        % -----------------------------------------------------------------
        function addDiscontinuityData(this,additionalData)
            this.addRelRotationMode = additionalData.addRelRotationMode;
            this.addStretchingMode  = additionalData.addStretchingMode;
        end

        %------------------------------------------------------------------
        function initializeDiscontinuitySegments(this)
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                nDiscontinuitySeg = this.discontinuitySet(i).getNumberOfDiscontinuitySegments();
                for j = 1:nDiscontinuitySeg
                    this.discontinuitySet(i).segment(j).t = this.t;
                    this.discontinuitySet(i).segment(j).addStretchingMode(this.addStretchingMode);
                    this.discontinuitySet(i).segment(j).addRelRotationMode(this.addRelRotationMode);
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
    end
    %% Static methods
    methods (Static)

        % -----------------------------------------------------------------
        function printResultsHeader()
            fprintf('\n  Node           ux        uy\n');
        end

    end
end