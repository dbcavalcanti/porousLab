%% Discontinuity Class
% This class represents a discontinuity in a finite element mesh. It 
% provides methods for defining, intersecting, and manipulating 
% discontinuities, as well as initializing discontinuity segments for 
% further analysis.
%
%% Methods
% * *setNodeRepel*: Sets the flag to enable/disable the repel process.
% * *setRepelTol*: Sets the node repel tolerance.
% * *setSavePerturbNodes*: Sets the flag to save the perturbed nodes.
% * *intersectMesh*: Intersects the discontinuity with the mesh and 
%                    optionally applies the repel process.
% * *initializeDiscontinuitySegments*: Initializes discontinuity segments
%                                      based on the intersected geometry 
%                                      and material properties.
% * *mat*: Creates the material data strcture
% * *createMaterialDataStructure*: Creates a material data structure for 
%                                  the discontinuity.
% * *getNumberOfDiscontinuitySegments*: Returns the number of discontinuity
%                                       segments.
% * *plotOriginalGeometry*: Plots the original polyline geometry.
% * *plotIntersectedGeometry*: Plots the intersected polyline geometry 
%                              (Xlin).
% * *plotPerturbNodes*: Plots the perturbed nodes.
% * *computeXlin*: Computes the linearized polyline (Xlin) by intersecting
%                  the discontinuity with the mesh.
% * *findElementIDsForXlinSegments*: Finds the element IDs for each segment
%                                    of the linearized polyline (Xlin).
% * *findElementContainingSegment*: Finds the element ID that contains the 
%                                   given segment.
% * *extractEdgesMesh*: Extracts edges from the element connectivity 
%                       matrix.
% * *extractEdgesElement*: Extracts edges from a single element's 
%                          connectivity.
% * *computeNormal*: Computes the normal vector to the discontinuity at a
%                    given point in Xlin.
% * *repelNodes*: Repels nodes in the mesh that are too close to the 
%                 discontinuity, skipping nodes near the first and last 
%                 points of the discontinuity.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef Discontinuity < handle    
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Geometry
        X    = [];  % Nodes defining the polyline
        Xlin = [];  % Nodes of the "linearized" polyline
        PERT = [];  % Nodes from the mesh that were perturbed

        % Geometry tools 
        useRepel = false;          % Flag to enable/disable the repel process
        repelTol = 1.0e-2;         % Node repel tolerance
        savePerturbNodes = false;  % Flag to save the perturbed nodes

        % Topology
        elemID  = [];  % Identification of the element where which discontinuity segment is located
        segment = [];  % Vector with the DiscontinuityElement objects

        % Properties:
        % The properties must be included in the data structure
        % constructed in the createMaterialDataStructure method
        cohesiveLaw         = [];
        fluid               = [];
        initialAperture     = [];
        normalStiffness     = [];
        shearStiffness      = [];
        contactPenalization = [];
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Discontinuity(X,useRepel)
            if nargin == 0, X = []; end
            this.X = X;
            if nargin > 1
                this.useRepel = useRepel;
            end
        end
    end

    %% Public methods
    methods
        %------------------------------------------------------------------
        % Sets the flag to enable/disable the repel process
        function setNodeRepel(this,flag)
            this.useRepel = flag;
        end

        %------------------------------------------------------------------
        % Sets the node repel tolerance
        function setRepelTol(this,tol)
            this.repelTol = tol;
        end

        %------------------------------------------------------------------
        % Sets the flag to save the perturbed nodes
        function setSavePerturbNodes(this,flag)
            this.savePerturbNodes = flag;
        end

        %------------------------------------------------------------------
        % Perform intersection and optionally apply the repel process.
        function intersectMesh(this,model)
            % Compute Xlin using current algorithm
            this.computeXlin(model);

            % Check if there is at least one segment
            if (size(this.Xlin,1) == 1)
                return;
            end

            % Apply the repel process if useRepel is true
            if this.useRepel
                this.repelNodes(model);

                % Recompute Xlin based on the updated mesh
                this.computeXlin(model);
            end

            % Find element IDs for each segment of Xlin
            this.findElementIDsForXlinSegments(model);

            % Check if discontinuity is not fully crossing any element
            if (isempty(this.elemID)) || (sum(this.elemID>0) == 0)
                return;
            end

            % Create discontinuity segments
            this.initializeDiscontinuitySegments(model);
        end

        %------------------------------------------------------------------
        % Initializes discontinuity segments based on the intersected 
        % geometry and material properties
        function initializeDiscontinuitySegments(this,model)
            n = this.getNumberOfDiscontinuitySegments();

            % Create material data structure
            mat = this.createMaterialDataStructure();
            
            % Initialize discontinuity segments according to the physics
            % seg = model.initializeDiscontinuitySegArray(n);
            k = 1;
            for i = 1:size(this.Xlin, 1)-1
                nodes  = [this.Xlin(i,:); this.Xlin(i+1,:)];
                if (this.elemID(i) > 0)
                    seg(k) = model.initializeDiscontinuitySegment(nodes,mat);
                    k = k + 1;
                end
            end
            for i=1:n
                seg(i).initializeIntPoints();
            end
            this.segment = seg;
        end

        %------------------------------------------------------------------
        % Creates the material data strcture
        function mat = createMaterialDataStructure(this)
            mat = struct( ...
                         'fluid',this.fluid,...
                         'cohesiveLaw',this.cohesiveLaw, ...
                         'initialAperture',this.initialAperture, ...
                         'normalStiffness',this.normalStiffness, ...
                         'shearStiffness',this.shearStiffness,...
                         'contactPenalization',this.contactPenalization);
        end

        %------------------------------------------------------------------
        % Gets the number of discontinuities
        function n = getNumberOfDiscontinuitySegments(this)
            n = sum(this.elemID > 0);
        end

        %------------------------------------------------------------------
        % Plot original polyline
        function plotOriginalGeometry(this)
            plot(this.X(:,1), this.X(:,2), '-.xk');
        end

        %------------------------------------------------------------------
        % Plot intersected polyline (Xlin)
        function plotIntersectedGeometry(this)
            for i = 1:size(this.Xlin, 1)-1
                if (this.elemID > 0)
                    seg = [this.Xlin(i,:); this.Xlin(i+1,:)];
                    plot(seg(:,1), seg(:,2), '-.r', 'Marker', 'o', 'MarkerSize', 1.0, 'LineWidth', 1.5);
                end
            end
        end

        %------------------------------------------------------------------
        % Plots perturbed nodes
        function plotPerturbNodes(this)
            if ~isempty(this.PERT)
                plot(this.PERT(:,1), this.PERT(:,2), 'sr');
            end
        end
    end

    %% Private methods
    methods (Access = private)
        %------------------------------------------------------------------
        % Computes the linearized polyline by intersecting the
        % discontinuity with the mesh
        function computeXlin(this,model)
            % Get mesh from model
            NODE = model.NODE;
            ELEM = model.ELEM;

            % Initialize list of intersection points
            intersectionPoints = [];

            % Extract edges from mesh
            edges = this.extractEdgesMesh(ELEM);

            % Iterate over each segment of the polyline
            for i = 1:size(this.X, 1)-1
                % Define current segment of the polyline
                polylineSegment = [this.X(i, :); this.X(i+1, :)];

                % Initialize list of intersection points of this segment
                intersectionPointsSegment = [];
                s = [];

                % Iterate over each edge of the mesh
                for j = 1:size(edges, 1)
                    % Define current edge of the mesh
                    edge = [NODE(edges(j, 1), :); NODE(edges(j, 2), :)];

                    % Compute intersection point between polyline segment and edge
                    [intersect,point] = intersectionSegment(polylineSegment,edge);

                    % If there is an intersection, add the point to list
                    if intersect
                        intersectionPointsSegment = [intersectionPointsSegment; point];
                        sp = sqrt((point(1) - this.X(i,1))^2 + (point(2) - this.X(i,2))^2);
                        s = [s;sp];
                    end
                end

                % Guarantee that the points are ordered
                [s,order] = sort(s);
                intersectionPointsSegment = intersectionPointsSegment(order,:);
                [~,order] = uniquetol(s,1e-9);
                
                % Order the intersection points of the segment
                intersectionPoints = [intersectionPoints;intersectionPointsSegment(order,:)];
            end

            % Store intersection points in Xlin  
            this.Xlin = uniquetol(intersectionPoints,1.0e-9,'ByRows',true);
        end

        %------------------------------------------------------------------
        % Finds the element IDs for each segment of the linearized polyline
        % (Xlin)
        function findElementIDsForXlinSegments(this,model)
            % Get mesh from model
            NODE = model.NODE;
            ELEM = model.ELEM;

            % Initialize list of element IDs for each segment
            elemIDs = [];

            % Iterate over each segment of Xlin
            for i = 1:size(this.Xlin, 1)-1
                % Define current segment of Xlin
                seg = [this.Xlin(i,:);this.Xlin(i+1,:)];

                % Find the element that contains this segment
                eID = this.findElementContainingSegment(NODE,ELEM,seg);

                % Store element ID
                elemIDs = [elemIDs;eID];
            end

            % Store element IDs in elemID
            this.elemID = elemIDs;
        end

        %------------------------------------------------------------------
        % Find the element ID that contains the given segment.
        % Input:
        %   NODE: nx2 matrix of node coordinates
        %   ELEM: mxk matrix of element connectivity
        %   segment: 2x2 matrix defining the segment (two consecutive points in Xlin)
        % Output:
        %   elemID: ID of the element containing the segment
        function eID = findElementContainingSegment(this,NODE,ELEM,segment)
            % Iterate over each element
            for i = 1:size(ELEM, 1)
                count = 0;
                edges = this.extractEdgesElement(ELEM{i});

                % Iterate over each edge of the mesh
                for j = 1:size(edges, 1)
                    % Define current edge of the mesh
                    edge = [NODE(edges(j,1),:); NODE(edges(j,2),:)];

                    % Compute intersection point between polyline segment and edge
                    intersect = intersectionSegment(segment,edge);

                    % If there is an intersection, add the point to the list
                    if intersect
                        count = count + 1;
                    end
                end
                if count == 2
                    eID = i;
                    return
                end
            end

            % If no element is found, return an error
            eID = 0;
        end

        %------------------------------------------------------------------
        % Extract edges from element connectivity matrix.
        function edges = extractEdgesMesh(this,ELEM)
            % Iterate over each element
            edges = [];
            for i = 1:size(ELEM,1)
                elemEdges = this.extractEdgesElement(ELEM{i});
                edges = [edges; elemEdges];
            end

            % Remove duplicate edges
            edges = unique(edges,'rows');
        end

        %------------------------------------------------------------------
        % Extract edges from element connectivity matrix
        function edges = extractEdgesElement(~,elem)
            edges = []; % List of edges
            numNodes = length(elem); % Number of nodes in the element

            for j = 1:numNodes
                % Define edge between node j and node j+1 (wrapping around the first node)
                node1 = elem(j);
                node2 = elem(mod(j,numNodes)+1);

                % Add edge to the list (ensure node1 < node2 to avoid duplicates)
                edges = [edges; sort([node1,node2])];
            end
        end

        %------------------------------------------------------------------
        % Compute normal vector to the discontinuity at a given point in Xlin.
        % Input:
        %   index: Index of the point in Xlin.
        % Output:
        %   normal: Normal vector (unit vector).
        function normal = computeNormal(this,index)
            % Handle cases where X has fewer than 3 points
            if size(this.Xlin, 1) == 2
                % For a straight line, compute the normal directly
                tangent = this.Xlin(2,:) - this.Xlin(1,:); % Tangent vector
                tangent = tangent/norm(tangent); % Normalize
                normal = [-tangent(2),tangent(1)]; % Rotate by 90 degrees
                return;
            end

            % For polylines with 3 or more points, compute normal based on segment
            if index == 1
                seg = this.Xlin(1:2,:); % First segment
            elseif index == size(this.Xlin,1)
                seg = this.Xlin(end-1:end,:); % Last segment
            else
                seg = this.Xlin(index-1:index+1, :); % Middle segment
            end

            % Compute tangent vector of the segment
            tangent = seg(2,:) - seg(1,:);
            tangent = tangent/norm(tangent); % Normalize

            % Compute normal vector (rotate tangent by 90 degrees)
            normal = [-tangent(2),tangent(1)];
        end

        %------------------------------------------------------------------
        % Repel nodes in the mesh that are too close to the discontinuity.
        % Skip repulsion for nodes close to the first and last nodes of the discontinuity.
        function repelNodes(this,model)
            % Get mesh from model
            NODE = model.NODE;
        
            % Get the mean characteristic lengths of the elements associated with each node
            Lc = model.getNodeCharacteristicLength();

            % Get the first and last nodes of the discontinuity
            firstNode = this.X(1,:);
            lastNode  = this.X(end,:);

            % Iterate over each node in the mesh
            for i = 1:size(NODE, 1)
                node = NODE(i,:); % Current mesh node

                % Distance to detect and perturn nodes
                repelDistance = this.repelTol * Lc(i);

                % Check if this node is close to any node in Xlin
                for j = 1:size(this.Xlin, 1)
                    xlinNode = this.Xlin(j,:);      % Current Xlin node
                    distance = norm(node-xlinNode); % Euclidean distance

                    % Skip repulsion if the Xlin node is the first or last node of the discontinuity
                    if isequal(xlinNode, firstNode) || isequal(xlinNode, lastNode)
                        continue;
                    end

                    % If the node is too close, repel it
                    if distance < repelDistance
                        % Compute normal direction to the discontinuity at this point
                        normal = this.computeNormal(j);

                        % Repel node in the normal direction
                        NODE(i,:) = node + repelDistance * normal;
                        if this.savePerturbNodes
                            this.PERT = [this.PERT;NODE(i,:)];
                        end
                        break; % Move to the next mesh node
                    end
                end
            end

            % Update data in model object
            model.NODE = NODE;
        end
    end
end
