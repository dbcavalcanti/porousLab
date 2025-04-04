function [NODE_D, FRACT, NODE_D_TIPS, W, LEAKOFF] = fractureNodesAndElements(NODE, ELEM, XD, SEGD,aperture,leakoff, TIP, FixedPressureJump, FixedPf, FixedDisplJump)
% 
% This function divides a fracture segment according to the continuum
% finite element mesh.
%
% It does not assume that the continuum mesh is structured.
%

% Initialize the output arrays
NODE_D      = [];
FRACT       = [];
NODE_D_TIPS = [];
W           = [];
LEAKOFF     = [];

nFract = size(SEGD,1);
nElem  = size(ELEM,1);

globalCounter = 1;

% Loop through the fracture segments --------------------------------------
for i = 1:nFract

    % Points that the define the fracture segment
    p1 = XD(SEGD(i,1),:);
    p2 = XD(SEGD(i,2),:);

    % Get the id if these points are crack tips
    id = [0 0];
    if isempty(TIP) ~= 1
        if isempty(find(TIP == SEGD(i,1))) ~= 1
            id(1) = 1;
        end
        if isempty(find(TIP == SEGD(i,2))) ~= 1
            id(2) = 1;
        end
    end
    
    % Initialize local counter
    localCounter = 1;
    isolateNode = [];

    % Loop through the continuum elements ---------------------------------
    for el = 1:nElem 

        % Get the number of edges of the element
        nEdges = size(ELEM,2);
        if (nEdges == 6) || (nEdges == 8)
            nEdges = nEdges/2;
        end

        % Get the coordinates of the element (repeat the first one to close
        % the polygon)
        cX = [NODE(ELEM(el,1:nEdges),1); NODE(ELEM(el,1),1)];
        cY = [NODE(ELEM(el,1:nEdges),2); NODE(ELEM(el,1),2)];

        % Initialize a matrix to store the edges where the discontinuity
        % nodes are located
        % The discontinuity always crosses an element in 2D in two points.
        % Each point is located at a different edge. The discontinuities
        % nodes will be marked when they belong to edges that share a node.
        dEdgeNodes = zeros(nEdges,1);

        % Loop through the edges of the element ---------------------------
        for j = 1:nEdges

            % Points that defined the element edge
            p3 = [cX(j)  , cY(j)];
            p4 = [cX(j+1), cY(j+1)];

            % Evaluate if the segments p1-p2 and p3-p4 intersect
            [flagInt,pint,t12] = intersectionSegment([p1; p2],[p3; p4]);

            % Update the intersection vector points
            if flagInt == true
                NODE_Di(localCounter,:) = pint;
                tNd(localCounter)      = t12;
                dEdgeNodes(j)          = 1;
                localCounter           = localCounter + 1;
            end
        end

        % if sum(dEdgeNodes) > 0
        %     % Find the indices of the 1's in the vector
        %     onesIndices = find(dEdgeNodes == 1);
        % 
        %     % Check if the 1's are in alternate positions
        %     checkIsolateNode = (onesIndices(2) - onesIndices(1)) ~= 2;
        %     isolateNode = [isolateNode;checkIsolateNode;checkIsolateNode];
        % end
    end

    % isolateNode(1) = false;
    % isolateNode(end) = false;

    % Get the unique nodes of the intersection of the fracture "i" with the
    % mesh
    [tNd,I] = sort(tNd);
    NODE_Di = NODE_Di(I,:);
    % isolateNode = isolateNode(I);
    [~,I]   = uniquetol(tNd,1e-9);
    NODE_Di = NODE_Di(I,:);
    % 
    % if length(isolateNode)>3
    %     isolateNodeUnique = zeros(2 + (length(isolateNode)-2)/2,1);
    %     isolateNodeUnique(1) = isolateNode(1);
    %     isolateNodeUnique(end) = isolateNode(end);
    %     ii = 2;
    %     for kk = 2:2:(length(isolateNode)-1)
    %         if isolateNode(kk) || isolateNode(kk+1)
    %             isolateNodeUnique(ii) = true; 
    %         end
    %         ii = ii + 1;
    %     end
    % end

    % Fill the id of the crack tip vector
    NODE_D_TIPSi = zeros(size(NODE_Di,1),4);
    if FixedDisplJump == true
        NODE_D_TIPSi(1,1:2)   = [id(1) , id(1)];
        NODE_D_TIPSi(end,1:2) = [id(2) , id(2)];
    end
    if FixedPressureJump == true
        NODE_D_TIPSi(1,3)   = id(1);
        NODE_D_TIPSi(end,3) = id(2);
    end
    if FixedPf == true
        NODE_D_TIPSi(1,end)   = id(1);
        NODE_D_TIPSi(end,end) = id(2);
    end

    nAddedNd = size(NODE_Di,1);

    if nAddedNd > 1
        FRACT_i = zeros(nAddedNd-1,2);
        for k = 1:nAddedNd-1
            FRACT_i(k,:)  = [globalCounter globalCounter+1];
            W             = [W; aperture(i)];
            LEAKOFF       = [LEAKOFF; leakoff(i)];
            globalCounter = globalCounter + 1;
        end
    end

    % Assemble the global arrays
    NODE_D      = [NODE_D; NODE_Di];
    FRACT       = [FRACT;  FRACT_i];
    NODE_D_TIPS = [NODE_D_TIPS; NODE_D_TIPSi];
    
    globalCounter = globalCounter + 1;

    NODE_Di = [];
    tNd = [];

end

end
