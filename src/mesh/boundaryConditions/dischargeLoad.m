function [LOAD] = dischargeLoad(p,Xref,dir,NODE,ELEM,LOAD)
% Apply pressure loads 
%
% Valid for quadrilateral domains
%
% Inputs:
%  p: pressure value (kPa)
%  Xref: reference point at the border that the pressure is being applied
%  dir: direction of the pressure load (1 = 'x', 2 = 'y')
%  NODE: nodes coordinates
%  ELEM: elements connectivity
%  LOAD: matrix with the nodal loads in the x and y direction
%
% Output: updated LOAD matrix
%
% -------------------------------------------------------------------------

% Number of elements
nElem = size(ELEM,1);

for el = 1:nElem 

    % Get the number of edges of the element
    nEdges = size(ELEM,2);

    % Get the coordinates of the element
    cX = [NODE(ELEM(el,:),1); NODE(ELEM(el,1),1)];
    cY = [NODE(ELEM(el,:),2); NODE(ELEM(el,1),2)];

    % Get the nodes of the borders
    NdBorders = [ELEM(el,:), ELEM(el,1)];

    % Loop through the edges of the element ---------------------------
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
        if norm(edge-Xref(dir)) < 1.0e-12
            
            % Compute the length of the edge
            dx = edgeX(2) - edgeX(1);
            dy = edgeY(2) - edgeY(1);
            l = sqrt(dx*dx + dy*dy);

            % Equivalent nodal load
            feq = 0.5*p*l;

            % id of the nodes of the edge
            idNds = [NdBorders(j); NdBorders(j+1)];

            % Add contribution to the LOAD matrix
            LOAD(idNds) = LOAD(idNds) + feq;
        end

    end
end