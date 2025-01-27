function [NODE_quad, ELEM_quad] = convertToQuadraticMesh(NODE, ELEM)
    % Function to convert a linear mesh into a quadratic mesh
    % Supports Q4 -> Q8 and T3 -> T6 elements
    %
    % Inputs:
    %   NODE - Matrix with x and y coordinates of the nodes (size Nx2)
    %   ELEM - Connectivity matrix (size Mx4 for Q4 or Mx3 for T3)
    %
    % Outputs:
    %   NODE_quad - Updated NODE matrix with mid-edge nodes
    %   ELEM_quad - Updated ELEM matrix with quadratic connectivity

    % Initialize
    num_nodes = size(NODE, 1);
    edge_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    new_nodes = [];

    % Determine element type (triangular or quadrilateral)
    num_nodes_per_elem = size(ELEM, 2);
    if num_nodes_per_elem == 3
        % Triangular element (T3 -> T6)
        num_edges = 3;
        edge_pairs = [1, 2; 2, 3; 3, 1];
    elseif num_nodes_per_elem == 4
        % Quadrilateral element (Q4 -> Q8)
        num_edges = 4;
        edge_pairs = [1, 2; 2, 3; 3, 4; 4, 1];
    else
        error('Unsupported element type. ELEM should have 3 or 4 columns.');
    end

    % Initialize new ELEM matrix
    ELEM_quad = zeros(size(ELEM, 1), num_nodes_per_elem + num_edges);

    % Loop through elements to add mid-edge nodes
    for elem_idx = 1:size(ELEM, 1)
        elem = ELEM(elem_idx, :);
        new_elem = elem;

        for edge_idx = 1:num_edges
            % Get the edge node indices
            n1 = elem(edge_pairs(edge_idx, 1));
            n2 = elem(edge_pairs(edge_idx, 2));

            % Sort to ensure consistent edge identification
            edge_key = sprintf('%d-%d', min(n1, n2), max(n1, n2));

            if isKey(edge_map, edge_key)
                % Mid-edge node already exists
                mid_node = edge_map(edge_key);
            else
                % Create new mid-edge node
                mid_coord = (NODE(n1, :) + NODE(n2, :)) / 2;
                new_nodes = [new_nodes; mid_coord];
                mid_node = num_nodes + size(new_nodes, 1);
                edge_map(edge_key) = mid_node;
            end

            % Add mid-edge node to the new element connectivity
            new_elem(num_nodes_per_elem + edge_idx) = mid_node;
        end

        % Update quadratic ELEM matrix
        ELEM_quad(elem_idx, :) = new_elem;
    end

    % Update NODE matrix with new mid-edge nodes
    NODE_quad = [NODE; new_nodes];
end
