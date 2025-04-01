%% Shape class
%
% This is an abstract class to define the type of the element
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version:January 2023
%
classdef Shape < handle

    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        type     = [];          
    end

    %% Constructor method
    methods
        function this = Shape(type)
            if nargin > 0
                this.type = type;
            end
        end
    end

    %% Abstract methods
    methods (Abstract)
         % Evaluate the shape function at a given point X
         N = shapeFnc(this,Xn)

         % Get the shape function matrix
         N = shapeFncMtrx(this,Xn)

         % Get the shape function derivatives matrix
         dNdxi = shapeFncDrv(this,Xn)

         % Compute the jacobian matrix
         J = JacobianMtrx(this,X,Xn)

         % Compute the determinant of the jacobian
         detJ = detJacobian(this,X,Xn)

         % Compute the strain-displacement matrix
         [B,detJ] = BMatrix(this,X,Xn)

         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         X = coordNaturalToCartesian(this,NODE,Xn)

         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         Xn = coordCartesianToNatural(this,NODE,X)

         % Get the integration points
         [X,W,n] = getIntegrationPoints(this,intOrder,elem)
    end
    
    %% Public methods
    methods
        function Xc = computeCentroid(~,NODE)
            x = NODE(1:3,1);
            y = NODE(1:3,2);
            polyin = polyshape({x},{y});
            [xc,yc] = centroid(polyin);
            Xc = [xc yc];
        end

        function af = axisSymmetricFactor(~,N,X)
            r = N*X(:,1);
            af = 2.0 * r * pi;
        end
    end
end
