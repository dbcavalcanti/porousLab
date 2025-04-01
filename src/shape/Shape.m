%% Shape Class
% Abstract class defining the type of the element and providing a 
% framework for implementing shape functions, jacobian computations, 
% and other related operations needed for finite element analysis.
%
%% Methods
% This class provides the following methods:
%
% * *shapeFnc*: Evaluates the shape function at a given point.
% * *shapeFncMtrx*: Returns the shape function matrix.
% * *shapeFncDrv*: Computes the shape function derivatives matrix.
% * *JacobianMtrx*: Computes the Jacobian matrix.
% * *detJacobian*: Computes the determinant of the Jacobian matrix.
% * *BMatrix*: Constructs the strain-displacement matrix.
% * *coordNaturalToCartesian*: Transforms a point to the global Cartesian 
%                              coordinate system.
% * *coordCartesianToNatural*: Transforms a point to the natural coordinate 
%                              system.
% * *getIntegrationPoints*: Retrieves the integration points.
% * *computeCentroid*: Computes the centroid of an element.
% * *axisSymmetricFactor*: Computes the axisymmetric factor.
%
%% Author
% - Danilo Cavalcanti
%
%% Version History
% - Version 1.00: Initial version (January 2023).
%
%% Class Definition
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
        % Compute the centroid af an element
        function Xc = computeCentroid(~,NODE)
            x = NODE(1:3,1);
            y = NODE(1:3,2);
            polyin = polyshape({x},{y});
            [xc,yc] = centroid(polyin);
            Xc = [xc yc];
        end

        % Compute the axisymmetric factor
        function af = axisSymmetricFactor(~,N,X)
            r = N*X(:,1);
            af = 2.0 * r * pi;
        end
    end
end
