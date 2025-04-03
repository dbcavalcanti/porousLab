%% DiscontinuityElement_H Class
% This class defines a hydraulic discontinuity element for modeling 
% discontinuities in porous media mechanics. It inherits from the 
% _DiscontinuityElement_ base class and provides specific implementations 
% for handling hydraulic discontinuities.
%
%% Methods
% * *initializeIntPoints*: Initializes the integration points for the 
%                          discontinuity element. Retrieves integration 
%                          points' coordinates and weights, and creates 
%                          _IntPoint_ objects with the associated material 
%                          model.
% * *elementData*: Computes the element stiffness matrix, internal force 
%                  vector, and other element data based on the input 
%                  displacement vector. Performs numerical integration 
%                  over the integration points.
% 
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class Definition
classdef DiscontinuityElement_H < DiscontinuityElement    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement_H(node, mat)
            this = this@DiscontinuityElement(node, mat)
            this.ndof = 2;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initializes the integration points for the element obtaining the
        % coordinates and weights
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(1);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = MaterialDiscontinuity_H(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % Computes the element stiffness matrix, internal force vector and
        % other optional outputs using numerical integration over the
        % element
        % 
        % Outputs:
        %   Ke    - Element stiffness matrix.
        %   Ce    - Element damping matrix.
        %   fi    - Internal force vector.
        %   fe    - External force vector.
        %   dfidu - Derivative of internal force with respect to 
        %           displacement.
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)
            
            % Declare output matrices that won't be used
            fi = []; fe = []; dfidu = [];

            % Initialize the matrices for the numerical integration
            Ke = zeros(this.ndof,this.ndof);
            Ce = zeros(this.ndof,this.ndof);

            % Initialize output matrices
            for i = 1:this.nIntPoints

                % Get the shape function matrix
                N  = this.shape.shapeFnc(this.intPoint(i).X);

                % Compute the B matrix at the int. point and the detJ
                [dN, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Compute the permeability matrix
                kh = this.intPoint(i).constitutiveMdl.longitudinalPermeability();

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibility();

                % Numerical integration term. The determinant is ld/2.
                c = detJ * this.intPoint(i).w * this.t;

                % Compute the permeability matrix
                Ke = Ke + dN' * kh * dN * c;

                % Compute the compressibility matrix
                Ce = Ce + N' * comp * N * c;

            end
        end
    end
end