%% Mechanical discontinuity element class
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef DiscontinuityElement_M < DiscontinuityElement    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        stretchingMode  = false;
        relRotationMode = false;
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement_M(node, mat)
            this = this@DiscontinuityElement(node, mat)
            this.ndof = 2;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        function addStretchingMode(this,flag)
            this.stretchingMode = flag;
            this.ndof = this.ndof + 1;
        end

        %------------------------------------------------------------------
        function addRelRotationMode(this,flag)
            this.relRotationMode = flag;
            this.ndof = this.ndof + 1;
        end

        %------------------------------------------------------------------
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints(1);

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints
                constModel = MaterialDiscontinuity_M(this.mat);
                intPts(i) = IntPoint(X(:,i),w(i), constModel);
                intPts(i).initializeMechanicalAnalysisModel('Interface');
            end
            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        function [Ke, Ce, fi, fe, dfidu] = elementData(this, ae)
            
            % Declare output matrices that won't be used
            Ce = []; fe = []; dfidu = [];

            % Initialize output matrices


        end
    end
end