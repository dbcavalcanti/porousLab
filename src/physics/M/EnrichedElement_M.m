%% EnrichedElement_M class
%
% This class defines a mechanical finite element 
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef EnrichedElement_M < RegularElement_M    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_M(type, node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress)
            this = this@RegularElement_M(type, node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % This function assembles the element matrices and vectors 
        %
        % Output:
        %   Ke : element "stiffness" matrix
        %   Ce : element "damping" matrix
        %   fe : element "internal force" vector
        %
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)

           % [Ke, Ce, fi, fe, dfidu] = this@elementData();
           if isempty(this.discontinuity)
               [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_M(this);
           end
            
        end

        %------------------------------------------------------------------

    end
end