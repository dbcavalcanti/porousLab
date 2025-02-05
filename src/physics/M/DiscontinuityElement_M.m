%% Mechanical discontinuity element class
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef DiscontinuityElement_M < DiscontinuityElement       
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = DiscontinuityElement_M(node, mat)
            this = this@DiscontinuityElement(node, mat)
        end
    end
    %% Abstract methods from the super-class
    methods
        function [Ke, Ce, fi, fe, dfidu] = elementData(this, ae)

            Ce = []; fe = []; dfidu = [];

        end
    end
end