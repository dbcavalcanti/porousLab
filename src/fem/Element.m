%% EnrichedElement class
%
% This class defines a finite element.
%
classdef Element < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = [];
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Element(type)
            if nargin > 0
                this.type = type;
            end
        end
    end
end
