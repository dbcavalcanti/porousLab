%% EnrichedElement class
%
% This class defines a finite element
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%%%
% Initially prepared for the course CIV 2801 - Fundamentos de Computação
% Gráfica, 2022, second term, Department of Civil Engineering, PUC-Rio.
%
classdef Element < handle

    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        type     = [];          
    end

    %% Constructor method
    methods
        function this = Element(type)
            if nargin > 0
                this.type = type;
            end
        end
    end
end