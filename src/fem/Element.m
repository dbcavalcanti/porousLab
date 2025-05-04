%% Element Class
% This class represents a finite element in a finite element method (FEM) 
% simulation. It is a handle class, meaning that instances of this class 
% are passed by reference.
%
%% Methods
% * *Element*: Initializes an instance of the _Element_ class.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef Element < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = [];
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        % Constructor method to initialize the _Element_
        function this = Element(type)
            if nargin > 0
                this.type = type;
            end
        end
    end
end
