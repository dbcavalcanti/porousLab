%% Anl Class
% This in an abstract class that defines an analysis object.
%
%% Methods
% * *Anl*: Constructor method that initializes the _type_ property if an 
%          input argument is provided.
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef Anl < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = [];
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        % Initializes the analysis object
        function this = Anl(type)
            if nargin > 0
                this.type = type;
            end
        end
    end
end
