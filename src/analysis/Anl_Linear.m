%% Anl_Linear Class
% This class implements the solution of a static linear analysis. It 
% inherits from the base class _Anl_ and provides functionality for 
% performing linear static analysis on a given model. 
%
%% Methods
% * *run*: Executes the linear analysis for the given model object _mdl_
%
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
% 
%% Class definition
classdef Anl_Linear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
    end

    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Linear()
            this = this@Anl('Linear');
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        function run(~,mdl)
            % Initialize model object
            mdl.preComputations();
            
            disp("*** Performing linear analysis...")

            % Compute global stiffness matrix
            [K,~,~,Fext] = mdl.globalMatrices(mdl.U);

            % Set linear system
            A = K(mdl.doffree,mdl.doffree);
            b = Fext(mdl.doffree) - K(mdl.doffree,mdl.doffixed) * mdl.U(mdl.doffixed);

            % Solve linear system
            mdl.U(mdl.doffree) = A\b;

            % Save final result
            for i = 1:mdl.nelem
                gle = mdl.element(i).type.gle;
                mdl.element(i).type.ue = mdl.U(gle);
            end

            % Call it again to update state variables
            mdl.globalMatrices(mdl.U);
            mdl.updateStateVar();

            disp("*** Analysis completed!");
        end
    end
end
