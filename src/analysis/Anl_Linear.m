%% Anl_Linear class
%
% This class implements the solution of a static linear analysis. 
%
% Authors:
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
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
