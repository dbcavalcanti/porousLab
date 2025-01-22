%% Anl_Nonlinear Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with linear-elastic analysis.
%
classdef Anl_Linear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Linear(result)
            this = this@Anl('Linear',result);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function process(~,mdl)

            % Compute the global stiffness matrix
            [K, ~, ~, Fext] = mdl.globalMatrices(mdl.U);

            % Set linear system
            A = K(mdl.doffree,mdl.doffree);
            b = Fext(mdl.doffree) - K(mdl.doffree,mdl.doffixed)*mdl.U(mdl.doffixed);

            % Solve linear system
            mdl.U(mdl.doffree) = A\b;

            % Save final result
            for el = 1:mdl.nelem
                gle = mdl.element(el).type.gle;
                mdl.element(el).type.ue = mdl.U(gle);
            end
            
            % Call it again to update the state variables
            mdl.globalMatrices(mdl.U);
            mdl.updateStateVar();
        end

    end
    
end