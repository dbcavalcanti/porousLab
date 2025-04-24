%% EnrichedElement_H class
% This class extends the _RegularElement_H_ class to define a finite 
% element for single-phase fluid flow that incorporates enriched elements 
% to handle discontinuities. It provides methods to compute element data, 
% manage discontinuities, and calculate enriched degrees of freedom.
%
%% Methods
% * *elementData*: Computes the element data (stiffness matrix, damping 
%                  matrix, internal force vector, external force vector, 
%                  and derivative of internal force vector) based on 
%                  whether the element has discontinuities.
% * *getDiscontinuitiesData*: Computes the stiffness matrix and force 
%                             vector contributions from the 
%                             discontinuities.
% * *getNumberEnrichedDofs*: Returns the total number of enriched degrees 
%                            of freedom.
% * *getNumberOfDiscontinuities*: Returns the number of discontinuities 
%                                 associated with the element.
% * *getNumberOfDofPerDiscontinuity*: Returns the number of degrees of 
%                                     freedom per discontinuity, 
%                                     considering the enabled enrichment 
%                                     modes.
% * *addDiscontinuitySegment*: Adds a discontinuity segment to the element.
%% Author
% Danilo Cavalcanti
%
%% Version History
% Version 1.00.
%
%% Class definition
classdef EnrichedElement_H < RegularElement_H   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_H(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric)
            this = this@RegularElement_H(node, elem, t, ...
                mat, intOrder, glu, massLumping, lumpStrategy, ...
                isAxisSymmetric);
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes the element data for the current element based on wether
        % the element contains a discontinuity or not.
        % 
        % Outputs:
        %   Ke    - Element stiffness matrix.
        %   Ce    - Element damping matrix.
        %   fi    - Internal force vector.
        %   fe    - External force vector.
        %   dfidu - Derivative of internal force with respect to 
        %           displacement.
        function [Ke, Ce, fi, fe, dfidu] = elementData(this)

           % Get the continuum contribution
           [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_H(this);
           if ~isempty(this.discontinuity)
               [Kd, Sd, Tcd] = getDiscontinuitiesData(this);
               Ke = Ke + Tcd' * Kd * Tcd;
               Ce = Ce + Tcd' * Sd * Tcd;   
           end
        end

        %------------------------------------------------------------------
        % Computes the stiffness matrix and force vector contributions
        % from discontinuities in the enriched element
        function [Kd, Sd, Tcd] = getDiscontinuitiesData(this)

            nEnrDofs          = this.getNumberEnrichedDofs();
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();

            % Initialize the output data 
            Kd = zeros(nEnrDofs,nEnrDofs);
            Sd = zeros(nEnrDofs,nEnrDofs);

            % Condensation matrix
            Tcd = zeros(2*nDiscontinuities,this.nglp);

            % Loop through the discontinuities
            k = 1;
            for i = 1:nDiscontinuities

                % Dofs associated with this discontinuity segment
                dofs = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;

                % Get the discontinuity data
                [Kdi,Sdi,~,~,~] = this.discontinuity(i).elementData();

                % Assemble the contribution of this discontinuity
                Kd(dofs,dofs) = Kdi;
                Sd(dofs,dofs) = Sdi;

                % Loop through the nodes of the discontinuity to fill the
                % condensation matrix
                for j = 1:2
                    X = this.discontinuity(i).node(j,:);
                    Xn = this.shape.coordCartesianToNatural(this.node,X);
                    Np = this.shape.shapeFncMtrx(Xn);
                    Tcd(k,:) = Np;
                    k = k + 1;
                end
                
            end
        end

        %------------------------------------------------------------------
        % Gets the number of enriched degrees of freedom
        function nEnrDof = getNumberEnrichedDofs(this)
            nEnrDof = this.getNumberOfDiscontinuities();
            nEnrDof = nEnrDof * this.getNumberOfDofPerDiscontinuity();
        end

        %------------------------------------------------------------------
        % Obtain the number of discontinuities
        function n = getNumberOfDiscontinuities(this)
            n = size(this.discontinuity,1);
        end

        %------------------------------------------------------------------
        % Calculates the number of degrees of freedom per discontinuity for
        % the enriched element
        function n = getNumberOfDofPerDiscontinuity(this)
            n = 2;  
        end

        %------------------------------------------------------------------
        % Adds a discontinuity segment to the element
        function addDiscontinuitySegment(this,dseg)
            this.discontinuity = [this.discontinuity; dseg];
        end
        %------------------------------------------------------------------
    end
end