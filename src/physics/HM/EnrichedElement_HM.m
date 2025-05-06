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
classdef EnrichedElement_HM < RegularElement_HM   
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        discontinuity = [];
        addStretchingMode  = false;
        addRelRotationMode = false
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = EnrichedElement_HM(node, elem, t, ...
                mat, intOrder, glu, glp, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress, ...
                addRelRotationMode,addStretchingMode)
            this = this@RegularElement_HM(node, elem, t, ...
                mat, intOrder, glu, glp, massLumping, lumpStrategy, ...
                isAxisSymmetric,isPlaneStress);
            this.addStretchingMode  = addStretchingMode;
            this.addRelRotationMode = addRelRotationMode;
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
           
           if isempty(this.discontinuity)
               % Get the continuum contribution
               [Ke, Ce, fi, fe, dfidu] = elementData@RegularElement_HM(this);
           else
               % Compute contribution of the continuum
               [Kuu, Kua, Kau, Kaa, Quc, Quj, Qac, Qaj, Qca, Qja, Hcc, Hcj, Hjc, Hjj, Scc, Scj, Sjc, Sjj, fiu, fia, fic, fij, feu, fea, fec, fej] = this.fillElementSubData();
               
               % Compute contribution of the discontinuity
               [fid, Kdd, Qad, Hdd, Sdd, Lcc, Lcj, Lcd, Ljc, Ljj, Ljd, Ldc, Ldj, Ldd, Tdc, Tdj] = getDiscontinuitiesData(this);
               
               % Condense the internal pressure dofs
               Hcc = Hcc + Tdc' * Hdd * Tdc;
               Hcj = Hcj + Tdc' * Hdd * Tdj;
               Hjc = Hjc + Tdj' * Hdd * Tdc;
               Hjj = Hjj + Tdj' * Hdd * Tdj;
               
               Scc = Scc + Tdc' * Sdd * Tdc;
               Scj = Scj + Tdc' * Sdd * Tdj;
               Sjc = Sjc + Tdj' * Sdd * Tdc;
               Sjj = Sjj + Tdj' * Sdd * Tdj;

               Lcc = Lcc + Tdc' * Ldd * Tdc - Lcd * Tdc - Tdc' * Ldc;
               Lcj = Lcj + Tdc' * Ldd * Tdj - Lcd * Tdj - Tdc' * Ldj;
               Ljc = Ljc + Tdj' * Ldd * Tdc - Ljd * Tdc - Tdj' * Ldc;
               Ljj = Ljj + Tdj' * Ldd * Tdj - Ljd * Tdj - Tdj' * Ldj;

               Qac = Qac + Qad * Tdc;
               Qaj = Qaj + Qad * Tdj;
               Qca = Qca + Tdc' * Qad';
               Qja = Qja + Tdj' * Qad';

               % Add contribution of the coupling matrices
               Hcc = Hcc + Lcc;
               Hcj = Hcj + Lcj;
               Hjc = Hjc + Ljc;
               Hjj = Hjj + Ljj;

               % Add contribution of the discontinuity stiffness
               fia = fia + fid;
               Kaa = Kaa + Kdd;

               % Number of displacement and pressure jump dofs
               ndofa = size(Kaa,1);
               ndofj = size(Hjj,1);

               % Create auxiliary zero sub-matrices
               Ouu = zeros(this.nglu, this.nglu);
               Ouc = zeros(this.nglu, this.nglp);
               Oua = zeros(this.nglu, ndofa);
               Ouj = zeros(this.nglu, ndofj);

               Ocu = zeros(this.nglp, this.nglu);
               Occ = zeros(this.nglp, this.nglp);
               Oca = zeros(this.nglp, ndofa);
               Ocj = zeros(this.nglp, ndofj);

               Oau = zeros(ndofa, this.nglu);
               Oac = zeros(ndofa, this.nglp);
               Oaa = zeros(ndofa, ndofa);
               Oaj = zeros(ndofa, ndofj);

               Oju = zeros(ndofj, this.nglu);
               Ojc = zeros(ndofj, this.nglp);
               Oja = zeros(ndofj, ndofa);
               Ojj = zeros(ndofj, ndofj);
               
               % Assemble the element matrices
               dfidu = [ Kuu, Ouc, Kua, Ouj;
                         Ocu, Occ, Oca, Ocj;
                         Kau, Oac, Kaa, Oaj;
                         Oju, Ojc, Oja, Ojj];

               Ke = [ Ouu, -Quc, Oua, -Quj;
                      Ocu,  Hcc, Oca,  Hcj;
                      Oau, -Qac, Oaa, -Qaj;
                      Oju,  Hjc, Oja,  Hjj];
               
               Ce = [ Ouu,  Ouc, Oua, Ouj;
                      Quc', Scc, Qca, Scj;
                      Oau,  Oac, Oaa, Oaj;
                      Quj', Sjc, Qja, Sjj];
               
               fi = [fiu; fic; fia; fij];
               
               fe = [feu; fec; fea; fej];
           end
        end

        %------------------------------------------------------------------
        % Computes and assembles the sub-matrices and sub-vectors for an
        % enriched finite element
        function [Kuu, Kua, Kau, Kaa, Quc, Quj, Qac, Qaj, Qca, Qja, Hcc, Hcj, Hjc, Hjj, Scc, Scj, Sjc, Sjj, fiu, fia, fic, fij, feu, fea, fec, fej] = fillElementSubData(this)

            % Number of displacement and pressure jump dofs
            nd = this.getNumberOfDiscontinuities();
            ndofa = nd * this.getNumberOfDisplacementDofPerDiscontinuity();
            ndofj = nd; % Each discontinuity has a pressure jump dof

            % Initialize stiffness sub-matrices
            Kuu = zeros(this.nglu, this.nglu);
            Kua = zeros(this.nglu, ndofa);
            Kau = zeros(ndofa, this.nglu);
            Kaa = zeros(ndofa, ndofa);

            % Initialize the hydro-mechanical coupling sub-matrices
            Quc = zeros(this.nglu, this.nglp);
            Quj = zeros(this.nglu, ndofj);
            Qac = zeros(ndofa, this.nglp);
            Qaj = zeros(ndofa, ndofj);
            Qca = zeros(this.nglp, ndofa);
            Qja = zeros(ndofj, ndofa);
            
            % Initialize fluid-flow sub-matrices
            Hcc = zeros(this.nglp, this.nglp);
            Hcj = zeros(this.nglp, ndofj);
            Hjc = zeros(ndofj, this.nglp);
            Hjj = zeros(ndofj, ndofj);

            % Initialize compressibity sub-matrices
            Scc = zeros(this.nglp , this.nglp);
            Scj = zeros(this.nglp , ndofj);
            Sjc = zeros(ndofj, this.nglp);
            Sjj = zeros(ndofj, ndofj);

            % Initialize external force vector
            feu = zeros(this.nglu, 1);
            fec = zeros(this.nglp, 1);
            fea = zeros(ndofa, 1);
            fej = zeros(ndofj, 1);

            % Initialize internal force vector
            fiu = zeros(this.nglu, 1);
            fic = zeros(this.nglp, 1);
            fia = zeros(ndofa, 1);
            fij = zeros(ndofj, 1);
            
            % Vector of the nodal dofs
            u  = this.getNodalDisplacement();
            pc = this.getNodalPressure();
            a  = this.getDisplacementJumpDofs(ndofa);
            pj = this.getPressureJumpDofs(ndofa);

            % Initialize 2D identity vector
            m = [1.0 ; 1.0 ; 1.0 ; 0.0];

            % Numerical integration of the sub-matrices
            for i = 1:this.nIntPoints

                % Shape function matrix
                Npc = this.shape.shapeFncMtrx(this.intPoint(i).X);

                % Cartesian coordinates of the integration point
                Xcar = this.shape.coordNaturalToCartesian(this.node,this.intPoint(i).X);

                % Enriched shape function matrix
                Npj = this.enrichedShapeFncMtrx(Npc, Xcar);
               
                % Compute the B matrix at the int. point and the detJ
                [Bp, detJ] = this.shape.dNdxMatrix(this.node,this.intPoint(i).X);

                % Assemble the B-matrix for the mechanical part
                Bu = this.shape.BMatrix(Bp);

                % Compute the G matrix
                Gp = this.Gmatrix(Bp);

                % Get kinematic enriched matrix
                Gr = this.kinematicEnrichment(Bu);

                % Get the static enriched matrix (TEMP)
                Gv = this.kinematicEnrichment(Bu);

                % Pressure values at the integration point
                pIP = Npc * pc + Npj * pj;
        
                % Compute the permeability matrix
                kh = this.intPoint(i).constitutiveMdl.permeabilityTensor();

                % Get compressibility coefficient
                comp = this.intPoint(i).constitutiveMdl.compressibilityCoeff();

                % Get Biot's coefficient
                biot = this.intPoint(i).constitutiveMdl.biotCoeff();

                % Compute the strain vector
                this.intPoint(i).strain = Bu * u + Gr * a;

                % Compute the stress vector and the constitutive matrix
                [stress,Duu] = this.intPoint(i).mechanicalLaw();
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
                if this.isAxisSymmetric
                    c = c * this.shape.axisSymmetricFactor(Np,this.node);
                end

                % Internal force vector
                fiu = fiu + Bu' * stress * c;
                fia = fia + Gv' * stress * c;

                % Compute the stiffness sub-matrices
                Kuu = Kuu + Bu' * Duu * Bu * c;
                Kua = Kua + Bu' * Duu * Gr * c;
                Kau = Kau + Gv' * Duu * Bu * c;
                Kaa = Kaa + Gv' * Duu * Gr * c;

                % Compute the hydro-mechanical coupling sub-matrices
                Quc = Quc + Bu'  * m  * biot * Npc * c;
                Quj = Quj + Bu'  * m  * biot * Npj * c;
                Qac = Qac + Gv'  * m  * biot * Npc * c;
                Qaj = Qaj + Gv'  * m  * biot * Npj * c;
                Qca = Qca + Npc' * m' * biot * Gr  * c;
                Qja = Qja + Npj' * m' * biot * Gr  * c;
        
                % Compute permeability sub-matrices
                Hcc = Hcc + Bp' * kh * Bp * c;
                Hcj = Hcj + Bp' * kh * Gp * c;
                Hjc = Hjc + Gp' * kh * Bp * c;
                Hjj = Hjj + Gp' * kh * Gp * c;

                % Compute compressibility matrices
                Scc = Scc + Npc' * comp * Npc * c;
                Scj = Scj + Npc' * comp * Npj * c;
                Sjc = Sjc + Npj' * comp * Npc * c;
                Sjj = Sjj + Npj' * comp * Npj * c;
                
                % Compute the gravity forces
                if (this.gravityOn)
                    [feu,fep] = this.addGravityForces(feu,fep,Npc,Bp,kh,pIP,c);
                end
            end

        end

        %------------------------------------------------------------------
        % Computes the stiffness matrix and force vector contributions
        % from discontinuities in the enriched element
        function [fid, Kdd, Qad, Hdd, Sdd, Lcc, Lcj, Lcd, Ljc, Ljj, Ljd, Ldc, Ldj, Ldd, Tdc, Tdj] = getDiscontinuitiesData(this)

            % Number of displacement and pressure jump dofs
            nDiscontinuities = this.getNumberOfDiscontinuities();
            ndofa_discontinuity = this.getNumberOfDisplacementDofPerDiscontinuity();
            ndofa = nDiscontinuities * ndofa_discontinuity;
            ndofj = nDiscontinuities;
            ndofd = nDiscontinuities * 2;

            % Initialize the output data
            fid = zeros(ndofa, 1);
            Kdd = zeros(ndofa, ndofa);
            Qad = zeros(ndofa, ndofd);
            Hdd = zeros(ndofd, ndofd);
            Sdd = zeros(ndofd, ndofd);

            % Flow coupling matrices
            Lcc = zeros(this.nglp  , this.nglp);
            Lcj = zeros(this.nglp  , ndofj);
            Lcd = zeros(this.nglp  , ndofd);
            Ljc = zeros(ndofj , this.nglp);
            Ljj = zeros(ndofj , ndofj);
            Ljd = zeros(ndofj , ndofd);
            Ldc = zeros(ndofd  , this.nglp);
            Ldj = zeros(ndofd  , ndofj);
            Ldd = zeros(ndofd  , ndofd);

            % Condensation matrices
            Tdc = zeros(ndofd,this.nglp);
            Tdj = zeros(ndofd,ndofj);
            
            % Get dof vector
            a = this.getDisplacementJumpDofs(ndofa);

            % Loop through the discontinuities
            for i = 1:nDiscontinuities

                % Displacement jump dofs associated with this discontinuity
                % segment (local numbering)
                dofs_a =  ndofa_discontinuity*(i-1)+1 : ndofa_discontinuity*i;

                % Internal pressure dofs associated with this discontinuity
                % segment (local numbering)
                dofs_d = 2*(i-1)+1 : 2*i;

                % Get the discontinuity data
                [fidi, Kddi, Qadi, Hddi, Sddi, Lcci, Lcji, Lcdi, Ljci, Ljji, Ljdi, Ldci, Ldji, Lddi,~,~,~] = this.discontinuity(i).elementData(a(dofs_a),this,i);

                % Assemble the contribution of this discontinuity
                fid(dofs_a) = fid(dofs_d) + fidi;
                Kdd(dofs_a, dofs_a) = Kdd(dofs_d, dofs_d) + Kddi;
                Qad(dofs_a, dofs_d) = Qad(dofs_a, dofs_d) + Qadi;
                Hdd(dofs_d, dofs_d) = Hdd(dofs_d, dofs_d) + Hddi;
                Sdd(dofs_d, dofs_d) = Sdd(dofs_d, dofs_d) + Sddi;
                Lcc = Lcc + Lcci;
                Lcj(:,i) = Lcj(:,i) + Lcji;
                Lcd(:,dofs_d) = Lcd(:,dofs_d) + Lcdi;
                Ljc(i,:) = Ljc(i,:) + Ljci;
                Ljj(i,i) = Ljj(i,i) + Ljji;
                Ljd(i,dofs_d) = Ljd(i,dofs_d) + Ljdi;
                Ldc(dofs_d,:) = Ldc(dofs_d,:) + Ldci;
                Ldj(dofs_d,i) = Ldj(dofs_d,i) + Ldji;
                Ldd(dofs_d,dofs_d) = Ldd(dofs_d,dofs_d) + Lddi;

                % Loop through the nodes of the discontinuity to fill the
                % condensation matrix
                for j = 1:2
                    X = this.discontinuity(i).node(j,:);
                    Xn = this.shape.coordCartesianToNatural(this.node,X);
                    Np = this.shape.shapeFncMtrx(Xn);
                    Tdc(dofs_d(j),:) = Np;
                    Nenr_i = this.enrichedShapeFncValues(i, Np, X);
                    Tdj(dofs_d(j),i) = 0.5 * (Nenr_i(2) + Nenr_i(3));
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
        function n = getNumberOfDisplacementDofPerDiscontinuity(this)
            n = 2;  
            if this.addStretchingMode
                n = n + 1;
            end
            if this.addRelRotationMode
                n = n + 1;
            end
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
        % Function to get the displacement jump dofs of this element
        function a = getDisplacementJumpDofs(this,ndofa)
            i = this.nglu + this.nglp + 1;
            j = this.nglu + this.nglp + ndofa;
            a = this.ue(i:j);
        end

        %------------------------------------------------------------------
        % Function to get the displacement jump dofs of this element
        function pj = getPressureJumpDofs(this,ndofa)
            i = this.nglu + this.nglp + ndofa + 1;
            pj = this.ue(i:end);
        end

        %------------------------------------------------------------------
        % Adds the discontinuities dofs to the element dof vector
        function addEnrichmentToDofVector(this)
            nDiscontinuities = this.getNumberOfDiscontinuities();
            for i = 1:nDiscontinuities
                this.gle = [this.gle, this.discontinuity(i).dof];
            end
            this.ngle = length(this.gle);
        end

        %------------------------------------------------------------------
        % Compute the enrichment shape function matrix
        function Nenr = enrichedShapeFncMtrx(this, N, Xcar) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            Nenr = zeros(1,nDiscontinuities);  % Each discontinuity has a pressure jump dof
            for i = 1:nDiscontinuities  
                Nenr_i = this.enrichedShapeFncValues(i, N, Xcar);
                Nenr(i) = Nenr_i(1);
            end
        end

        %------------------------------------------------------------------
        % Compute the enrichment shape function values of a discontinuity
        function Nenr_i = enrichedShapeFncValues(this, id, N, Xcar) 
            phi = 0.0;
            h = this.discontinuity(id).heaviside(Xcar);
            for j = 1:this.nnd_el
                Xj = this.node(j,:);
                hj = this.discontinuity(id).heaviside(Xj);
                if (hj > 0.0)
                    phi = phi + N(:, j);
                end
            end
            Nenr_i = [h-phi;
                       -phi;        % Nbot
                      1.0-phi];     % Ntop
        end

        %------------------------------------------------------------------
        % Compute the gradient enrichment matrix
        function G = Gmatrix(this, Bu) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            G = zeros(2,nDiscontinuities);  % Each discontinuity has a pressure jump dof
            for i = 1:nDiscontinuities    
                for j = 1:this.nnd_el
                    Xj = this.node(j,:);
                    h = this.discontinuity(i).heaviside(Xj);
                    if (h > 0.0)
                        G(:,i) = G(:,i) - Bu(:, j);
                    end
                end
            end
        end

        %------------------------------------------------------------------
        % Computed the kinematic enrichment matrix for an enriched finite
        % element.
        % It calculated the enrichment matrix by considering the
        % contributions of discontinuities in the element. The enrichment
        % includes translation, stretching and relative rotation modes
        % depending on the configuration
        function Gc = kinematicEnrichment(this, Bu) 
            nDiscontinuities  = this.getNumberOfDiscontinuities();
            nDofDiscontinuity = this.getNumberOfDofPerDiscontinuity();
            Gc = zeros(4,nDofDiscontinuity * nDiscontinuities);
            for i = 1:nDiscontinuities    
                Gci = zeros(4,nDofDiscontinuity);
                % Get the discontinuity orientation vectors
                m = this.discontinuity(i).tangentialVector();
                n = this.discontinuity(i).normalVector();
                % Get the discontinuity reference point
                Xr = this.discontinuity(i).referencePoint();
                for j = 1:this.nnd_el
                    Xj = this.node(j,:);
                    h = this.discontinuity(i).heaviside(Xj);
                    if (h > 0.0)
                        % Columns of the B-matrix associated with this node
                        Buj = Bu(:, 2*(j-1) + 1 : 2*j);
                        % Add translation modes
                        Gci(:,1) = Gci(:,1) - Buj * m;
                        Gci(:,2) = Gci(:,2) - Buj * n;
                        % Add stretching mode
                        c = 3;
                        if this.addStretchingMode
                            Gci(:,c) = Gci(:,c) - Buj * (m * m') * (Xj' - Xr');
                            c = c + 1;
                        end
                        % Add relative rotation mode
                        if this.addRelRotationMode
                            mn = (n * m') - (m * n');
                            Gci(:,c) = Gci(:,c) - Buj * mn * (Xj' - Xr');
                        end
                    end
                end
                % Assemble the matrix associated with discontinuity i
                cols = nDofDiscontinuity*(i-1)+1 : nDofDiscontinuity*i;
                Gc(:,cols) = Gci;
            end
        end
        %------------------------------------------------------------------
        
    end
end