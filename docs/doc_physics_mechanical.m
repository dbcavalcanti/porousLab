%% Mechanical physics
%
%% Theorical background
%
% *Primary field*: displacement field
%
% $\mathbf{u} = \left[u_x, u_y\right]^T$
% 
% *Governing equation*: Linear momentum balance
%
% $\nabla \cdot \mathbf{\sigma} + \rho \, \mathbf{g} = \mathbf{0}$
%
% where:
%
% * $\mathbf{\sigma}$: stress tensor
% * $\rho$: density
% * $\mathbf{g}$: gravity acceleration vector
% 
% *Weak form*: 
%
% $\int_{\Omega}{\nabla\delta\mathbf{u} :\mathbf{\sigma}\, d\Omega}-\int_{\Gamma_t}{\delta\mathbf{u} \cdot \mathbf{t}\, d\Omega} -\int_{\Omega}{\rho \,\delta \mathbf{u} \cdot\mathbf{g}\, d\Omega} = 0$
%
%% Main classes associated
%
% * <Model_M.html Model_M>: initialize the finite element
% model with the elements of the <RegularElement_M.html RegularElement_M>
% class.
% * <RegularElement_M.html RegularElement_M>: initialize the <intPoint.html intPoint> object with <Material_M.html Material_M>
