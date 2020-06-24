function newValue = applyTransport(point, dt, D, t)
%% APPLYTRANSPORT() calculates and applies the transport operator
%
% applyTransport calculates the Laplace-Beltrami operator for diffusion
% along a curved surface using the available interpolation of the surface
% and density along it
%
% INPUTS
%   point           ===     struct: single struct object of point being evaluated, see CBPM function header for properties
%   D               ===     scalar: diffusivity of cells
%   t               ===     scalar: current simulation time
%   dt              ===     scalar: time discretisation
%
%
% Author: Solene Hegarty-Cremer
%%
%Extract necessary information
alpha = point.alpha;
beta = point.beta;
zeta = point.zeta;
advecValue = point.val;
localPointX = point.localPoint(1);

%Solve diffusion by applying 1D Laplace Beltrami operator
g = sqrt(1 + (alpha(2) + 2*alpha(3)*localPointX)^2);
dss = 2*beta(3);
tau = [1, (alpha(2) + 2*alpha(3)*localPointX)]/g;
gamma_ss = [0, 2*alpha(3)];
ds = beta(2) + 2*beta(3)*localPointX;


%Laplace Beltrami operator
lpb = g^(-2)*dss - g^(-3)*dot(tau, gamma_ss)*ds;

newValue = advecValue + D*dt*lpb;
end