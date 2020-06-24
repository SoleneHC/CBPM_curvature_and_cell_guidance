function newNormalLocal = getNewNormal(newLocalPoint, localGridPoint, oldNormalLocal, alpha)
%% Gets new normal either by approximation or local interpolation
%
%
%
%
%
%
% Author: Solene Hegarty-Cremer
%%
norm_tol = 1e-4;


tangent = [1, 2*alpha(3)*newLocalPoint(1) + alpha(2)];
newNormalLocal = [-tangent(2); tangent(1)];
newNormalLocal = newNormalLocal/norm(newNormalLocal);


normalSign = sign(dot(oldNormalLocal, newNormalLocal));

newNormalLocal = normalSign*newNormalLocal;
end