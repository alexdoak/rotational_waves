function [J] = aux_jacobian_lin(unvec,pa,me,fe)
% Unpack structs
M=pa.M; N=pa.N; g=pa.g; k=me.k; Sx=me.Sx; Ex=me.Ex; Sy=me.Sy; Ey=me.Ey; 
ds=me.ds; dt=me.dt; 
Ds=me.Ds; Dss=me.Dss; Dt=me.Dt; Dtt=me.Dtt; 
kt=me.kt;

MN=M*N;
J=sparse(numel(unvec),numel(unvec));

% ----------------------------------------------------------------------- %
% This is the Jacobian for the steamfunction, in the first field equation %
% ----------------------------------------------------------------------- %
% Bottom and top boundary conditions
J(Sy,Sy) = J(Sy,Sy) + speye(M);
J(Ey,Ey) = J(Ey,Ey) + speye(M);

% Interior mesh points
J(kt,1:MN) = Dtt(kt,1:MN) + Dss(kt,1:MN);

% ----------------------------------------------------------------------- %
%         This is the Jacobian for Y, in the second field equation        %
% ----------------------------------------------------------------------- %
if fe.Freesurface==1
% Bottom and top boundary conditions
J(MN+Sy,MN+Sy) = J(MN+Sy,MN+Sy) + speye(M);
%J(MN+Ey,MN+Ey) = J(MN+Ey,MN+Ey) + g*speye(M);

% Interior mesh points
J(MN+kt,MN+(1:MN)) = Dtt(kt,:) + Dss(kt,:);
end














