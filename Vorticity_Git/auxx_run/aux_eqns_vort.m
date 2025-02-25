function [Hi,ja] = aux_eqns_vort(un,pa,me,fe,de,ja,li)

%%% unpack struct
L=pa.L; d=pa.d; M=pa.M; N=pa.N; v1=pa.v1; v2=pa.v2; MN=M*N; g=pa.g; 
Amp=un.Amp; H=pa.H;
s=me.s; t=me.t; ds=me.ds; dt=me.dt; S=me.S; T=me.T; k=me.k; Sx=me.Sx; 
Sy=me.Sy; Ex=me.Ex; Ey=me.Ey; kt=me.kt; ks=me.ks;
Y=un.Y; Psi=un.Psi; y=un.y; psi=un.psi; B=un.B; Q=un.Q; c=un.c;
dyds=de.dyds; dydss=de.dydss; dydt=de.dydt; dydtt=de.dydtt;
dpsids=de.dpsids; dpsidss=de.dpsidss; dpsidt=de.dpsidt; 
dpsidtt=de.dpsidtt;
               

% ----------------------------------------------------------------------- %
%               Equations for the streamfunction psi                      %
% ----------------------------------------------------------------------- %
%%% Boundary conditions for psi
Hi(Sy) = psi(Sy);
Hi(Ey) = psi(Ey)-Q;

%%% Governing equation for psi
jac=zeros(M*N,1);
jac(kt) =  (dyds(kt).^2+dydt(kt).^2);
if fe.Stratified == 0
    Hi(kt) = dpsidss(kt) + dpsidtt(kt) ...
          + jac(kt).*pa.vort_fun(psi(kt),Q);
else
    c=un.c;
    Hi(kt) = dpsidss(kt) + dpsidtt(kt) ...
          +  jac(kt).*g.*pa.strat_fun(psi(kt),c).*(y(kt)-psi(kt)/c);
end
% ----------------------------------------------------------------------- %
%                     Equations for the mapping Y                         %
% ----------------------------------------------------------------------- %
if fe.Freesurface==1 | fe.Bathymetry==1
    %%% Boundary conditions for Y
    if fe.Bathymetry==0
        Hi(MN+Sy) = y(Sy);
    else 
        Hi(MN+Sy) = y(Sy) - pa.bath_fun(me.s);
    end
    %%% Top wall condiion
    if fe.Bathymetry==1
        Hi(MN+Ey) = y(Ey)-pa.H;
    end
    if fe.Freesurface==1
        if fe.Ghostpoints==1
            %%% Bernoulli condition, where we use dy/ds(Sx)=dy/ds(Ex) = 0
            [Hi,yGP,dydt,jac,psiGP,dpsidt] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
        else 
            dyds(Ey(1)) = 0; dyds(Ey(end)) = 0;
            jac(Ey) = (dyds(Ey).^(2) + dydt(Ey).^(2));
            Hi(MN+Ey) = .5*(dpsidt(Ey).^2)./jac(Ey) + g*y(Ey) - B;
        end
    end
    %%% Governing equation for y
    Hi(MN+kt) = dydss(kt) + dydtt(kt);
end

% ----------------------------------------------------------------------- %
%                         Parameter equations                             %
% ----------------------------------------------------------------------- %
[Hi,ja] = aux_parameter_equation_vort(Hi,un,pa,me,fe,de,ja,li);




if fe.Solitary==1
    Hi(Ex(2:end-1)) = un.psi(Ex(2:end-1))/un.c - un.y(Ex(2:end-1));
end
