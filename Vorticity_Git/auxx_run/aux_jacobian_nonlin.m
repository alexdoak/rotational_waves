%********************************************************************%
%                Jacobian Computation 
% The method here is to use the pre-existing J as a base, applying
% modyifications only to terms which change every iteration
%********************************************************************% 
function [J,ja] = aux_jacobian_nonlin(J,Hi,unvec,un,pa,me,fe,de,ja,li,delta)

% Unpack structs
L=pa.L; d=pa.d; M=pa.M; N=pa.N; v1=pa.v1; v2=pa.v2; MN=M*N; g=pa.g; 
Amp=un.Amp;
s=me.s; t=me.t; ds=me.ds; dt=me.dt; S=me.S; T=me.T; k=me.k; Sx=me.Sx; 
Sy=me.Sy; Ex=me.Ex; Ey=me.Ey; kt=me.kt; ks=me.ks;
Y=un.Y; Psi=un.Psi; y=un.y; psi=un.psi; B=un.B; Q=un.Q;
dyds=de.dyds; dydss=de.dydss; dydt=de.dydt; dydtt=de.dydtt;
dpsids=de.dpsids; dpsidss=de.dpsidss; dpsidt=de.dpsidt; 
dpsidtt=de.dpsidtt;

Ds=me.Ds; Dss=me.Dss; Dt=me.Dt; Dtt=me.Dtt;

if fe.Stratified==1
    c=un.c;
end

if fe.Mapping==1

F=zeros(M*N,1);
F_delta = zeros(M*N,1);
if fe.Stratified == 0
    F(kt) =  pa.vort_fun(psi(kt),Q);
    F_delta(kt) = pa.vort_fun(psi(kt)+delta,Q) - F(kt);
else
    F(kt) =  g*pa.strat_fun(psi(kt),c).*(y(kt)-psi(kt)/c);
    F_delta(kt) = g*pa.strat_fun((psi(kt)+delta),c).*(y(kt)-(psi(kt)+delta)/c) - F(kt);
end

FE_psi = (dyds(kt).^2+dydt(kt).^2).*(F_delta(kt))/ delta;
J(kt,kt) =  J(kt,kt) + spdiags(FE_psi,0,numel(kt),numel(kt));
ja.FE_psi=FE_psi;


if fe.Stratified == 1 & fe.Freesurface == 1
    FE_y =  (dyds(kt).^2+dydt(kt).^2).*g.*pa.strat_fun(psi(kt),c);
    J(kt,M*N+kt) =  J(kt,M*N+kt) + spdiags(FE_y,0,numel(kt),numel(kt));

    ja.FE_y=FE_y;
end

%  dy/ds = D_{i,i+1}y_{i+1} + D_{i,i-1}y_{i-1}
%  (dy/ds)^2 = D_{i,i+1}^2 y_{i+1}^2 + 2D_{i,i-1}D_{i,i+1}y_{i-1}y_{i+1}
%            + D_{i,i-1}^2 y_{i-1}^2 
%%% Note that d/ds=0 on boundaries. So these modifications hit k for
%%% (d/ds)^2 and kt for (dy/dt)^2

if fe.Freesurface==1
FE_1 =  ((2*diag(Ds(k,k+1)).^2 .* y(k+1)) + (2*diag(Ds(k,k+1)).*diag(Ds(k,k-1)) .* y(k-1)))  .* F(k);
FE_m1 =  ((2*diag(Ds(k,k-1)).^2 .* y(k-1)) + (2*diag(Ds(k,k+1)).*diag(Ds(k,k-1)) .* y(k+1)))  .* F(k);
FE_M = ((2*diag(Dt(kt,kt+M)).^2 .* y(kt+M)) + (2*diag(Dt(kt,kt+M)).*diag(Dt(kt,kt-M)) .* y(kt-M)))  .* F(kt);
FE_mM = ((2*diag(Dt(kt,kt-M)).^2 .* y(kt-M)) + (2*diag(Dt(kt,kt+M)).*diag(Dt(kt,kt-M)) .* y(kt+M)))  .* F(kt);

J(k,MN+k+1) = J(k,MN+k+1) + spdiags(FE_1,0,numel(k),numel(k));
J(k,MN+k-1) = J(k,MN+k-1) + spdiags(FE_m1,0,numel(k),numel(k));
J(kt,MN+kt+M) = J(kt,MN+kt+M) + spdiags(FE_M,0,numel(kt),numel(kt));
J(kt,MN+kt-M) = J(kt,MN+kt-M) + spdiags(FE_mM,0,numel(kt),numel(kt));

ja.FE_1=FE_1;ja.FE_m1=FE_m1;ja.FE_M=FE_M;ja.FE_mM=FE_mM;
end



    
%%%%%%%%%%%%%%%%%%%%%
% BERNouli equation:%
%%%%%%%%%%%%%%%%%%%%%
%  dy/ds = D_{i,i+1}y_{i+1} + D_{i,i-1}y_{i-1}
%  (dy/ds)^2 = D_{i,i+1}^2 y_{i+1}^2 + 2D_{i,i-1}D_{i,i+1}y_{i-1}y_{i+1}
%            + D_{i,i-1}^2 y_{i-1}^2 
%  dy/dt =  D_{i,i}y_{i} + D_{i,i-M}y_{i-M} + D_{i,i-2M}y_{i-2M)
% (dy/dt)^2 = D_{i,i}^2y_{i}^2 + D_{i,i-M}^2y_{i-M}^2 
%             + D_{i,i-2M}^2y_{i-2M)^2 + 2*crossterms



%  dy/ds = (y_{i+1}-y_{i-1})/(2ds)
% (dy/ds)^2 = (y_{i+1}^2 + y_{i-1}^2 - 2y_{i+1}y_{i-1}}/(4ds^2)
%  dy/dt =  (3y_{i} - 4y_{i-M} + y_{i-2M)) /(2dt)
% (dy/dt)^2 = (9y{i}^2 + 16y{i-M}^2 + y{i-2M}^2 ...
%              -24 y{i}y_{i-M} + 6 y{i}y{i-2M} - 8y_{i-M}y_{i-2M})/(4dt^2)
% Clearly, this calculation is rather tedious. We instead modify the
% deriavites and approximate Jacobian elemnts at first order.
%0.5*(dpsidt(Ey).^2)./(dyds(Ey).^(2) + dydt(Ey).^(2)) + g*y(Ey) - B;


Hi_D=Hi;
dpsidt_d=zeros(M*N,1);
dyds_store=dyds;
dydss_store=dydss;
dte = t(end)-t(end-1);
jac = zeros(MN,1);
jac(kt) = dyds(kt).^2+dydt(kt).^2;
if fe.Freesurface==1 
    if fe.Ghostpoints==1
        %%% First compute terms we need to modify
        [Hi,yGP,dydt,jac,psiGP,dpsidt] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
        %%%%%%% Variations in psi
        %%% Ey
        psi(Ey)=psi(Ey)+delta;
        [Hi_D] =  aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
        BERN_p0 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_p0=BERN_p0'; 
        psi(Ey)=psi(Ey)-delta;
        
        %%% Ey-M
        psi(Ey-M)=psi(Ey-M)+delta;
        [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
        BERN_pmM = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_pmM=BERN_pmM';
        psi(Ey-M)=psi(Ey-M)-delta;
        
        %%% Ey-2M
        BERN_pm2M=zeros(M,1);
        
        
        %%%%%%% Variations in y
        %%% MN + Ey
        y(Ey) = y(Ey)+delta;
        dyds(Ey) = dyds(Ey) + delta*diag(Ds(Ey,Ey));
        dydss(Ey) = dydss(Ey) + delta*diag(Dss(Ey,Ey));
        [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
        BERN_y0 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_y0=BERN_y0';
        y(Ey) = y(Ey)-delta;
        dyds = dyds_store;
        dydss = dydss_store;
        
        %%% MN + Ey-M
        y(Ey-M) = y(Ey-M)+delta;
        [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
        BERN_ymM = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_ymM=BERN_ymM';
        y(Ey-M) = y(Ey-M)-delta;
        dyds = dyds_store;
        dydss = dydss_store;
            
        %%% MN + Ey-2M
        BERN_ym2M=zeros(M,1);
        
        % MN + 1
        Ey2 = Ey(1:M-1);
        dyds(Ey2) = dyds(Ey2) + delta*diag(Ds(Ey2,Ey2+1));
        dydss(Ey2) = dydss(Ey2) + delta*diag(Dss(Ey2,Ey2+1));
        [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
        BERN_y1 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_y1=BERN_y1';
        dyds = dyds_store;
        dydss = dydss_store;
        
        % MN - 1
        Ey2 = Ey(2:M);
        dyds(Ey2) = dyds(Ey2) + delta*diag(Ds(Ey2,Ey2-1));
        dydss(Ey2) = dydss(Ey2) + delta*diag(Dss(Ey2,Ey2-1));
        [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
        BERN_ym1 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_ym1=BERN_ym1';
        dyds = dyds_store;
        dydss = dydss_store;
        

    else
        Ey2 = Ey(2:end-1);
        %%%%%%% Variations in psi
        % Ey
        dpsidt_d(Ey) = dpsidt(Ey) + delta*Dt(Ey(1),Ey(1));
        Hi_D(MN+Ey2) = .5*(dpsidt_d(Ey2).^2)./(dyds(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt_d(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt_d(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_p0 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; 
        BERN_p0=BERN_p0';
        
        % Ey-M
        dpsidt_d(Ey) = dpsidt(Ey) + delta*Dt(Ey(1),Ey(1)-M);
        Hi_D(MN+Ey2) = .5*(dpsidt_d(Ey2).^2)./(dyds(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt_d(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt_d(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_pmM = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_pmM=BERN_pmM';
        
        % Ey-2M
        dpsidt_d(Ey) = dpsidt(Ey) + delta*Dt(Ey(1),Ey(1)-2*M);
        Hi_D(MN+Ey2) = .5*(dpsidt_d(Ey2).^2)./(dyds(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt_d(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt_d(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_pm2M = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_pm2M=BERN_pm2M';
        
        %%%%%%% Variations in y
        % MN + Ey
        dydt_d(Ey) = dydt(Ey) + delta*Dt(Ey(1),Ey(1));
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds(Ey2).^(2) + dydt_d(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt_d(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt_d(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_y0 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_y0=BERN_y0';
        
        % MN + Ey-M
        dydt_d(Ey) = dydt(Ey) + delta*Dt(Ey(1),Ey(1)-M);
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds(Ey2).^(2) + dydt_d(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt_d(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt_d(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_ymM = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_ymM=BERN_ymM';
        
        % MN + Ey-2M
        dydt_d(Ey) = dydt(Ey) + delta*Dt(Ey(1),Ey(1)-2*M);
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds(Ey2).^(2) + dydt_d(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt_d(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt_d(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_ym2M = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_ym2M=BERN_ym2M';
        
        % MN + 1
        dyds_d(Ey2) = dyds(Ey2) + delta*diag(Ds(Ey2,Ey2+1));
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds_d(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_y1 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_y1=BERN_y1';
        
        % MN - 1
        dyds_d(Ey2) = dyds(Ey2)  + delta*diag(Ds(Ey2,Ey2-1));
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds_d(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_ym1 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_ym1=BERN_ym1';
    end

    %%%
    J(MN+Ey,Ey) = J(MN+Ey,Ey) + spdiags(BERN_p0,0,M,M);
    J(MN+Ey,Ey-M) = J(MN+Ey,Ey-M) + spdiags(BERN_pmM,0,M,M);
    J(MN+Ey,Ey-2*M) = J(MN+Ey,Ey-2*M) + spdiags(BERN_pm2M,0,M,M);
    
    J(MN+Ey,MN+Ey) = J(MN+Ey,MN+Ey) + spdiags(BERN_y0,0,M,M);
    J(MN+Ey,MN+Ey-M) = J(MN+Ey,MN+Ey-M) + spdiags(BERN_ymM,0,M,M);
    J(MN+Ey,MN+Ey-2*M) = J(MN+Ey,MN+Ey-2*M) + spdiags(BERN_ym2M,0,M,M);
    J(MN+Ey,MN+Ey+1) = J(MN+Ey,MN+Ey+1) + spdiags(BERN_y1,0,M,M);
    J(MN+Ey,MN+Ey-1) = J(MN+Ey,MN+Ey-1) + spdiags(BERN_ym1,0,M,M);

    %%% Pack into structs
    ja.BERN_p0=BERN_p0;ja.BERN_pmM=BERN_pmM;ja.BERN_pm2M=BERN_pm2M;
    ja.BERN_y0=BERN_y0;ja.BERN_ymM=BERN_ymM;ja.BERN_ym2M=BERN_ym2M;
    ja.BERN_y1=BERN_y1;ja.BERN_ym1=BERN_ym1;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Ampltiude equation:%
%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'amplitude')==1 | isequal(fe.Fix,'Q & amplitude')==1
    if fe.Freesurface==1
    %Hi = y(M*N)-y((N-1)*M+1)-amp;   
        J(ja.parameqnloc,MN+fe.Amploc) = 1;
        J(ja.parameqnloc,MN+fe.Amploc-M+1) = -1;
    else
    %Hi = [Hi,-psi(fe.Amploc-M+1)/un.c+un.Amp];   
      J(ja.parameqnloc,fe.Amploc-M+1) = -1/un.c;
      J(ja.parameqnloc,fe.Amploc)     =  1/un.c;
    end
end
% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'area')==1
    if fe.Freesurface==1

    else
     % sum((psi/c-y).^2)/(M*N) - un.Area];   
     J(ja.parameqnloc,1:M*N) = 2/c*(psi/c-y)*(ds*dt)*pa.alpha^2 ;
    end
end
if isequal(fe.Fix,'area2')==1
     % sum((psi/c-y).^2)/(M*N) - un.Area]; 
     eqnloc = fe.Embedloc-M+1:fe.Embedloc;
     J(ja.parameqnloc,eqnloc) = 2/c*(psi(eqnloc)/c-y(eqnloc))*(ds*dt)*pa.alpha^2 ;
end

% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'Q')==1 
    if fe.Freesurface==1
        J(ja.parameqnloc,2*MN+2) = 1;
    else
        J(ja.parameqnloc,MN+1) = 1;
    end
end
% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'cont')==1
p1_old = li.H3_list(end);
p2_old = li.amp_list(end);

C = y(fe.Amploc-pa.M+1)-y(fe.Amploc)-p2_old;

J(ja.parameqnloc,fe.Amploc-M+1) =2*(un.psi(fe.Amploc-M+1)-un.psi(fe.Amploc))+2*C;
J(ja.parameqnloc,fe.Amploc) = -2*(un.psi(fe.Amploc-M+1)-un.psi(fe.Amploc))+2*C;
end
% ----------------------------------------------------------------------- %


%%%%%%%%%%%%%%%%%%%%%%%
% Wavelength equation:%
%%%%%%%%%%%%%%%%%%%%%%%
% Hi = trapz(s,dydt(1:M)) - pi;
% dydt(Sx(1)) = (-y(Sx(3))+4*y(Sx(2))-3*y(Sx(1)))/(2*dt);    
% TRAPZ hits f like ds*(.5*f1 + f2 + f3 + ... + fM-1+.5fM) 

if fe.Freesurface==1
if fe.Fixwavelength==1
    
    dfds = me.dfds(me.a,pa.L);
    Dt = me.Dt;

    WL_1 = ds * Dt(1,1) * sparse(ones(1,M))./dfds;
    WL_2 = ds * Dt(1,1+M)  * sparse(ones(1,M))./dfds; %* (y(M+1:2*M));
    WL_3 = ds * Dt(1,1+2*M) * sparse(ones(1,M))./dfds;% * (y(2*M+1:3*M));
    
    WL_1(1)=WL_1(1)/2; WL_1(end)=WL_1(end)/2;
    WL_2(1)=WL_2(1)/2; WL_2(end)=WL_2(end)/2;
    WL_3(1)=WL_3(1)/2; WL_3(end)=WL_3(end)/2;
    
    J(ja.wavelengtheqnloc,MN+1:MN+M) = J(ja.wavelengtheqnloc,MN+1:MN+M) + WL_1;
    J(ja.wavelengtheqnloc,MN+M+1:MN+2*M) = J(ja.wavelengtheqnloc,MN+M+1:MN+2*M) + WL_2;
    J(ja.wavelengtheqnloc,MN+2*M+1:MN+3*M) = J(ja.wavelengtheqnloc,MN+2*M+1:MN+3*M) + WL_3;
    
    
    ja.WL_1 = WL_1;ja.WL_2 = WL_2; ja.WL_3 = WL_3;
end

%%%%%%%%%%%%%%%%%%
% Depth equation:%
%%%%%%%%%%%%%%%%%%
if fe.Fixdepth==1
    J(ja.deptheqnloc,:)=0;
    J(ja.deptheqnloc,MN+MN) = 1;
end
end

%%% Final unknowns, recompute residuals completely due to frequency with
%%% which we change these equations
STORE=unvec;
for p=ja.extra
    unvec=STORE;
    unvec(p)=unvec(p)+delta;
    [un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);
   

    if fe.Fixdepth==1 & fe.Fixwavelength==1 |  isequal(fe.Fix,'Q & amplitude')==1
        [me] = aux_mesh_vort(pa,me,fe);
    end
    
    pa.vort_fun = eval(pa.vort_fun_str);

    %%% Differentiate the variables with second-order finite differences
    [de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = aux_diff_vort(un.psi,me,pa);
    if fe.Freesurface==1
        [de.dyds,de.dydss,de.dydt,de.dydtt] = aux_diff_vort(un.y,me,pa);
    else
        de.dyds = zeros(1,M*N); de.dydss = zeros(1,M*N); 
        de.dydtt = zeros(1,M*N); de.dydt = ones(1,M*N)*pa.alpha; 
    end
    
    %%% Evaluate Residuals of the system
    [Hi_D,ja] = aux_eqns_vort(un,pa,me,fe,de,ja,li);

    J(:,p) = (Hi_D-Hi)/delta;
  %  unvec(p)=unvec(p)-delta;
end
unvec=STORE;




%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% 
elseif fe.Mapping==0

    
%%%%%%%%%%%%%%%%%%%%%%%%% 
% Field equation for psi%
%%%%%%%%%%%%%%%%%%%%%%%%%
F=zeros(M*N,1);
F_delta = zeros(M*N,1);
if fe.Stratified == 0
    F(kt) =  pa.vort_fun(psi(kt),Q);
    F_delta(kt) = pa.vort_fun(psi(kt)+delta,Q) - F(kt);
else
    F(kt) =  g*pa.strat_fun(psi(kt),c).*(y(kt)-psi(kt)/c);
    F_delta(kt) = g*pa.strat_fun((psi(kt)+delta),c).*(y(kt)-(psi(kt)+delta)/c) - F(kt);
end

FE_psi = (dyds(kt).^2+dydt(kt).^2).*(F_delta(kt))/ delta;
J(kt,kt) =  J(kt,kt) + spdiags(FE_psi,0,numel(kt),numel(kt));
ja.FE_psi=FE_psi;


if fe.Stratified == 1 & fe.Freesurface == 1
    FE_y =  (dyds(kt).^2+dydt(kt).^2).*g.*pa.strat_fun(psi(kt),c);
    J(kt,M*N+kt) =  J(kt,M*N+kt) + spdiags(FE_y,0,numel(kt),numel(kt));

    ja.FE_y=FE_y;
end

%  dy/ds = (y_{i+1}-y_{i-1})/(2ds)
% (dy/ds)^2 = (y_{i+1}^2 + y_{i-1}^2 - 2y_{i+1}y_{i-1}}/(4ds^2)
%%% Note that d/ds=0 on boundaries. So these modifications hit k for
%%% (d/ds)^2 and kt for (dy/dt)^2

if fe.Freesurface==1
FE_1 =  (2*y(k+1) - 2*y(k-1))/(4*ds^2) .* F(k);
FE_m1 = (2*y(k-1) - 2*y(k+1))/(4*ds^2) .* F(k);
FE_M =  (2*y(kt+M) - 2*y(kt-M))/(4*dt^2) .* F(kt);
FE_mM = (2*y(kt-M) - 2*y(kt+M))/(4*dt^2) .* F(kt);

J(k,MN+k+1) = J(k,MN+k+1) + spdiags(FE_1,0,numel(k),numel(k));
J(k,MN+k-1) = J(k,MN+k-1) + spdiags(FE_m1,0,numel(k),numel(k));
J(kt,MN+kt+M) = J(kt,MN+kt+M) + spdiags(FE_M,0,numel(kt),numel(kt));
J(kt,MN+kt-M) = J(kt,MN+kt-M) + spdiags(FE_mM,0,numel(kt),numel(kt));

ja.FE_1=FE_1;ja.FE_m1=FE_m1;ja.FE_M=FE_M;ja.FE_mM=FE_mM;
end



    
%%%%%%%%%%%%%%%%%%%%%
% BERNouli equation:%
%%%%%%%%%%%%%%%%%%%%%
%  dy/ds = (y_{i+1}-y_{i-1})/(2ds)
% (dy/ds)^2 = (y_{i+1}^2 + y_{i-1}^2 - 2y_{i+1}y_{i-1}}/(4ds^2)
%  dy/dt =  (3y_{i} - 4y_{i-M} + y_{i-2M)) /(2dt)
% (dy/dt)^2 = (9y{i}^2 + 16y{i-M}^2 + y{i-2M}^2 ...
%              -24 y{i}y_{i-M} + 6 y{i}y{i-2M} - 8y_{i-M}y_{i-2M})/(4dt^2)
% Clearly, this calculation is rather tedious. We instead modify the
% deriavites and approximate Jacobian elemnts at first order.
%0.5*(dpsidt(Ey).^2)./(dyds(Ey).^(2) + dydt(Ey).^(2)) + g*y(Ey) - B;


Hi_D=Hi;
dpsidt_d=zeros(M*N,1);
dydt_d=zeros(M*N,1);
dyds_d=zeros(M*N,1);
if fe.Freesurface==1
    if fe.Ghostpoints==1
            Hi_D=Hi;
            dpsidt_d=zeros(M*N,1);
            dyds_store=dyds;
            dydss_store=dydss;
            dte = t(end)-t(end-1);
            jac = zeros(MN,1);
            jac(kt) = dyds(kt).^2+dydt(kt).^2;

      %%% First compute terms we need to modify
            [Hi,yGP,dydt,jac,psiGP,dpsidt] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
            %%%%%%% Variations in psi
            %%% Ey
            psi(Ey)=psi(Ey)+delta;
            [Hi_D] =  aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
            BERN_p0 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_p0=BERN_p0'; 
            psi(Ey)=psi(Ey)-delta;
            
            %%% Ey-M
            psi(Ey-M)=psi(Ey-M)+delta;
            [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
            BERN_pmM = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_pmM=BERN_pmM';
            psi(Ey-M)=psi(Ey-M)-delta;
            
            %%% Ey-2M
            BERN_pm2M=zeros(M,1);
            
            
            %%%%%%% Variations in y
            %%% MN + Ey
            y(Ey) = y(Ey)+delta;
            dyds(Ey) = dyds(Ey) + delta*diag(Ds(Ey,Ey));
            dydss(Ey) = dydss(Ey) + delta*diag(Dss(Ey,Ey));
            [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
            BERN_y0 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_y0=BERN_y0';
            y(Ey) = y(Ey)-delta;
            dyds = dyds_store;
            dydss = dydss_store;
            
            %%% MN + Ey-M
            y(Ey-M) = y(Ey-M)+delta;
            [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
            BERN_ymM = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_ymM=BERN_ymM';
            y(Ey-M) = y(Ey-M)-delta;
            dyds = dyds_store;
            dydss = dydss_store;
                
            %%% MN + Ey-2M
            BERN_ym2M=zeros(M,1);
            
            % MN + 1
            Ey2 = Ey(1:M-1);
            dyds(Ey2) = dyds(Ey2) + delta*diag(Ds(Ey2,Ey2+1));
            dydss(Ey2) = dydss(Ey2) + delta*diag(Dss(Ey2,Ey2+1));
            [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
            BERN_y1 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_y1=BERN_y1';
            dyds = dyds_store;
            dydss = dydss_store;
            
            % MN - 1
            Ey2 = Ey(2:M);
            dyds(Ey2) = dyds(Ey2) + delta*diag(Ds(Ey2,Ey2-1));
            dydss(Ey2) = dydss(Ey2) + delta*diag(Dss(Ey2,Ey2-1));
            [Hi_D] = aux_Berneqn(y,dyds,dydss,dydt,psi,dpsidt,dpsidss,Ey,t,M,N,Q,g,B,Hi,pa);
            BERN_ym1 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; BERN_ym1=BERN_ym1';
            dyds = dyds_store;
            dydss = dydss_store;
    
    
    
    
    else
    
    
        Ey2 = Ey(2:end-1);
        %%%%%%% Variations in psi
        % Ey
        dpsidt_d(Ey) = dpsidt(Ey) + 3*delta/(2*dt);
        Hi_D(MN+Ey2) = .5*(dpsidt_d(Ey2).^2)./(dyds(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt_d(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt_d(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_p0 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta; 
        BERN_p0=BERN_p0';
        
        % Ey-M
        dpsidt_d(Ey) = dpsidt(Ey) - 4*delta/(2*dt);
        Hi_D(MN+Ey2) = .5*(dpsidt_d(Ey2).^2)./(dyds(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt_d(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt_d(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_pmM = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_pmM=BERN_pmM';
        
        % Ey-2M
        dpsidt_d(Ey) = dpsidt(Ey) + 1*delta/(2*dt);
        Hi_D(MN+Ey2) = .5*(dpsidt_d(Ey2).^2)./(dyds(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt_d(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt_d(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_pm2M = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_pm2M=BERN_pm2M';
        
        %%%%%%% Variations in y
        % MN + Ey
        dydt_d(Ey) = dydt(Ey) + 3*delta/(2*dt);
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds(Ey2).^(2) + dydt_d(Ey2).^(2)) + g*(y(Ey2)+delta) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt_d(Ey(1)).^(2)) + g*(y(Ey(1))+delta) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt_d(Ey(end)).^(2)) + g*(y(Ey(end))+delta) - B;
        BERN_y0 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_y0=BERN_y0';
        
        % MN + Ey-M
        dydt_d(Ey) = dydt(Ey) - 4*delta/(2*dt);
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds(Ey2).^(2) + dydt_d(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt_d(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt_d(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_ymM = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_ymM=BERN_ymM';
        
        % MN + Ey-2M
        dydt_d(Ey) = dydt(Ey) + 1*delta/(2*dt);
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds(Ey2).^(2) + dydt_d(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt_d(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt_d(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_ym2M = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_ym2M=BERN_ym2M';
        
        % MN + 1
        dyds_d(Ey2) = dyds(Ey2) + delta/(2*ds); 
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds_d(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_y1 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_y1=BERN_y1';
        
        % MN - 1
        dyds_d(Ey2) = dyds(Ey2) - delta/(2*ds); 
        Hi_D(MN+Ey2) = .5*(dpsidt(Ey2).^2)./(dyds_d(Ey2).^(2) + dydt(Ey2).^(2)) + g*y(Ey2) - B;
        Hi_D(MN+Ey(1)) = .5*(dpsidt(Ey(1)).^2)./(dydt(Ey(1)).^(2)) + g*y(Ey(1)) - B;
        Hi_D(MN+Ey(end)) = .5*(dpsidt(Ey(end)).^2)./(dydt(Ey(end)).^(2)) + g*y(Ey(end)) - B;
        BERN_ym1 = (Hi_D(MN+Ey)-Hi(MN+Ey))/delta;
        BERN_ym1=BERN_ym1';
    end
    
        
    J(MN+Ey,Ey) = J(MN+Ey,Ey) + spdiags(BERN_p0,0,M,M);
    J(MN+Ey,Ey-M) = J(MN+Ey,Ey-M) + spdiags(BERN_pmM,0,M,M);
    J(MN+Ey,Ey-2*M) = J(MN+Ey,Ey-2*M) + spdiags(BERN_pm2M,0,M,M);
    
    J(MN+Ey,MN+Ey) = J(MN+Ey,MN+Ey) + spdiags(BERN_y0,0,M,M);
    J(MN+Ey,MN+Ey-M) = J(MN+Ey,MN+Ey-M) + spdiags(BERN_ymM,0,M,M);
    J(MN+Ey,MN+Ey-2*M) = J(MN+Ey,MN+Ey-2*M) + spdiags(BERN_ym2M,0,M,M);
    J(MN+Ey,MN+Ey+1) = J(MN+Ey,MN+Ey+1) + spdiags(BERN_y1,0,M,M);
    J(MN+Ey,MN+Ey-1) = J(MN+Ey,MN+Ey-1) + spdiags(BERN_ym1,0,M,M);
    
    %%% Pack into structs
    ja.BERN_p0=BERN_p0;ja.BERN_pmM=BERN_pmM;ja.BERN_pm2M=BERN_pm2M;
    ja.BERN_y0=BERN_y0;ja.BERN_ymM=BERN_ymM;ja.BERN_ym2M=BERN_ym2M;
    ja.BERN_y1=BERN_y1;ja.BERN_ym1=BERN_ym1;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Ampltiude equation:%
%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'amplitude')==1 | isequal(fe.Fix,'Q & amplitude')==1
    if fe.Freesurface==1
    %Hi = y(M*N)-y((N-1)*M+1)-amp;   
        J(ja.parameqnloc,MN+fe.Amploc) = 1;
        J(ja.parameqnloc,MN+fe.Amploc-M+1) = -1;
    else
    %Hi = [Hi,-psi(fe.Amploc-M+1)/un.c+un.Amp];   
      J(ja.parameqnloc,fe.Amploc-M+1) = -1/un.c;
      J(ja.parameqnloc,fe.Amploc)     =  1/un.c;
    end
end
% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'area')==1
    if fe.Freesurface==1

    else
     %sum((psi/c-y).^2)*(ds*dt)*pa.alpha^2 - un.Area];  
     J(ja.parameqnloc,1:M*N) = 2/c*(psi/c-y)*(ds*dt)*pa.alpha^2;
    end
end

% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'Q')==1 
    if fe.Freesurface==1
        J(ja.parameqnloc,2*MN+2) = 1;
    else
        J(ja.parameqnloc,MN+1) = 1;
    end
end
% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'cont_H3_a')==1
    psi1 = psi(fe.Amploc-M+1); psi2=psi(fe.Amploc);
    y1=y(fe.Amploc-M+1);  y2=y(fe.Amploc);
    C = y1-y2 - li.amp_list(end);
    J(ja.parameqnloc,fe.Amploc-M+1) = -2/c^2*(psi2-psi1) - 2*C/c;
    J(ja.parameqnloc,fe.Amploc) = +2/c^2*(psi2-psi1) + 2*C/c;
end
% ----------------------------------------------------------------------- %
 if isequal(fe.Fix,'cont_H3_Q')==1
    % Do nothing, let the final loop pick it up
end
% ----------------------------------------------------------------------- %
if isequal(fe.Fix,'cont_Q_amp')==1
    p1_old = li.amp_list(end);
    p2_old = li.Q_list(end);
    
    p1 = un.y(fe.Amploc-M+1)-un.y(fe.Amploc);
    p2 = un.Q;
    

    J(ja.parameqnloc,M*N+fe.Amploc-M+1) = 2*(p1-p1_old);
    J(ja.parameqnloc,M*N+fe.Amploc) = -2*(p1-p1_old);
    J(ja.parameqnloc,2*MN+2) = 2*(p2-p2_old);
end
% ----------------------------------------------------------------------- %


%%%%%%%%%%%%%%%%%%%%%%%
% Wavelength equation:%
%%%%%%%%%%%%%%%%%%%%%%%
% Hi = trapz(s,dydt(1:M)) - pi;
% dydt(Sx(1)) = (-y(Sx(3))+4*y(Sx(2))-3*y(Sx(1)))/(2*dt);    
% TRAPZ hits f like ds*(.5*f1 + f2 + f3 + ... + fM-1+.5fM) 

if fe.Freesurface==1
if fe.Fixwavelength==1
   

    WL_1 = ds * (-3/(2*dt)) * sparse(ones(1,M));
    WL_2 = ds * (4/(2*dt))  * sparse(ones(1,M)); %* (y(M+1:2*M));
    WL_3 = ds * (-1/(2*dt)) * sparse(ones(1,M));;% * (y(2*M+1:3*M));
    
    WL_1(1)=WL_1(1)/2; WL_1(end)=WL_1(end)/2;
    WL_2(1)=WL_2(1)/2; WL_2(end)=WL_2(end)/2;
    WL_3(1)=WL_3(1)/2; WL_3(end)=WL_3(end)/2;
    
    J(ja.wavelengtheqnloc,MN+1:MN+M) = J(2*MN+2,MN+1:MN+M) + WL_1;
    J(ja.wavelengtheqnloc,MN+M+1:MN+2*M) = J(2*MN+2,MN+M+1:MN+2*M) + WL_2;
    J(ja.wavelengtheqnloc,MN+2*M+1:MN+3*M) = J(2*MN+2,MN+2*M+1:MN+3*M) + WL_3;
    
    
    ja.WL_1 = WL_1;ja.WL_2 = WL_2; ja.WL_3 = WL_3;
end

%%%%%%%%%%%%%%%%%%
% Depth equation:%
%%%%%%%%%%%%%%%%%%
if fe.Fixdepth==1
    J(ja.deptheqnloc,:)=0;
    J(ja.deptheqnloc,MN+MN) = 1;
end
end

%%%%%%%%%%%%%%%%%%
% Embed equation %
%%%%%%%%%%%%%%%%%%
if fe.Embed==1
    J(ja.embedeqnloc,fe.Embedloc)=-2/ds^2;
    J(ja.embedeqnloc,fe.Embedloc-1)=2/ds^2;
%   J(ja.embedeqnloc,fe.Embedloc)=1;
 %  J(ja.embedeqnloc,fe.Embedloc-1)=-1;
end




%%%%%%%%%%%%%%%%%%%%%
% Solitary equation %
%%%%%%%%%%%%%%%%%%%%%
if fe.Solitary==1
   % Hi(Ex(2:end-1)) = un.psi(Ex(2:end-1))/un.c - un.y(Ex(2:end-1));
   %un.c = un.Q/pa.H;
    J(Ex(2:end-1),:)=0;
    J(Ex(2:end-1),Ex(2:end-1)) = speye(N-2)/un.c;

end




%%% Final unknowns, recompute residuals completely due to frequency with
%%% which we change these equations
STORE=unvec;
for p=ja.extra
    unvec=STORE;
    unvec(p)=unvec(p)+delta;
    %%% update unknowns following a Newton iteration
    [un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);

    pa.vort_fun = eval(pa.vort_fun_str);
    pa.strat = eval(pa.strat_str);
    pa.strat_fun = eval(pa.strat_fun_str);

    % Pre-allocate
    if fe.Fixdepth==1 & fe.Fixwavelength==1 & fe.Freesurface==1 | isequal(fe.Fix,'Q & amplitude')==1
        [me] = aux_mesh_vort(pa,me,fe);
    end


    %%% Differentiate the variables with second-order finite differences
    [de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = aux_diff_vort(un.psi,me,pa);
    if fe.Freesurface==1
        [de.dyds,de.dydss,de.dydt,de.dydtt] = aux_diff_vort(un.y,me,pa);
    else
        de.dyds = zeros(M*N,1); de.dydss = zeros(M*N,1); 
        de.dydtt = zeros(M*N,1); de.dydt = ones(M*N,1)*pa.alpha; 
    end

    %%% Evaluate Residuals of the system
    [Hi_D,ja] = aux_eqns_vort(un,pa,me,fe,de,ja,li);

    J(:,p) = (Hi_D-Hi)/delta;
  %  unvec(p)=unvec(p)-delta;
end
unvec=STORE;



















end
%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%% 




%}
%{
J2=J;
for p=1:numel(unvec)
    unvec=STORE;
    unvec(p)=unvec(p)+delta;
    [un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);
    
    if fe.Fixdepth==1 & fe.Fixwavelength==1 |  isequal(fe.Fix,'Q & amplitude')==1
        [me] = aux_mesh_vort(pa,me,fe);
    end

    pa.strat = eval(pa.strat_str);
    pa.strat_fun = eval(pa.strat_fun_str);

    %%% Differentiate the variables with second-order finite differences
    [de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = aux_diff_vort(un.psi,me,pa);
    if fe.Freesurface==1
        [de.dyds,de.dydss,de.dydt,de.dydtt] = aux_diff_vort(un.y,me,pa);
    else
        de.dyds = zeros(1,M*N); de.dydss = zeros(1,M*N); 
        de.dydtt = zeros(1,M*N); de.dydt = ones(1,M*N)*pa.alpha; 
    end
    
    %%% Evaluate Residuals of the system
    [Hi_D,ja] = aux_eqns_vort(un,pa,me,fe,de,ja);

    J2(:,p) = (Hi_D-Hi)/delta;
  %  unvec(p)=unvec(p)-delta;
end
%}











