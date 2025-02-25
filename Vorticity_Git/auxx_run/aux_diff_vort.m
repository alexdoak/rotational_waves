function [dyds,dydss,dydt,dydtt] = aux_diff_vort(y,me,pa)
M=pa.M; N=pa.N;
k=me.k; ds=me.ds; dt=me.dt; Sx=me.Sx; Sy=me.Sy; Ex=me.Ex; Ey=me.Ey;

dyds = me.Ds*y;
dydss = me.Dss*y;
dydt = me.Dt*y;
dydtt = me.Dtt*y;

dyds=dyds;
dydss=dydss;
dydt=dydt;
dydtt=dydtt;

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Difference equations %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central difference equations
dyds(k)  = (y(k+1)-y(k-1))/(2*ds);
dydss(k) = (y(k+1) - 2.*y(k) +  y(k-1)) ./(ds^2);
dydt(k) =  (y(k+M) - y(k-M))/(2*dt);
dydtt(k) = (y(k+M) - 2.*y(k) +  y(k-M)) ./(dt^2);

% Differences at s=0,
dyds(Sx) = (-3*y(Sx)+4*y(Sx+1)-y(Sx+2))/(2*ds);
dydt(Sx(2:end-1)) = (y(Sx(3:end))-y(Sx(1:end-2)))/(2*dt);
dydt(Sx(1)) = (-y(Sx(3))+4*y(Sx(2))-3*y(Sx(1)))/(2*dt);
dydt(Sx(end)) = (+y(Sx(end-2))-4*y(Sx(end-1))+3*y(Sx(end)))/(2*dt);

% Differences at s=pi,
dyds(Ex) = (3*y(Ex)-4*y(Ex-1)+y(Ex-2))/(2*ds);
dydt(Ex(2:end-1)) = (y(Ex(3:end))-y(Ex(1:end-2)))/(2*dt);
dydt(Ex(1)) = (-y(Ex(3))+4*y(Ex(2))-3*y(Ex(1)))/(2*dt);
dydt(Ex(end)) = (+y(Ex(end-2))-4*y(Ex(end-1))+3*y(Ex(end)))/(2*dt);

% Differences at  t=0
dydt(Ey) =  (3*y(Ey) - 4*y(Ey-M) + y(Ey-2*M))/(2*dt); 
dyds(Ey(2:M-1)) = (y(Ey(2:M-1)+1)-y(Ey(2:M-1)-1))/(2*ds); 


% Differences at  t=-d
dydt(Sy) =  (-3*y(Sy) + 4*y(Sy+M) - y(Sy+2*M))/(2*dt); 
dyds(Sy(2:M-1)) = (y(Sy(2:M-1)+1)-y(Sy(2:M-1)-1))/(2*ds); 

% Impose dy/ds(Ey(1))=0, same for psi, to put into Bernoulli
% equation
dyds(Ey(1))=0;
dyds(Ey(end))=0;
%}

