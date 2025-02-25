 function [Hi,ja] = aux_parameter_equation_vort(Hi,un,pa,me,fe,de,ja,li)
%%% Unpack struct
M=pa.M; N=pa.N; y=un.y; c=un.c; psi=un.psi; Fix = fe.Fix; 
Freesurface = fe.Freesurface; Stratified = fe.Stratified; 
dydt = de.dydt; dyds=de.dyds; s=me.s; t=me.t; Ey=me.Ey; ds=me.ds; dt=me.dt;


%%% Parameter equation
if isequal(fe.Fix,'amplitude')==1 | isequal(fe.Fix,'Q & amplitude')==1
if fe.Freesurface==1
Hi = [Hi,y(fe.Amploc)-y(fe.Amploc-M+1)+un.Amp];   
ja.parameqnloc = numel(Hi);
else
Hi = [Hi,(-un.psi(fe.Amploc-M+1)/un.c+un.y(fe.Amploc-pa.M+1)) - (-un.psi(fe.Amploc)/un.c+un.y(fe.Amploc))  - un.Amp];   
ja.parameqnloc = numel(Hi);
end
end


if isequal(fe.Fix,'Q')
Hi=[Hi,un.Q-pa.fixQ];
ja.parameqnloc = numel(Hi);
end

if isequal(fe.Fix,'area')==1
if fe.Freesurface==1
Hi = [Hi,trapz(me.s,y(Ey)-y(Ey(end))) - un.Area];   
ja.parameqnloc = numel(Hi);
else
Hi = [Hi,sum((psi/c-y).^2)*(ds*dt)*pa.alpha^2 - un.Area];   
ja.parameqnloc = numel(Hi);
end
end 

if isequal(fe.Fix,'area2')==1
eqnloc = fe.Embedloc-M+1:fe.Embedloc;
Hi = [Hi,sum((psi(eqnloc)/c-y(eqnloc)).^2)*(ds*dt)*pa.alpha^2 - un.Area2];   
ja.parameqnloc = numel(Hi);
end

if isequal(fe.Fix,'cont_H3_amp')==1
p1_old = li.H3_list(end);
p2_old = li.amp_list(end);

p1 = pa.H3;
p2 = -un.psi(fe.Amploc-M+1)/un.c+un.y(fe.Amploc-pa.M+1) - (-un.psi(fe.Amploc)/un.c+un.y(fe.Amploc));
Hi = [Hi, ((p1-p1_old)^2 + (p2-p2_old)^2 - pa.contparam)];
ja.parameqnloc = numel(Hi);
end

if isequal(fe.Fix,'cont_H3_Q')==1
p1_old = li.H3_list(end);
p2_old = li.Q_list(end);

p1 = pa.H3;
p2 = un.Q;
Hi = [Hi, (p1-p1_old)^2 + 10000*(p2-p2_old)^2 - pa.cont_ds^2];
ja.parameqnloc = numel(Hi);
end

if isequal(fe.Fix,'cont_Q_amp')
    p1_old = li.amp_list(end);
    p2_old = li.Q_list(end);
    
    p1 = un.y(fe.Amploc-M+1)-un.y(fe.Amploc);
    p2 = un.Q;
    Hi = [Hi, (p1-p1_old)^2 + (p2-p2_old)^2 - pa.cont_ds^2];
    ja.parameqnloc = numel(Hi);
end

%%% If we have a free-surface. we need to fix some length scales of the
%%% mapping
if fe.Freesurface==1
    
    if fe.Fixwavelength==1
        % Fix wavelength
        dxds = dydt;
        dxdt = -dyds;
        if fe.Mapping==0
            Hi = [Hi,trapz(s,dxds(1:M)) - pa.wavelength/2];
        else
            dfds = me.dfds(me.a,pa.L);
            Hi = [Hi,trapz(me.a,dxds(1:M)./dfds')  - pa.wavelength/2];
        end
        ja.wavelengtheqnloc = numel(Hi);
    end
    if fe.Fixdepth==1
        % Fix Depth
        Hi=[Hi,y(N*M)-pa.H];
        ja.deptheqnloc = numel(Hi);
    end

end

%if fe.Stratified==1
%    Hi = [Hi, un.c * y((N-1)*M+1) - psi((N-1)*M+1))];

% Embedded solitary waves
if fe.Embed==1
 %   E = fe.Embedloc 
 % %  Hiembed = (-1 * psi(Ex-3) + 4*psi(Ex-2) - 5*psi(Ex-1) + psi*r(Ex))/ds^2;
   %   Hi = [Hi,psi(fe.Embedloc)-psi(fe.Embedloc-1)];
    Hi=[Hi,de.dpsidss(fe.Embedloc)-pa.K_f];%psi(fe.Embedloc)-psi(fe.Embedloc-1)];
    ja.embedeqnloc = numel(Hi);
end

