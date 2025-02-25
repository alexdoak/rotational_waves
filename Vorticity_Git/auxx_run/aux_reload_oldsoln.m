




unvec = li.a_list{end};
%H3 = li.H3_list(end);
%H2 = li.H2_list(end);

[un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);
pa.strat = eval(pa.strat_str);
pa.strat_fun = eval(pa.strat_fun_str);


y=un.y; psi=un.psi; M=pa.M; N=pa.N; c=un.c; ds=me.ds; dt=me.dt;

%%% Parameters
if fe.Freesurface==1
un.Amp=-(y(fe.Amploc)-y(fe.Amploc-M+1));   
else
un.Amp=(-un.psi(fe.Amploc-M+1)/un.c+un.y(fe.Amploc-pa.M+1)) - (-un.psi(fe.Amploc)/un.c+un.y(fe.Amploc)) ;  
end
Amp=un.Amp;

if fe.Freesurface==1
un.Area=trapz(me.s,y(Ey)-y(Ey(end))); 
else
un.Area=sum((psi/c-y).^2)*(ds*dt)*pa.alpha^2   ;
end
Area=un.Area;

un.Area2=sum((psi(eqnloc)/c-y(eqnloc)).^2)*(ds*dt)*pa.alpha^2; Area2=un.Area2;



%dp=dp/2;
%if 

%end