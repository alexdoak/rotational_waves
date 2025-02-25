function [unvec,me,de,fe,pa,un,li,ja,psi,Amp,Area,delta] = aux_load_reducedsize(filename,fe)

Loaddata=fe.Loaddata;
me=struct; de=struct; pa=struct; un=struct; li=struct;ja=struct;
load(filename);
fe.Ghostpoints=0;
fe.Loaddata=Loaddata;
M=pa.M; N=pa.N;


%%% This fixes old workspaces, shifting location of some data for aecestic
%%% choice

if isfield(pa,'K_f')==0
    pa.K_f=0
end

if isfield(fe,'Mapping')==0
    fe.Mapping=0;
end
if isfield(pa,'Amp')
un.Amp=pa.Amp;  pa=rmfield(pa,'Amp'); 
end
if isfield(pa,'Area')
un.Area=pa.Area; pa=rmfield(pa,'Area');
end

if isfield(fe,'Soliary') == 0
fe.Solitary=0;
end
if isfield(fe,'Fixscale')==1
    if fe.Fixscale=='depth'
        fe.Fixdepth=1;
        fe.Fixwavelength=0;
    else
        fe.Fixdepth=0;
        fe.Fixwavelength=1;
    end
    fe=rmfield(fe,'Fixscale');
end

if isfield(pa,'vort_fun_str')==1
    pa.vort_fun = eval(pa.vort_fun_str);
    pa.strat = eval(pa.strat_str);
    pa.strat_fun = eval(pa.strat_fun_str);
end

if fe.Mapping==1
    me.finv=eval(me.finv_str);
    me.ginv=eval(me.ginv_str);
    me.dfds=eval(me.dfds_str);
    me.dfdss=eval(me.dfdss_str);
    me.dgdt=eval(me.dgdt_str);
    me.dgdtt=eval(me.dgdtt_str);
end

if isfield(fe,'amploc')==1
    fe.Amploc=fe.amploc;
    fe=rmfield(fe,'amploc');
end


%%% Reshaping, due to inconsistency in old data files
un.y=reshape(un.y,[],1);
un.psi=reshape(un.psi,[],1);


%%% Mesh, reshape matrix from vectors
[me] = aux_mesh_vort(pa,me,fe);  
un.Y = reshape(un.y,[],pa.N); un.Psi = reshape(un.psi,[],pa.N);

%%% Store unknowns in a vector unvec
[unvec,un,pa] = aux_pack(pa,un,me,fe);


%%% update unknowns following a Newton iteration
[un,ja] = aux_unpack(unvec,pa,un,fe,ja);

%%% Differentiate the variables with second-order finite differences
[de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = aux_diff_vort(un.psi,me,pa);
if fe.Freesurface==1
    [de.dyds,de.dydss,de.dydt,de.dydtt] = aux_diff_vort(un.y,me,pa);
else
    de.dyds = zeros(M*N,1); de.dydss = zeros(M*N,1); 
    de.dydtt = zeros(M*N,1); de.dydt = ones(M*N,1)*pa.alpha; 
end

%%% Evaluate Residuals of the system
[Hi,ja] = aux_eqns_vort(un,pa,me,fe,de,ja,li);
Residuals = max(abs(Hi)); 
fprintf('Residuals = %s \n \n', Residuals)

psi=un.psi; Amp=un.Amp; Area=un.Area; delta=10^(-11);



%%% Post-process: compute some flow properties
%[pa,un,de,li] = aux_postprocess(pa,me,un,de,fe,li);

%%% Plot the flow
pa.wavelength_true = trapz(me.s,de.dydt(1:pa.M))*2;
%}




