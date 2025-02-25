clear all

PATH_TO_FOLDER = '/home/adoak/Documents/MATLAB/Vorticity_new';
addpath([PATH_TO_FOLDER '/auxx_run'])
addpath([PATH_TO_FOLDER '/auxx_load'])
addpath([PATH_TO_FOLDER '/auxx_plot'])
addpath([PATH_TO_FOLDER '/auxx_IG'])
addpath([PATH_TO_FOLDER '/auxx_disprel'])

% ----------------------------------------------------------------------- %
%  Pre-amble before solving: pick variables, create mesh, and such        %
% ----------------------------------------------------------------------- %
%%% Pre-allocate structures for mesh (me), unknowns (un), derivatives (de),
%%% features (fe), parameters (pa), lists (li), jacobian terms (ja)
me=struct; de=struct; fe=struct; pa=struct; un=struct; li=struct;ja=struct;
addpath  /home/adoak/Documents/MATLAB/Vorticity/vort/affine

%%% Load data from .mat file. Turn fe.Loaddata=0 if you do not wish to load
%%% a solution
fe.Loaddata=0;
if fe.Loaddata==1
   fe.Ghostpoints=0;
   filename = 'vort=10.mat';
   [unvec,me,de,fe,pa,un,li,ja,psi,Amp,Area,delta] = aux_load_reducedsize(filename,fe);

    
    %%% This will load data from the list, element CNT
    CNT=numel(li.Q_list);
    [unvec,me,de,fe,pa,un,li,ja,psi]= aux_listloader(me,fe,pa,un,li,ja,CNT);
    Amp=un.Amp; Q=un.Q; Area=un.Area;
end


%%% Decide on some features of the solver: what do we fix, what do we vary?
fe.Newlists=1; %New lists for storing solution branch?
if isfield(li,'a_list')==0
fe.Newlists=1;
end

fe.Fix =  'amplitude';%'Q';%'amplitude';  'cont_Q_amp'; % 'amplitude'%'amplitude';
fe.Fixdepth = 1;
fe.Fixwavelength = 1;
fe.Scont = 0;
fe.Freesurface = 1; 
fe.Stratified = 0; 
fe.Bathymetry = 0;
fe.Solitary = 0;
fe.Mapping=0;
fe.Ghostpoints = 0;
fe.Embed=0; fe.Embedloc = 1;
fe.Plotwave = 0; % Plot the wave for every solution computed?


%%% Choose your coordinate mapping. alpha=f(s), s=finv(alpha). dfds=f'
gammaf=1; 
gammag=1;
pa.gammaf=gammaf;
pa.gammag=gammag;
if fe.Mapping==1
    me.finv = @(a,L) a.^pa.gammaf * (L/2)^(1-pa.gammaf);
    me.dfds = @(a,L) 1/pa.gammaf * (L/2)^(pa.gammaf-1)*a.^(1-pa.gammaf);
    me.dfdss = @(a,L) 1/pa.gammaf*(1/pa.gammaf-1) * (L/2)^(2*pa.gammaf-2)*a.^(1-2*pa.gammaf);

    me.ginv = @(b,d) -(-b).^pa.gammag * (d)^(1-pa.gammag);
    me.dgdt = @(b,d) 1/pa.gammag*d^(pa.gammag-1)*(-b).^(1-pa.gammag);
    me.dgdtt = @(b,d) -1/pa.gammag*(1/pa.gammag-1) * (d)^(2*pa.gammag-2)*(-b).^(1-2*pa.gammag);
end

%%% Create an initial guess from linear theory
if exist("psi",'var') == 0
    pa.g=1; pa.H=1; pa.wavelength=2*pi;
    pa.v1=50; pa.v2=0;  pa.wavelength=2*pi;  un.Q=0.035; 
    pa.vort_fun = @(psi,Q) pa.v1*psi + pa.v2;
  
    [pa,un,me,fe,delta] = aux_create_IG_vort(pa,un,me,fe);
end

%%% Create mesh 
if exist("s") == 0
    [me] = aux_mesh_vort(pa,me,fe);
end


%%% Modify initial guess
%[un,pa,me,fe] = aux_IG_mappingfit(un,pa,me,fe,1,1);
%[un,pa,me,fe] = aux_vary_wavelength(un,pa,me,fe,1.5);
%[un,pa,me,fe] = aux_vary_MN_vort(un,pa,me,fe,500,500);


%%% Store unknowns in a vector 
[unvec,un,pa] = aux_pack(pa,un,me,fe);

%%% Reset data_storage
if fe.Newlists == 1
    [li] = aux_newlists(unvec,pa,un,fe);
end

% ----------------------------------------------------------------------- %
%                     Newton-Raphson solver
% ----------------------------------------------------------------------- %
[J] = aux_jacobian_lin(unvec,pa,me,fe);

%%% Update bifurcation paramaters, since can't loop through un.""
Area=un.Area;
wavelength=pa.wavelength;  d=pa.d; 
Q=un.Q;


for Amp=0.001
    pa.TOL = (10^(-11)/min([me.dt,me.ds]).^2);
    %%% Update Bifurcation parameteres    
    un.Amp=Amp;
    un.Area=Area;
    un.Q=Q;
    pa.wavelength=wavelength;
    if fe.Freesurface==0
      pa.d = pa.H*pa.L/wavelength;
    end
    pa.fixQ = un.Q;

    % Mesh and pack a
    [me] = aux_mesh_vort(pa,me,fe);
    [unvec,un,pa] = aux_pack(pa,un,me,fe);

    % Linear Jacobian and vorticity function
    [J] = aux_jacobian_lin(unvec,pa,me,fe);

    % Save vorticity, stratification, and mappings, as strings,
    % for matlab reasons
    pa.vort_fun_str = func2str(pa.vort_fun);
    pa.strat_str = func2str(pa.strat);
    pa.strat_fun_str = func2str(pa.strat_fun);
    if fe.Mapping==1
    me.finv_str = func2str(me.finv);
    me.dfds_str = func2str(me.dfds);
    me.dfdss_str = func2str(me.dfdss);
    me.ginv_str = func2str(me.ginv);
    me.dgdt_str = func2str(me.dgdt);
    me.dgdtt_str = func2str(me.dgdtt);
    end

    % Pre-allocate
    Hi = 1;
    
    % Solve the flow
    [de,Hi,Residuals] = aux_solve_single(unvec,un,pa,me,fe,de,ja,li);
    run aux_solve
    '*******************************************************************'
    
    %%% Post-process: compute some flow properties
    [pa,un,de,li] = aux_postprocess(pa,me,un,de,fe,li);
    
    %%% Plot the flow
    pa.wavelength_true = trapz(me.s,de.dydt(1:pa.M))*2;
    if fe.Plotwave == 1
        aux_plotstream_vort(un,pa,de,fe,50,[],[],0.03)
    end
    
    %%% Add variables to list
    [li] = aux_listadd(unvec,pa,me,un,de,fe,li);
    
    %%% Update bifurcation paramaters
    Amp=un.Amp;
    Area=un.Area;
    Q=un.Q; 
    wavelength=pa.wavelength;
    me=orderfields(me);
    de=orderfields(de);
    fe=orderfields(fe);
    pa=orderfields(pa);
    un=orderfields(un);
    li=orderfields(li);
    ja=orderfields(ja);


    save('safety_delete.mat')
      
% filename=['A=' num2str(Amp) '_crapper_M=' num2str(pa.M) '_N=' num2str(pa.N) '_d=' num2str(pa.d) '.mat']; 
% aux_save_reducedsize(filename,pa,un,me,fe,li,delta,1);

end


