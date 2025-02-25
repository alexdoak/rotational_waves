

function [unvec,me,de,fe,pa,un,li,ja,psi]= aux_listloader(me,fe,pa,un,li,ja,CNT)

delta=10^(-11);
%%% Read saved data from list
unvec = li.a_list{CNT};
qint = li.qint_list(CNT);
Amp = li.amp_list(CNT); un.Amp=Amp;
M = li.M_list(CNT); pa.M=M;
N = li.N_list(CNT); pa.N=N;
Q = li.Q_list(CNT); un.Q=Q;
L = li.L_list(CNT); pa.L=L;
d = li.d_list(CNT); pa.d=d;


if isfield(li,'v1_list')==1
    v1=li.v1_list(CNT); pa.v1=v1;
end
if isfield(li,'v2_list')==1
    v2=li.v2_list(CNT); pa.v2=v2;
end
if isfield(li,'v3_list')==1
    v3=li.v3_list(CNT); pa.v3=v3;
end
if isfield(li,'v4_list')==1
    v4=li.v4_list(CNT); pa.v4=v4;
end


Area = li.area_list(CNT); un.Area = Area;
wavelength = li.wavelength_list(CNT); pa.wavelength=wavelength;
%fe.Amploc = pa.M*pa.N;

pa.vort_fun = eval(pa.vort_fun_str);
pa.strat = eval(pa.strat_str);
pa.strat_fun = eval(pa.strat_fun_str);


[un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);
%%% Compute everything else
[me] = aux_mesh_vort(pa,me,fe);
[unvec,un,pa] = aux_pack(pa,un,me,fe);

%%% Linear Jacobian and vorticity function
[J] = aux_jacobian_lin(unvec,pa,me,fe);

% Pre-allocate
Hi = 1;

% Solve the flow
M=pa.M; N=pa.N;
%%% Start timer
tic


pa.vort_fun = eval(pa.vort_fun_str);
pa.strat = eval(pa.strat_str);
pa.strat_fun = eval(pa.strat_fun_str);

% Pre-allocate
Hi = 1;
%%% If we are varying depth and wavelength, we need to allow the mesh
%%% to vary. Hence, we must re-mesh, and recompute linear Jacobian
%%% Only matter for case of free-surface.
if fe.Fixdepth==1 & fe.Fixwavelength==1 & fe.Freesurface==1 | isequal(fe.Fix,'Q & amplitude')==1
    [me] = aux_mesh_vort(pa,me,fe);
    [J] = aux_jacobian_lin(unvec,pa,me,fe);
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
[Hi,ja] = aux_eqns_vort(un,pa,me,fe,de,ja,li);
Residuals = max(abs(Hi))

%%% Post-process: compute some flow properties
%[pa,un,de,li] = aux_postprocess(pa,me,un,de,fe,li);


if isfield(li,'maxy_list')==0
p=1;
li.maxy_list(p) =  max(li.a_list{p}(li.M_list(p)*li.N_list(p)+1:2*li.M_list(p)*li.N_list(p)));
for p=2:numel(li.Q_list)
    li.maxy_list(p) = max(li.a_list{p}(li.M_list(p)*li.N_list(p)+1:2*li.M_list(p)*li.N_list(p)));
end
end


[li] = aux_list_truncate(li,CNT);
[pa,un,de,li] = aux_postprocess(pa,me,un,de,fe,li);

psi=un.psi;
