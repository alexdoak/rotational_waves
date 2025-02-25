
function [de,Hi,Residuals] = aux_solve_single(unvec,un,pa,me,fe,de,ja,li)
M=pa.M; N=pa.N;


%%% update unknowns following a Newton iteration
[un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);
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
Residuals = max(abs(Hi)); 
fprintf('Residuals = %s \n \n', Residuals)

end