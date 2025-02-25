its=0;
flag=0;
while (max(abs(Hi)))>pa.TOL
%for its=1
    its=its+1;
    M=pa.M; N=pa.N;
    %%% Start timer
    tic
    
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
    %fprintf('Residuals = %s \n \n', Residuals)
       
    if (max(abs(Hi))) > 10^(7) | its>150
        'Residuals too high, divergence likely'
        flag = 1;
        break
    end
    % Compute the nonlinear contributions to the Jacobian
    [J,ja] = aux_jacobian_nonlin(J,Hi,unvec,un,pa,me,fe,de,ja,li,delta);
    
    %---------------------------------------------------------------------%
    % This code snippet is used to check our Jacobian is computed
    % correctly: the code below is very slow, and should not be used for
    % actual computations. It's a brute force linear approx of the Jacobian
    %---------------------------------------------------------------------%
    %{
    J2=J;
    for p=1:numel(unvec)
        unvec(p)=unvec(p)+delta;
        [un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);

        pa.vort_fun = eval(pa.vort_fun_str);
        pa.strat = eval(pa.strat_str);
        pa.strat_fun = eval(pa.strat_fun_str);

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
        [Hi_D,ja] = aux_eqns_vort(un,pa,me,fe,de,ja,li);
    
        J(:,p) = (Hi_D-Hi)/delta;
        unvec(p)=unvec(p)-delta;
    end
    %}
    
    % ------------------------------------------------------------------- %
    %        Newton iteration
    % ------------------------------------------------------------------- %
    
    unvec = unvec + J\transpose(-Hi);
    [un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);

   
    %{
    % Do another iteration with the same Jacobian if error is small
    if (max(abs(Hi)))<10^(-4) & (max(abs(Hi)))>pa.TOL
        for ITS=1:4
        [un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja);
        [de.dyds,de.dydss,de.dydt,de.dydtt] = aux_diff_vort(un.y,me,pa);
        [de.dpsids,de.dpsidss,de.dpsidt,de.dpsidtt] = aux_diff_vort(un.psi,me,pa);
        [Hi,ja] = aux_eqns_vort(un,pa,me,fe,de,ja);

        Residuals = max(abs(Hi));
        fprintf('Residuals = %s \n \n', Residuals)
        unvec = unvec + transpose(J\transpose(-Hi));
        'hi'
        end
    end
 
    %}

    % ------------------------------------------------------------------- %
    %  Remove nonlinear corrections to the Jacobian
    % ------------------------------------------------------------------- %
    [J] = aux_jacobian_nonlin_remove(J,pa,me,fe,ja);
    toc

    
    
end