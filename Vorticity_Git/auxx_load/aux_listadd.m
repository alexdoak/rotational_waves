function [li] = aux_listadd(a,pa,me,un,de,fe,li)

if numel(li.a_list)==0
    li.a_list{1}=a;
  
    li.qint_list(1) = pa.qint;
    li.amp_list(1) = un.Amp;%un.y(fe.Amploc)-un.y(fe.Amploc-pa.M+1);
    li.M_list(1) = pa.M;
    li.N_list(1) = pa.N;
    li.Q_list(1) = un.Q;
    li.area_list(1) = un.Area;
    li.wavelength_list(1) = pa.wavelength;
    li.L_list(1) = pa.L;
    li.d_list(1) = pa.d;
    li.speed_list(1) = -un.U(end,1);
    li.Burns_list(1) = trapz(un.Y(end,:),(1./(un.U(end,:)).^2));
    
    meandepth = trapz(un.X(:,end),un.Y(:,end))/(pa.wavelength/2);
    li.meandepth_list(1) = meandepth;
    li.Bernhead_list(1) = un.B*2; 
   
    li.K_list(1) = de.dpsidss(fe.Embedloc);

    li.min_q_list(1) = min(un.V(:,end).^2 + un.U(:,end).^2).^.5;

    li.maxy_list(1) = max(un.y);
 
    
    if isfield(pa,'v1')
        li.v1_list = pa.v1;
    end
    if isfield(pa,'v2')
        li.v2_list = pa.v2;
    end
    if isfield(pa,'v3')
        li.v3_list = pa.v3;
    end
    if isfield(pa,'v4')
        li.v4_list = pa.v4;
    end
    if isfield(pa,'core')
        li.stag_list=pa.core;
    end
    
  
else
    li.a_list{end+1,:} = a;
    
    li.qint_list(end+1) = pa.qint;
    li.amp_list(end+1) = un.Amp;%un.y(fe.Amploc)-un.y(fe.Amploc-pa.M+1);
    li.M_list(end+1) = pa.M;
    li.N_list(end+1) = pa.N;
    li.Q_list(end+1) = un.Q;
    li.area_list(end+1) = un.Area;
    li.wavelength_list(end+1) = pa.wavelength;
    li.L_list(end+1) = pa.L;
    li.d_list(end+1) = pa.d;
    li.speed_list(end+1) = -un.U(end,1);
    li.Burns_list(end+1) = trapz(un.Y(end,:),(1./(un.U(end,:)).^2));
    
    li.maxy_list(end+1) = max(un.y);

    meandepth = trapz(un.X(:,end),un.Y(:,end))/(pa.wavelength/2);
    li.meandepth_list(end+1) = meandepth;
    li.Bernhead_list(end+1) = un.B*2; 
    M=pa.M; N=pa.N; 

    li.K_list(end+1) = de.dpsidss(fe.Embedloc);
    li.min_q_list(end+1) = min(un.V(:,end).^2 + un.U(:,end).^2).^.5;

    if isfield(pa,'v1')
         li.v1_list(end+1) = pa.v1;
    end
    if isfield(pa,'v2')
         li.v2_list(end+1) = pa.v2;
    end
    if isfield(pa,'v3')
         li.v3_list(end+1) = pa.v3;
    end
    if isfield(pa,'v4')
         li.v4_list(end+1) = pa.v4;
    end
     if isfield(pa,'core')
         li.stag_list(end+1)=pa.core;
     end


    if fe.Freesurface==1
  %  for p=1:numel(li.Q_list)
  %     li.amp_list_KS(p) = li.a_list{p}(2*M*N-M+1) - li.meandepth_list(p);
  %  end
    
    elseif fe.Freesurface==0
        li.TE_list(end+1)=un.TE;
        li.KE_list(end+1)=un.KE;
        li.PE_list(end+1)=un.PE;
        li.mass_list(end+1) = un.mass;

        
    end

       

end






