function [unvec,un,pa] = aux_pack(pa,un,me,fe)

if fe.Freesurface==1
    unvec = [un.psi;un.y;un.B];
   % if isequal(fe.Fix,'amplitude')==1 | isequal(fe.Fix,'area')==1 
     unvec=[unvec;un.Q];
    %elseif isequal(fe.Fix,'Q & amplitude')==1
    % unvec=[unvec;pa.d];
   % end
    
    if   fe.Fixwavelength==1 & fe.Fixdepth==1 & fe.Freesurface==1
     unvec=[unvec;pa.d];
    end

end







if fe.Freesurface==0
    unvec = [un.psi;un.Q];
    
    un.c = un.Q/pa.H;
    pa.alpha = pa.wavelength/pa.L;
    un.y = reshape(me.T,[],1)*pa.H/pa.d + pa.H;
    
    if fe.Embed==1
        if isequal(fe.Embedvary,'H3')==1
            unvec=[unvec;pa.H3];
        elseif isequal(fe.Embedvary,'H2')==1
            unvec=[unvec;pa.H2];   
        elseif isequal(fe.Embedvary,'H1')==1
            unvec=[unvec;pa.H1];         
        elseif isequal(fe.Embedvary,'area')==1
            unvec=[unvec;un.Area];
        elseif isequal(fe.Embedvary,'amplitude')==1
            unvec=[unvec;un.Amp];
        elseif isequal(fe.Embedvary,'s3')==1
            unvec=[unvec;pa.s3];
        elseif isequal(fe.Embedvary,'s2')==1
            unvec=[unvec;pa.s2]; 
        end
    end



end
