function [un,pa,ja] = aux_unpack(unvec,pa,un,fe,ja)


MN=pa.M*pa.N;
un.psi = unvec(1:MN);

if fe.Freesurface==1 
    un.y = unvec(MN+1:2*MN);
    un.B = unvec(2*MN+1);

    if isequal(fe.Fix,'Q & amplitude')==1
    pa.d=unvec(2*MN+2);
    else
    un.Q = unvec(2*MN+2);
    end
    
    if  fe.Fixwavelength==1 & fe.Fixdepth==1 &  fe.Freesurface==1
        pa.d=unvec(2*MN+3);
        ja.extra = 2*MN+1:2*MN+3;
    else
        ja.extra = 2*MN+1:2*MN+2;
    end

    un.c = un.Q/pa.H;


end

if fe.Freesurface==0
    un.Q = unvec(MN+1);
 


    ja.extra = MN+1;
    un.c = un.Q/pa.H;

    if fe.Embed==1
    if isequal(fe.Embedvary,'H3')==1
        pa.H3=unvec(MN+2);
    elseif isequal(fe.Embedvary,'H2')==1
        pa.H2=unvec(MN+2); %pa.H3=(1-pa.H2)/2-0.005;
    elseif isequal(fe.Embedvary,'H1')==1
        pa.H1=unvec(MN+2); %pa.H3=(1-pa.H2)/2-0.005;
    elseif isequal(fe.Embedvary,'area')==1
        un.Area=unvec(MN+2);
    elseif isequal(fe.Embedvary,'amplitude')==1
         un.Amp=unvec(MN+2);
    elseif isequal(fe.Embedvary,'s3')==1
        pa.s3=unvec(MN+2);
    elseif isequal(fe.Embedvary,'s2')==1
         pa.s2=unvec(MN+2); pa.s5=pa.s2;
    end
    ja.extra=MN+1:MN+2;
    end


end

