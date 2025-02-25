function [me] = aux_mesh_vort(pa,me,fe)
%%% unpack struct
L=pa.L; d=pa.d; M=pa.M; N=pa.N;

%%% Build a mesh
a = linspace(0,L/2,M); da=a(2)-a(1);  
b = linspace(-d,0,N); db=b(2)-b(1); 
[A,B]=meshgrid(a,b); A=A'; B=B'; 
me.a=a; me.da=da; me.b=b; me.db=db; me.A=A; me.B=B;


%%% Index storing
k = []; kt = []; ks = []; 
for j=1:N-2
    k  = [k, j*M+2:j*M+M-1];
    kt = [kt,j*M+1:j*M+M];
end
for j=0:N-1
    ks = [ks,j*M+2:j*M+M-1];
end
k=reshape(k,[],1);
ks=reshape(ks,[],1);
kt=reshape(kt,[],1);



Sx = 1:M:M*N; Sx=reshape(Sx,[],1);
Ex = M:M:M*N; Ex=reshape(Ex,[],1);
Sy = 1:M; Sy=reshape(Sy,[],1);
Ey = (N-1)*M+1:N*M; Ey=reshape(Ey,[],1);

% Store in struct
me.k=k; me.Sx=Sx; me.Ex=Ex; me.Sy=Sy; me.Ey=Ey;


%%% Differentiation matrices
MN=M*N; Mk=numel(k); Mks=numel(ks); Mkt = numel(kt);
Ds = sparse(MN,MN);
Dss = sparse(MN,MN);
Dt = sparse(MN,MN);
Dtt = sparse(MN,MN);
Db = sparse(MN,MN);

if fe.Mapping==1

    %%% Build a mesh in mapped variables
    a = linspace(0,L/2,M); da=a(2)-a(1);  
    b = linspace(-d,0,N); db=b(2)-b(1); 
    [A,B]=meshgrid(a,b); A=A'; B=B'; 
    me.a=a; me.da=da;me.b=b; me.db=db; me.A=A; me.B=B;
    
    s = me.finv(a,L); me.s=s; 
    t = me.ginv(b,d); me.t=t; 
    [S,T]=meshgrid(s,t); S=S'; T=T'; me.S=S; me.T=T;
    
    % Junk 
    me.ds=me.da; me.dt=me.db;
    
    % Load mapping
    dfds=me.dfds(a,L); dfdss=me.dfdss(a,L); dgdt=me.dgdt(b,d); dgdtt=me.dgdtt(b,d); 
   
    % Build a vector to put into diagonals regarding the mapping s=f(alpha)
    dfdsmat = spdiags(repmat(dfds(2:end-1)',N,1),0,Mks,Mks);
    dfdssmat = spdiags(repmat(dfdss(2:end-1)',N,1),0,Mks,Mks);
 %   dgdtmat = spdiags(repmat(dgdt(2:end-1)',M,1),0,Mkt,Mkt);
 %   dgdttmat = spdiags(repmat(dgdtt(2:end-1)',M,1),0,Mkt,Mkt); repelem(dgdt(2:end-1),M)
    dgdtmat = spdiags(repelem(dgdt(2:end-1),M)',0,Mkt,Mkt);
    dgdttmat = spdiags(repelem(dgdtt(2:end-1),M)',0,Mkt,Mkt);
    
    % d/ds
    Ds(ks,ks+1) = Ds(ks,ks+1) + 1/(2*da)*dfdsmat;
    Ds(ks,ks-1) = Ds(ks,ks-1) - 1/(2*da)*dfdsmat;

    % d/dss
    Dss(ks,ks+1) = Dss(ks,ks+1) + 1/(da^2)*dfdsmat.^2 + 1/(2*da)*dfdssmat;
    Dss(ks,ks-1) = Dss(ks,ks-1) + 1/(da^2)*dfdsmat.^2 - 1/(2*da)*dfdssmat;
    Dss(ks,ks) = Dss(ks,ks) - 2/(da^2)*dfdsmat.^2;

    ds0 = me.s(2)-me.s(1);
    Dss(Sx,Sx+1) = Dss(Sx,Sx+1) +  speye(N) * 2/(ds0^2);
    Dss(Sx,Sx) = Dss(Sx,Sx) -  speye(N) * 2/(ds0^2);
    
    dse = me.s(M)-me.s(M-1);
    Dss(Ex,Ex-1) = Dss(Ex,Ex-1) +  speye(N) * 2/(dse^2);
    Dss(Ex,Ex) = Dss(Ex,Ex) -  speye(N) * 2/(dse^2);



    % d/dt
    Dt(kt,kt+M) = Dt(kt,kt+M) + 1/(2*db)*dgdtmat;
    Dt(kt,kt-M) = Dt(kt,kt-M) - 1/(2*db)*dgdtmat;
    
    [DtSycoefs] = coeffs_alpha(2,0,t(1:3)-t(1),1);
    Dt(Sy,Sy+2*M) = Dt(Sy,Sy+2*M) + speye(M) * DtSycoefs(3);
    Dt(Sy,Sy+M) =  Dt(Sy,Sy+M) + speye(M) * DtSycoefs(2);
    Dt(Sy,Sy) =  Dt(Sy,Sy) + speye(M) * DtSycoefs(1);
    
    [DtEycoefs] = coeffs_alpha(2,0,t(end-2:end)-t(end),1);
    Dt(Ey,Ey-2*M) = Dt(Ey,Ey-2*M) + speye(M) * DtEycoefs(1);
    Dt(Ey,Ey-M) =  Dt(Ey,Ey-M) + speye(M) * DtEycoefs(2);
    Dt(Ey,Ey) =  Dt(Ey,Ey) + speye(M) * DtEycoefs(3);
    
    % d/dtt
    Dtt(kt,kt+M) = Dtt(kt,kt+M) + 1/(db^2)*dgdtmat.^2 + 1/(2*db)*dgdttmat;
    Dtt(kt,kt-M) = Dtt(kt,kt-M) + 1/(db^2)*dgdtmat.^2 - 1/(2*db)*dgdttmat;
    Dtt(kt,kt) = Dtt(kt,kt) - 2/(db^2)*dgdtmat.^2;
    
    
    [DttSycoefs] = coeffs_alpha(3,0,t(1:4)-t(1),2);
    Dtt(Sy,Sy+3*M) = Dtt(Sy,Sy+3*M) + speye(M) * DttSycoefs(4);
    Dtt(Sy,Sy+2*M) = Dtt(Sy,Sy+2*M) + speye(M) * DttSycoefs(3);
    Dtt(Sy,Sy+M) =  Dtt(Sy,Sy+M) + speye(M) * DttSycoefs(2);
    Dtt(Sy,Sy) =  Dtt(Sy,Sy) + speye(M) * DttSycoefs(1);
    
    [DttEycoefs] = coeffs_alpha(3,0,t(end-3:end)-t(end),2);
    Dtt(Ey,Ey-3*M) = Dtt(Ey,Ey-3*M) + speye(M) * DttEycoefs(1);
    Dtt(Ey,Ey-2*M) = Dtt(Ey,Ey-2*M) + speye(M) * DttEycoefs(2);
    Dtt(Ey,Ey-M) =  Dtt(Ey,Ey-M) + speye(M) * DttEycoefs(3);
    Dtt(Ey,Ey) =  Dtt(Ey,Ey) + speye(M) * DttEycoefs(4);
    %}
elseif fe.Mapping==0
    % No mapping, let s=a and b=t
    s=a; t=b; ds=da; dt=db; S=A; T=B;
    me.s=s; me.t=t; me.ds=ds; me.dt=dt; me.S=S; me.T=T;

    % d/ds
    Ds(ks,ks+1) = Ds(ks,ks+1) + speye(Mks)/(2*ds);
    Ds(ks,ks-1) = Ds(ks,ks-1) - speye(Mks)/(2*ds);
    
    % d/dss
    Dss(ks,ks+1) = Dss(ks,ks+1) + speye(Mks) / (ds^2);
    Dss(ks,ks-1) = Dss(ks,ks-1) + speye(Mks) / (ds^2);
    Dss(ks,ks) = Dss(ks,ks) - speye(Mks) * 2/(ds^2);
    
    Dss(Sx,Sx+1) = Dss(Sx,Sx+1) +  speye(N) * 2/(ds^2);
    Dss(Sx,Sx) = Dss(Sx,Sx) -  speye(N) * 2/(ds^2);
    
    Dss(Ex,Ex-1) = Dss(Ex,Ex-1) +  speye(N) * 2/(ds^2);
    Dss(Ex,Ex) = Dss(Ex,Ex) -  speye(N) * 2/(ds^2);

    % d/dt
    Dt(kt,kt+M) = Dt(kt,kt+M) + speye(Mkt)/(2*dt);
    Dt(kt,kt-M) = Dt(kt,kt-M) - speye(Mkt)/(2*dt);
    
    Dt(Sy,Sy+2*M) = Dt(Sy,Sy+2*M) + speye(M) * (-1)/(2*dt);
    Dt(Sy,Sy+M) =  Dt(Sy,Sy+M) + speye(M) * (4)/(2*dt);
    Dt(Sy,Sy) =  Dt(Sy,Sy) + speye(M) * (-3)/(2*dt);
    
    Dt(Ey,Ey-2*M) = Dt(Ey,Ey-2*M) + speye(M) * (1)/(2*dt);
    Dt(Ey,Ey-M) =  Dt(Ey,Ey-M) + speye(M) * (-4)/(2*dt);
    Dt(Ey,Ey) =  Dt(Ey,Ey) + speye(M) * (3)/(2*dt);
    
    % d/dtt
    Dtt(kt,kt+M) = Dtt(kt,kt+M) + speye(Mkt)/(dt^2);
    Dtt(kt,kt-M) = Dtt(kt,kt-M) + speye(Mkt)/(dt^2);
    Dtt(kt,kt) = Dtt(kt,kt) - speye(Mkt) * 2/(dt^2);
    
    Dtt(Sy,Sy+3*M) = Dtt(Sy,Sy+3*M) + speye(M) * (-1)/(dt^2);
    Dtt(Sy,Sy+2*M) = Dtt(Sy,Sy+2*M) + speye(M) * (4)/(dt^2);
    Dtt(Sy,Sy+M) =  Dtt(Sy,Sy+M) + speye(M) * (-5)/(dt^2);
    Dtt(Sy,Sy) =  Dtt(Sy,Sy) + speye(M) * (2)/(dt^2);
    
    Dtt(Ey,Ey-3*M) = Dtt(Ey,Ey-3*M) + speye(M) * (-1)/(dt^2);
    Dtt(Ey,Ey-2*M) = Dtt(Ey,Ey-2*M) + speye(M) * (4)/(dt^2);
    Dtt(Ey,Ey-M) =  Dtt(Ey,Ey-M) + speye(M) * (-5)/(dt^2);
    Dtt(Ey,Ey) =  Dtt(Ey,Ey) + speye(M) * (2)/(dt^2);
    
end



me.Dt=Dt;
me.Dtt=Dtt;
me.Ds=Ds;
me.Dss=Dss;
me.kt=kt;
me.ks=ks;