function[pa,un,de,li] = aux_postprocess(pa,me,un,de,fe,li)

%%% Unpack structs
L=pa.L; d=pa.d; M=pa.M; N=pa.N; g=pa.g; wavelength=pa.wavelength;
v1=pa.v1; v2=pa.v2; Amp=un.Amp; vort_fun=pa.vort_fun;
y=un.y; psi=un.psi; 
s=me.s; t=me.t; ds=me.ds; dt=me.dt; S=me.S; T=me.T; k=me.k; Sx=me.Sx; 
Sy=me.Sy; Ex=me.Ex; Ey=me.Ey;
c=un.c;
k=me.k; kt=me.kt; ks=me.ks; Mkt=numel(kt);
db=me.db; a=me.a; b=me.b;

%%% Boundary derivatives, unpack de.
dyds=de.dyds; dydss=de.dydss; dydt=de.dydt; dydtt=de.dydtt;
dpsids=de.dpsids; dpsidss=de.dpsidss; dpsidt=de.dpsidt; 
dpsidtt=de.dpsidtt;




%%% Compute X
if fe.Mapping==0
    un.X=zeros(M,N);
    for j=1:N
        un.X(:,j) = cumtrapz(s,dydt((j-1)*M+1:j*M));
    end
else
    
    dfds=me.dfds(me.a,L); 
    un.X=zeros(M,N);
    for j=1:N
        un.X(:,j) = cumtrapz(a,dydt((j-1)*M+1:j*M)./dfds');
    end
   %  un.X2=zeros(M,N);
   % for j=1:N
   %     un.X2(:,j) = cumtrapz(s,dydt((j-1)*M+1:j*M));
   % end
end


%%% Reshape Y, Psi
un.Y = reshape(y,[],N); 
un.Psi = reshape(psi,[],N);

% Velocity field 
un.u = (dpsids.*dyds + dpsidt.*dydt)./(dyds.^(2) + dydt.^(2));
un.v = (-dpsids.*dydt + dpsidt.*dyds)./(dyds.^(2) + dydt.^(2));
un.U = reshape(un.u,[],N); 
un.V = reshape(un.v,[],N);

% Vorticity
jac =  (dyds.^2+dydt.^2);
un.vort = zeros(M,N);
un.vort = -1./jac .* (dpsidtt + dpsidss); 
%un.vort(Ey)=0;
un.vort = reshape(un.vort,[],N);


%%% Mean depth, mean 'speed', wavelengths etc.
pa.H_mean=mean(un.Y(:,end));

q = ((dpsids(Ey).^2 + dpsidt(Ey).^2)./(dyds(Ey).^(2) + dydt(Ey).^(2))).^.5;
pa.qint = sign(un.U((N-1)*M+1))*mean(q);


%%% Stratification, either through direct calculation or numerical
%%% integration 
Psi_base = c*un.Y;
if fe.Stratified==1
    drhodpsi = pa.strat_fun(un.Psi,c);
    drhodpsi_base = pa.strat_fun(Psi_base,c);
    un.rho = zeros(M,N);
   % rho_base = zeros(M,N);
    for i=1:M
        un.rho(i,:) = cumtrapz(un.Psi(i,:),drhodpsi(i,:));
        rho_base(i,:) = cumtrapz(Psi_base(i,:),drhodpsi_base(i,:));
    end

    % This is knowing the function of stratification
    % and integrating analytically.
    c
    rho2 = pa.strat(un.Psi,un.c);

    
    %rho_base = rho_base-rho_base(1,1);
    rho2 = rho2-rho2(1,1);
    un.rho2=rho2;
    %un.rho_base = rho_base;
end


% ----------------------------------------------------------------------- %
%            Energetics of solution - ONLY FOR RIGID LID ATM
% ----------------------------------------------------------------------- %
if fe.Freesurface==0
X=un.X; Y=un.Y; U=un.U-un.c; V=un.V; rho=un.rho; Psi=un.Psi;
velocity_magnitude = (U.^2 + V.^2).^.5;

KE_local = 1/2 .* velocity_magnitude.^2;
PE_local = zeros(M,N);
for i=1:M;
    for j=2:N-1
       %['i=' num2str(i) '_j=' num2str(j)]
        % Take values at this point
        x0 = X(i,j); y0 = Y(i,j);
        rho0 = rho2(i,j); psi0 = Psi(i,j);
        % Mesh between y0 and psi0/c, and interpolate rho for these values
        y_pe = linspace(y0,psi0/c,50);
        rho_pe = interp1(Y(i,:),-rho2(i,:)+rho0,y_pe,'cubic');
        % PE is the integral of rho(z)-rho0 over this range 
        PE_local(i,j) = trapz(y_pe,rho_pe);
    end
end
KE_Y=zeros(1,M);
PE_Y=zeros(1,M);
for i=1:M
    KE_Y(i) = KE_Y(i) + trapz(Y(i,:),KE_local(i,:));
    PE_Y(i) = PE_Y(i) + trapz(Y(i,:),PE_local(i,:));
end
KE = trapz(X(:,1),KE_Y);
PE = trapz(X(:,1),PE_Y);
TE=KE+PE;
un.KE=KE; un.PE=PE; un.TE=TE; un.PE_local=PE_local; un.KE_local=KE_local;
%KE = trapz(X(:,1),trapz(Y(1,:),KE_local));
%PE = trapz(X(:,1),trapz(Y(1,:),PE_local));
%TE=KE+PE;
%}

%%% Total mass of the wave?
un.mass = trapz(t,(trapz(s,rho2.*pa.alpha.^2))) - trapz(t,(trapz(s,rho_base.*pa.alpha.^2)));
end



%%% Update some of the parameters that may have not been fixed/edited
%%% thoughout the solver
if fe.Freesurface==1
un.Area = trapz(me.s,y(Ey)-y(Ey(end)));
un.Amp = -(y(fe.Amploc)-y(fe.Amploc-M+1)); 

else
un.Area= sum((psi/c-y).^2)*(me.ds*me.dt)*pa.alpha^2  ;
un.Amp = -un.psi(fe.Amploc-M+1)/un.c+un.y(fe.Amploc-pa.M+1) - ...
        (-un.psi(fe.Amploc)/un.c+un.y(fe.Amploc));
end








