function [pa,un,me,fe,delta] = aux_create_IG_vort(pa,un,me,fe)


%%%% PARAMETERS OF FLOW
pa.M=250; pa.N=200; 
%pa.g=1; 
%pa.wavelength=2*pi;
pa.strat_fun =  @(psi,c) 0;
pa.strat = @(psi,c) 0;
if fe.Bathymetry==1
pa.b1 = 0.02;
pa.b2 = 3;
pa.bath_fun = @(s) pa.b1*exp(-pa.b2*(s).^2) + 1/2*pa.b1*exp(-pa.b2*(s-2).^2) 
end

if fe.Freesurface==0
fe.Amploc = pa.M*round(pa.s3/pa.H*pa.N);
else 
fe.Amploc=pa.M*pa.N;
end

%%% Unpack structs
M=pa.M; N=pa.N; g=pa.g; H=pa.H; wavelength=pa.wavelength; v1=pa.v1; v2=pa.v2;

%%% Assuming the mapping is from box to box: given lambda, L, H, what is d?
pa.L=2*pi; pa.d = H*pa.L/wavelength; pa.alpha = pa.wavelength/pa.L;

%%% Construct mesh for IG
MN=M*N;
[me] = aux_mesh_vort(pa,me,fe);  
T=me.T; s=me.s;

%%% Recover base flow, where we choose Q at top streamline.
Y = (T+pa.d)*H/pa.d; 
[psi,B] = aux_baseflow_solve(pa,Y(1,:),un.Q,1);
Psi = repmat(psi,M,1);


for j=1:N
    Y((j-1)*M+1:j*M) = Y((j-1)*M+1:j*M) + (j-1)/(10000*N) * cos(2*pi/(2*pi)*s);
    Psi((j-1)*M+1:j*M) = Psi((j-1)*M+1:j*M) + (j-1)/(10000*N) * cos(2*pi/(2*pi)*s);
end

y = reshape(Y,[],1);
psi = reshape(Psi,[],1);

%%% Pack into structs
un.Y=Y; un.Psi=Psi; un.y=y; un.psi=psi;  un.B=B;

un.Area = 0;
un.Amp=0.01;





Amp=un.Amp; Area=un.Area; Q=un.Q;
clear L; clear d; clear g;
delta=10^(-11); 


