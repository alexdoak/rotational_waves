
H=1;
N=100;
Q=0.01;


 for v1=0:.1:10
    
     v2 = 25
     v3 = 1/10;
vort_fun = @(psi,Q) v1 * exp(-v2/abs(Q)*(psi-v3*Q).^2)
y = linspace(0,H,N)'; dy=y(2)-y(1);


psi = Q/H*y;
psi_yy = zeros(N,1);


%%% Solve
Hi=ones(N,1); Hi_D=Hi;
J=zeros(N,N);
while max(abs(Hi))>10^(-6)

    psi_yy(2:end-1) = (psi(3:end)+psi(1:end-2)-2*psi(2:end-1))/dy^2;

    Hi(1) = psi(1);
    Hi(2:end-1) = psi_yy(2:end-1) + vort_fun(psi(2:end-1),Q);
    Hi(N) = psi(N)-Q;

    STORE=psi;
    for p=1:N
        psi=STORE;
        psi(p)=psi(p)+delta;
        
        psi_yy(2:end-1) = (psi(3:end)+psi(1:end-2)-2*psi(2:end-1))/dy^2;

        Hi_D(1) = psi(1);
        Hi_D(2:end-1) = psi_yy(2:end-1) + vort_fun(psi(2:end-1),Q);
        Hi_D(N) = psi(N)-Q;

        J(:,p) = (Hi_D-Hi)/delta;
    end
    
    psi=STORE;
    psi = psi - J\Hi;
    
    Residuals=max(abs(Hi))
end

psi_yy(2:end-1) = (psi(3:end)+psi(1:end-2)-2*psi(2:end-1))/dy^2;
psi_y = zeros(N,1);
psi_y(1) = [-3     4    -1]/(2*dy) * psi(1:3);
psi_y(2:end-1) = (psi(3:end)-psi(1:end-2))/(2*dy);
psi_y(N) = [1    -4     3]/(2*dy) * psi(end-2:end);



 end