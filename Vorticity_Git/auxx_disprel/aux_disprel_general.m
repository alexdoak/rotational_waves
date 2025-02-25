%{
k = (2*pi)/pa.wavelength;
a = pa.v1+0.000000000000001;
b = pa.v2;
H=pa.H;
g=pa.g;

Qe=un.Q;
ce=un.U(end,end);
%}
clear all

%{
a=-20
b=0;
H=1;
g=1;
N=100;

vort_fun = @(psi,Q) a*psi+b;
vort_fun_diff = @(psi,Q) a;
%}


%{
a=10; b=25; d=1/10;
vort_fun = @(psi,Q) a*exp(-b/abs(Q)*(psi-d*Q).^2);
vort_fun_diff = @(psi,Q) -2*a*b/abs(Q)*(psi-d*Q)*exp(-b/abs(Q)*(psi-d*Q).^2);
%}



H=.6; g=1; N=200;
a=0; b=40; c=.5; d=10;
vort_fun =  @(psi,Q) a*tanh(b*(psi-Q*c)) + d;
vort_fun_diff = @(psi,Q) a*b./(cosh(b*(psi-Q*c )).^2);
%}

y=linspace(0,H,N)
k_list = 1;
Q_list=zeros(2,numel(k_list));

for p=1:numel(k_list)
k=k_list(p);


Q_l = linspace(-10,15,100);
ERR_l = Q_l*0;
Qits=0;
BREAK_l=[];;
for Q=Q_l
    Qits=Qits+1;
    Qits
    % Solve for base flow
    [psi,B,Ue,Uye,BREAK] = aux_disprel_general_BF_bisection(vort_fun,y,N,H,g,Q,0);
    
    % Use baseflow to solve f(y) with KBD either end
    [f,f_ye,f_yye,BREAK2] = aux_disprel_general_VF(vort_fun_diff,k,psi,y,N,H,g,Q,Ue,0);
    
    % Check error in DBC.
    Error_DBC = Ue*f_ye + Ue*Uye + g;
    ERR_l(Qits) =   Error_DBC ;
    BREAK_l = [BREAK_l,BREAK+BREAK2];

    % Check against affine vorticity gamma(psi) = a*psi + b
    % d = (a-k^2)^.5;
    % Cd = cot(d*H);
    % S = sin(a^.5 * H);
    % C = cos(a^.5 * H);
    % E = b/a*(1-cos(a^.5*H));
    % alpha = b/a;
    % beta = (Q + alpha*(1-C))/S;
    % M = a^.5 * (alpha*S - beta*C)/sin(d*H) ;
    % fe = M*sin(d*y);
    % 

  %  Error_DBC= 


end






end
plot(Q_l,ERR_l)

%Q_list(imag(Q_list)~=0)=nan;
%plot(k_list,Q_list,'ko')
hold on;

%clin = roots([1,b*tanh(k*H)/k ,-g*tanh(k*H)/k])';
%Qlin2 = clin*H + b*H^2 /2



%{
k2 = a^.5;
a0 = C*(C/H+a^.5*S);
a1 = -(b/a)*(2*S*C/H + a^(1/2)*(S^2-C^2));
a2 = (b/a)^2/H*S^2 - a^(1/2)*(b/a)^2*S*C - g/a;


a0 = C*(C/H+a^.5*S);
a1 = -(b/a)*(2*S*C/H + a^(1/2)*(S^2-C^2));
a2 = (b/a)^2*S^2/H - a^(1/2)*(b/a)^2*S*C - g/a;


beta2 = roots([a0,a1,a2]);
Qlin2 = beta2*S - (b/a)*(1-C);

plot(k2,Qlin2,'x')
%}



% Comparing with constant vorticity!
% clin = roots([1,b*tanh(k*H)/k ,-g*tanh(k*H)/k])';
% Qlin2 = clin*H + b*H^2 /2



%{
beta=beta(1);
Q=Qlin(1);
alpha=b/a;


M = a^.5*(alpha*S-beta*C)/sin(d*H);
f=M*sin(d*H);
f_y = M*d*cos(d*H);

psi_y = a^.5 * (-alpha*S + beta*C);
psi_yy = -a*(alpha*C + beta*S);

DBC0 = psi_y * f_y + (psi_y*psi_yy+g)
KBC = f + psi_y


DBC1 = M*d*cos(d*H)*a^.5*(beta*C-alpha*S) + a^(3/2)*(alpha^2*S*C + alpha*(S^2-C^2)*beta - S*C*beta^2) + g

DBC2 =  a* (-d*cot(d*H)*(beta*C-alpha*S)^2 +  a^(1/2)*(alpha^2*S*C + alpha*(S^2-C^2)*beta - S*C*beta^2) + g/a)

DBC3 = a * (-d*cot(d*H)*(C^2*beta^2-2*alpha*S*C*beta+ alpha^2*S^2)  +  a^(1/2)*(alpha^2*S*C + alpha*(S^2-C^2)*beta - S*C*beta^2) + g/a)

DBC4 = a * (beta^2 * (-d*cot(d*H)*C^2 + -S*C*a^(1/2)) ... 
          + beta^1 * (2*d*alpha*S*C*cot(d*H) + a^(1/2)*alpha*(S^2-C^2)) ...
          + beta^0 * (-d*cot(d*H)*alpha^2*S^2 + a^(1/2)*alpha^2*S*C + g/a)   )

% 
% 
% 
% y = linspace(0,H,100); dy=y(2)-y(1);
% f = M*sin(d*(y+H));
% 
% f_yy(2:N-1) = (f(3:N)+f(1:N-2)-2*f(2:N-1))/dy^2
% 
% 
% plot(f_yy(2:N-1)+d^2*f(2:N-1),y(2:N-1))
% 
% 
% 
%}







