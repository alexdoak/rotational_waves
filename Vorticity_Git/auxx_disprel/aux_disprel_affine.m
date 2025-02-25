
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig = figure();
set(fig,'Units','centimeters','Position',[2 2 12 7]);


FONTSIZE = 10;


kend=10
delta=10^(-11);



a=120;
b=5;
H=1;
g=1;


a0_list = [];
a1_list=[];
a2_list=[];

%%%%%% Find where asymptotes are, so that we can define decent k ranges.
%syms kk 
%d = (a-kk^2)^.5;
%Cd = cot(d*H);
S = sin(a^.5 * H);
C = cos(a^.5 * H);
E = b/a*(1-cos(a^.5*H));
a2 = C*(d*C*Cd+a^.5*S);

fun = @(kk) (a-kk^2)^.5*C*cot((a-kk^2)^.5*H)+a^.5*S;
k_l = [];
K_list = dictionary;

%for its=1:50
%its
%k_l=[k_l,vpasolve(a2==0,kk,Random=true)];
%end

for k0 = 0.01:0.01:10
    K = fzero(fun,k0,optimset('FunValCheck', 'off', 'Display', 'off'));
    if isnan(K)==0
    k_l=[k_l,K];
    end
end

k_l=sort(k_l(find(double(k_l)>0)))

if numel(k_l)>0
    k_l2=k_l(1);
    for p=2:numel(k_l)
        if k_l(p)-k_l2(end)>10^(-5)
            k_l2=[k_l2,k_l(p)];
        end
    end
    
    
    k_l=k_l2; clear k_l2; 
    k_l(end+1)=kend;
    
    
    P = numel(k_l);
    for p = 1:P
        if p==1
            K_list{num2str(p)} = linspace(0,k_l(p)-delta,1000);
        else
            K_list{num2str(p)} = linspace(k_l(p-1)+delta,k_l(p)-delta,1000);
        end
    end
else
    P=1;
    K_list{"1"} = linspace(0,kend,1000)
end
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

Q_list = dictionary;
A0_list = dictionary;
A1_list = dictionary;
A2_list = dictionary;

for I=1:P
    I
    k_list=K_list{num2str(I)};
    q_list=zeros(2,numel(k_list));


for p=1:numel(k_list)


k=k_list(p);






d = (a-k^2)^.5;
Cd = cot(d*H);
S = sin(a^.5 * H);
C = cos(a^.5 * H);
E = b/a*(1-cos(a^.5*H));

a2 = C*(d*C*Cd+a^.5*S);
a1 = -(b/a)*(2*d*S*C*Cd + a^(1/2)*(S^2-C^2));
a0 = (b/a)^2*d*S^2*Cd - a^(1/2)*(b/a)^2*S*C - g/a;

%a1 = 2*E + b/(a^(1/2)*d*C)*S
%a2 = E^2-g*S^2/(a*d*C) + b/(a^.5*d*C) * E*S;

beta = roots([a2,a1,a0]);
Qlin = beta*S - (b/a)*(1-C);

a0_list=[a0_list,a0];
a1_list=[a1_list,a1];
a2_list=[a2_list,a2];
q_list(:,p) = sort(Qlin);
end

q_list(imag(q_list)~=0)=nan;

Q_list{num2str(I)} = q_list;
A0_list{num2str(I)} = a0_list;
A1_list{num2str(I)} = a1_list;
A2_list{num2str(I)} = a2_list;


plot(K_list{num2str(I)},Q_list{num2str(I)},'k')
hold on;
title(['$a=$' num2str(a) ' $b=$' num2str(b)],'FontSize',FONTSIZE)

end


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







