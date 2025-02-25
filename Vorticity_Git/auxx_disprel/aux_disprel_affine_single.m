




function [Qlin] = aux_disprel_affine_single(a,b,H,g,k)



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
Qlin=sort(Qlin);



vort_fun = @(psi,Q)a*psi+b;
y = linspace(0,H,100); N=numel(y); dy=y(2)-y(1);
psi(1) = 0;
psi_d(1) = 0;
psi_yy = zeros(N,1);


figure()
Q=Qlin(1)
"Solving base-flow, Q-"
%%% Solve Using shooting method
Hi=ones(N,1); Hi_D=Hi;
J=zeros(N,N);

% Range to try bisection of y'(0) in
alpha=5;
delta=10^(-8);
error=1;
its=0;
while abs(error)>10^(-7) & its<100
    its=its+1;

    % Integrate for psi(y)
    psi(2) = dy*alpha(1);
    for p=3:N
        psi(p) = 2*psi(p-1) - psi(p-2) - dy^2*vort_fun(psi(p-1),Q);
    end

    % Add delta to psi'(0), integrate again for psi
    psi_d(2) = dy*(alpha(1)+delta);
    for p=3:N
        psi_d(p) = 2*psi_d(p-1) - psi_d(p-2) - dy^2*vort_fun(psi_d(p-1),Q);
    end

    error=(psi(N)-Q);
    error_d=(psi_d(N)-Q);
    J=(error_d-error)/delta;
    alpha = alpha - error/J;

end

if abs(error)>10^(-7) | its>99
   "Failed to solve for baseflow using shooting method, consider continuation methods from irrotational flow"
end

psi_y(1) =  [-3,4,-1]/(2*dy) * psi(1:3)';
psi_y(2:N-1) = (psi(3:N)-psi(1:N-2))/(2*dy);
psi_y(N) = [1,-4,3]/(2*dy) * psi(N-2:N)';
psi_yy(2:N-1) = (psi(3:end)+psi(1:end-2)-2*psi(2:end-1))/(dy^2);
psi_yy(N) =  [-1,4,-5,2]/(dy^2)*psi(N-3:N)';


subplot(1,2,1)
plot(psi_y,y)
hold on;
plot([0,0],[0,1],'k--')
xlabel('u')
ylabel('y')
title(['baseflow Q-=' num2str(round(Q,2))  ' for a=' num2str(a) ' b=' num2str(b)])




Q=Qlin(2)
"Solving base-flow, Q-"
%%% Solve Using shooting method
Hi=ones(N,1); Hi_D=Hi;
J=zeros(N,N);

% Range to try bisection of y'(0) in
alpha=5;
delta=10^(-8);
error=1;
its=0;
while abs(error)>10^(-7) & its<100
    its=its+1;

    % Integrate for psi(y)
    psi(2) = dy*alpha(1);
    for p=3:N
        psi(p) = 2*psi(p-1) - psi(p-2) - dy^2*vort_fun(psi(p-1),Q);
    end

    % Add delta to psi'(0), integrate again for psi
    psi_d(2) = dy*(alpha(1)+delta);
    for p=3:N
        psi_d(p) = 2*psi_d(p-1) - psi_d(p-2) - dy^2*vort_fun(psi_d(p-1),Q);
    end

    error=(psi(N)-Q);
    error_d=(psi_d(N)-Q);
    J=(error_d-error)/delta;
    alpha = alpha - error/J;

end

if abs(error)>10^(-7) | its>99
   "Failed to solve for baseflow using shooting method, consider continuation methods from irrotational flow"
end

psi_y(1) =  [-3,4,-1]/(2*dy) * psi(1:3)';
psi_y(2:N-1) = (psi(3:N)-psi(1:N-2))/(2*dy);
psi_y(N) = [1,-4,3]/(2*dy) * psi(N-2:N)';
psi_yy(2:N-1) = (psi(3:end)+psi(1:end-2)-2*psi(2:end-1))/(dy^2);
psi_yy(N) =  [-1,4,-5,2]/(dy^2)*psi(N-3:N)';


subplot(1,2,2)
plot(psi_y,y)
hold on;
plot([0,0],[0,1],'k--')
xlabel('u')
ylabel('y')
title(['baseflow Q+=' num2str(round(Q,2))  ' for a=' num2str(a) ' b=' num2str(b)])

end

