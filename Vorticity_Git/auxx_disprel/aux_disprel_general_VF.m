function  [f,f_yend,f_yyend,BREAK] = aux_disprel_general_VF(vort_fun_diff,k,psi,y,N,H,g,Q,Ue,PLOT)
BREAK=0;
dy=y(2)-y(1);

f(1) = 0;
f_d(1) = 0;
f_yy = zeros(N,1);


%%% Solve Using shooting method
Hi=ones(N,1); Hi_D=Hi;
J=zeros(N,N);

% Range to try bisection of y'(0) in
alpha=0;
delta=10^(-8);
error=1;
its=0;
%"Solving vertical function f(y)"
while abs(error)>10^(-12) & its<100
    its=its+1;

    % Integrate for psi(y)
    f(2) = dy*alpha(1);
    for p=3:N
        f(p) = 2*f(p-1) - f(p-2) - dy^2*(vort_fun_diff(psi(p-1),Q)-k^2)*f(p-1);
    end

    % Add delta to f'(0), integrate again for f
    f_d(2) = dy*(alpha(1)+delta);
    for p=3:N
        f_d(p) = 2*f_d(p-1) - f_d(p-2) -  dy^2*(vort_fun_diff(psi(p-1),Q)-k^2)*f_d(p-1);
    end

    error=(f(N)+Ue);
    error_d=(f_d(N)+Ue);
    J=(error_d-error)/delta;
    alpha = alpha - error/J;

end

if error>10^(-7) | its>99
   "Failed to solve for baseflow using shooting method, consider continuation methods from irrotational flow"
   BREAK=1;
end

if PLOT==1
f_y(1) =  [-3,4,-1]/(2*dy) * f(1:3)';
f_y(2:N-1) = (f(3:N)-f(1:N-2))/(2*dy);
f_y(N) = [1,-4,3]/(2*dy) * f(N-2:N)';

plot(f_y,y)
xlabel('u')
ylabel('y')
title('baseflow')

f_yy(2:N-1) = (f(3:end)+f(1:end-2)-2*f(2:end-1))/(dy^2)
f_yy(N) =  [-1,4,-5,2]/(dy^2)*f(N-3:N)'
end


f_y(N) = [1,-4,3]/(2*dy) * f(N-2:N)'; 
f_yy(N) =  [-1,4,-5,2]/(dy^2)*f(N-3:N)';
f_yend=f_y(N); 
f_yyend = f_yy(N);



end

% 
%{
% ----------------------------------------------------------------------- %
% Below is an exact solution for vort_fun = v1*psi + v2, to test method 
% ----------------------------------------------------------------------- %
v1=pa.v1; v2=pa.v2; H=pa.H;
A = v2/v1
B = (Q+v2/v1*(1-cos(v1^.5*H)))/sin(v1^.5*H)

psie = A * cos(v1^.5 * y) + B * sin(v1^.5 * y) - v2/v1;
psie_yye = -A * v1* cos(v1^.5 * y) - B * v1 * sin(v1^.5 * y)
psie_ye = -A * v1^.5* sin(v1^.5 * y) + B * v1^.5 * cos(v1^.5 * y)


psie_yy=zeros(1,N);
psie_yy(2:end-1) = (psie(3:end)+psie(1:end-2)-2*psie(2:end-1))/dy^2;




psi_yy=zeros(1,N);
psi_yy(2:end-1) = (psi(3:end)+psi(1:end-2)-2*psi(2:end-1))/dy^2;

%}

% 







