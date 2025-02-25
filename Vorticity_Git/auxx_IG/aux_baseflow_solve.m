function  [psi,B] = aux_baseflow_solve(pa,ylin,Q,PLOT)

   


vort_fun = pa.vort_fun; H=pa.H; 
y = ylin; N=numel(ylin); dy=y(2)-y(1);

psi(1) = 0;
psi_d(1) = 0;
psi_yy = zeros(N,1);


%%% Solve Using shooting method
Hi=ones(N,1); Hi_D=Hi;
J=zeros(N,N);

% Range to try bisection of y'(0) in
alpha=5;
delta=10^(-8);
error=1;
its=0;
"Solving base-flow"
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

if PLOT==1

plot(psi_y,y)
xlabel('u')
ylabel('y')
title('baseflow')

end


B = .5*psi_y(N)^2 + pa.g;

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







