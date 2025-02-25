function  [psi,B,Ue,Uye,BREAK] = aux_disprel_general_BF_bisection(vort_fun,y,N,H,g,Q,PLOT)
BREAK=0;
dy=y(2)-y(1);

psi1(1) = 0;
psi2(1) = 0;
psi3(1) = 0;
psi_yy = zeros(N,1);


%%% Solve Using shooting method
Hi=ones(N,1); Hi_D=Hi;
J=zeros(N,N);

% Range to try bisection of y'(0) in
alpha=[-100,0,100];
delta=10^(-8);
error1=1; error2=1; error3=1;
its=0;
%"Solving base-flow"
while (abs(error1)>10^(-12) & abs(error2)>10^(-12) & abs(error3)>10^(-12))   & its<100
    its=its+1;
    
    % Integrate for psi(y) with lower alpha
    psi1(2) = dy*alpha(1);
    for p=3:N
        psi1(p) = 2*psi1(p-1) - psi1(p-2) - dy^2*vort_fun(psi1(p-1),Q);
    end

    % Integrate for psi(y) with middle alpha 
    psi2(2) = dy*(alpha(2)+delta);
    for p=3:N
        psi2(p) = 2*psi2(p-1) - psi2(p-2) - dy^2*vort_fun(psi2(p-1),Q);
    end

  % Integrate for psi(y) with upper alpha 
    psi3(2) = dy*(alpha(3)+delta);
    for p=3:N
        psi3(p) = 2*psi3(p-1) - psi3(p-2) - dy^2*vort_fun(psi3(p-1),Q);
    end

    

    error1=(psi1(N)-Q);
    error2=(psi2(N)-Q);
    error3=(psi3(N)-Q);

    if error1 * error2 < 0
        alpha=[alpha(1),(alpha(1)+alpha(2))/2,alpha(2)];
    elseif error2 * error3 < 0
        alpha=[alpha(2),(alpha(2)+alpha(3))/2,alpha(3)];
    end
   
end

if  (abs(error1)>10^(-7) & abs(error2)>10^(-7) & abs(error3)>10^(-7)) |  its>99
   "Failed to solve for baseflow using shooting method, consider continuation methods from irrotational flow"
   BREAK=1;
end

[~,I] = min(abs([error1,error2,error3]));
psi = eval(['psi' num2str(I)]);


if PLOT==1
psi_y(1) =  [-3,4,-1]/(2*dy) * psi(1:3)';
psi_y(2:N-1) = (psi(3:N)-psi(1:N-2))/(2*dy);
psi_y(N) = [1,-4,3]/(2*dy) * psi(N-2:N)';

plot(psi_y,y)
xlabel('u')
ylabel('y')
title('baseflow')

psi_yy(2:N-1) = (psi(3:end)+psi(1:end-2)-2*psi(2:end-1))/(dy^2)
psi_yy(N) =  [-1,4,-5,2]/(dy^2)*psi(N-3:N)'
end


psi_y(N) = [1,-4,3]/(2*dy) * psi(N-2:N)'; 
psi_yy(N) =  [-1,4,-5,2]/(dy^2)*psi(N-3:N)';
Ue = psi_y(N); 
Uye = psi_yy(N);

B = .5*psi_y(N)^2 + g;

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







