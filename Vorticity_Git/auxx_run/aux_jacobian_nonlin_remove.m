%********************************************************************%
%                Jacobian Computation 
% This code removes nonlinear contributions from the Jacobian. The
% motivation being that we don't then have to recompute the linear part.
%********************************************************************%

function [J] = aux_jacobian_nonlin_remove(J,pa,me,fe,ja)

% Pack necessary data into a struct
M=pa.M; N=pa.N; MN=M*N;
k=me.k; Sx=me.Sx; Sy=me.Sy; Ex=me.Ex; Ey=me.Ey; kt=me.kt;

FE_psi=ja.FE_psi;


% Remove nonlinear contributions to the Jacobian
J(kt,kt) =  J(kt,kt) - spdiags(FE_psi,0,numel(kt),numel(kt));

if fe.Freesurface==1 | fe.Bathymetry==1
    FE_1=ja.FE_1;FE_m1=ja.FE_m1;FE_M=ja.FE_M;FE_mM=ja.FE_mM;
    
    J(k,MN+k+1) = J(k,MN+k+1) - spdiags(FE_1,0,numel(k),numel(k));
    J(k,MN+k-1) = J(k,MN+k-1) - spdiags(FE_m1,0,numel(k),numel(k));
    J(kt,MN+kt+M) = J(kt,MN+kt+M) - spdiags(FE_M,0,numel(kt),numel(kt));
    J(kt,MN+kt-M) = J(kt,MN+kt-M) - spdiags(FE_mM,0,numel(kt),numel(kt));
    
    if fe.Freesurface==1
        BERN_p0=ja.BERN_p0;BERN_pmM=ja.BERN_pmM;BERN_pm2M=ja.BERN_pm2M;
        BERN_y0=ja.BERN_y0;BERN_ymM=ja.BERN_ymM;BERN_ym2M=ja.BERN_ym2M;
        BERN_y1=ja.BERN_y1;BERN_ym1=ja.BERN_ym1;
        
        J(MN+Ey,Ey) = J(MN+Ey,Ey) - spdiags(BERN_p0,0,M,M);
        J(MN+Ey,Ey-M) = J(MN+Ey,Ey-M) - spdiags(BERN_pmM,0,M,M);
        J(MN+Ey,Ey-2*M) = J(MN+Ey,Ey-2*M) - spdiags(BERN_pm2M,0,M,M);
        J(MN+Ey,MN+Ey) = J(MN+Ey,MN+Ey) - spdiags(BERN_y0,0,M,M);
        J(MN+Ey,MN+Ey-M) = J(MN+Ey,MN+Ey-M) - spdiags(BERN_ymM,0,M,M);
        J(MN+Ey,MN+Ey-2*M) = J(MN+Ey,MN+Ey-2*M) - spdiags(BERN_ym2M,0,M,M);
        J(MN+Ey,MN+Ey+1) = J(MN+Ey,MN+Ey+1) - spdiags(BERN_y1,0,M,M);
        J(MN+Ey,MN+Ey-1) = J(MN+Ey,MN+Ey-1) - spdiags(BERN_ym1,0,M,M);
    end
    
    if fe.Stratified == 1 
        FE_y=ja.FE_y;
        J(k,M*N+k) =  J(k,M*N+k) - spdiags(FE_y,0,numel(k),numel(k));
    end
    
    if fe.Fixwavelength == 1 
        WL_1=ja.WL_1;WL_2=ja.WL_2;WL_3=ja.WL_3;
        J(2*MN+2,MN+1:MN+M) = J(2*MN+2,MN+1:MN+M) - WL_1;
        J(2*MN+2,MN+M+1:MN+2*M) =  J(2*MN+2,MN+M+1:MN+2*M) - WL_2;
        J(2*MN+2,MN+2*M+1:MN+3*M) =  J(2*MN+2,MN+2*M+1:MN+3*M) - WL_3;
    end
end


