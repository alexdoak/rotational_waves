% ----------------------------------------------------------------------- %
%                Create new lists to store data                           %
% ----------------------------------------------------------------------- %
function [li] = aux_newlists(a,pa,un,fe)

li.a_list = {};
li.qint_list = [];
li.amp_list = [];
li.M_list = [];
li.N_list=[];
li.Q_list=[];
li.L_list = [];
li.d_list = [];
li.speed_list=[];
li.Burns_list=[];

if fe.Freesurface==0
    li.TE_list=[];
    li.KE_list=[];
    li.PE_list=[];
    li.mass_list=[];
end