function [li] = aux_list_truncate(li,CNT)

F = fieldnames(li);
for p = 1:numel(F)
    F = fieldnames(li);
    F=F(p);
    F=F{1};

    call = ['li.', F '=li.', F, '(1:', num2str(CNT) ');' ];
    eval(call);

    %if isequal(F, 'a_list')
    %    call
    %    exec(call)
    %else 
    %end



end
clear call
clear F
