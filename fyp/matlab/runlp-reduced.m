function [optval, x, y, z] = runlp()

%% run intlinprog
    options = optimoptions('intlinprog','Display','iter',...
                           'MaxNodes', 10^16,...
                           'LPMaxIterations', Inf,...
                           'MaxTime', 14400);
    [res,~,exitflag] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);

    x = transpose(reshape(res(1:end_x), size_j, size_i))
    y = transpose(reshape(res((end_x+1):end_y), size_j, size_i))
    z = transpose(reshape(res((end_y+1):end_z), size_k, size_j))
    s = permute(reshape(res((end_z+1):numVars), size_j, size_i, size_k), [2 1 3])
    
    exitflag
    
%     x = reshape(res(1:size_i*size_j)', size_i,size_j)
%     for iter_v = 1:length(vars)
%         fprintf('%12.2f \t%s\n',res(iter_v),vars{iter_v})    
%     end
end
