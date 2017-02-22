function [optval, x, y, z] = runlp()
%% import data
    filename = 'importdata.xlsx';

    iInfo = readtable(filename, 'Sheet', 1);%, 'ReadRowNames', true);
    jInfo = readtable(filename, 'Sheet', 2);%, 'ReadRowNames', true);
    kInfo = readtable(filename, 'Sheet', 3);%, 'ReadRowNames', true);
    mInfo = readtable(filename, 'Sheet', 4);%, 'ReadRowNames', true);
    miscInfo = readtable(filename, 'Sheet', 5);

    pipeCostPerMetre = miscInfo.pipeCostPerMetre(1);
    flowPerCapitaPerDay = miscInfo.flowPerCapitaPerDay(1);
    
%% variable definition

    size_i = height(iInfo);
    size_j = height(jInfo);
    size_k = height(kInfo);
	size_m = height(mInfo);
        
    for iter_i = 1:size_i 
        for iter_j = 1:size_j
            var_x(iter_i, iter_j) = strcat('x_',num2str(iInfo.name(iter_i)), '_', jInfo.name(iter_j));
            var_y(iter_i, iter_j) = strcat('y_',num2str(iInfo.name(iter_i)), '_', jInfo.name(iter_j));
        end
    end
    
    for iter_j = 1:size_j
        for iter_k = 1:size_k
            var_z(iter_j, iter_k) = strcat('z_',jInfo.name(iter_j), '_', num2str(kInfo.name(iter_k)));
        end
    end
    var_x = reshape(var_x', 1,[]);
    var_y = reshape(var_y', 1,[]);
    var_z = reshape(var_z', 1,[]);
    vars = cat(2, var_x, var_y, var_z);
%     for iter_v = 1:length(vars)
%         eval([vars{iter_v}, ' = ', num2str(iter_v),';']);
%     end
    
    % determine dij: calculate distance between wastewater sources and potential cw sites
    dij = zeros(size_i, size_j);
    
    for iter_i = 1:size_i
        for iter_j = 1:size_j
            dij(iter_i,iter_j) = lldistkm([iInfo.latitude(iter_i) iInfo.longitude(iter_i)],...
                                          [jInfo.latitude(iter_j) jInfo.longitude(iter_j)]);
        end
    end
    
    dij = dij*1000; %convert from km to m

    % define fi: flow rate of wastewater sources
    flowRate = zeros(1,size_i);

    for iter_i = 1:size_i
        flowRate(iter_i) = iInfo.population(iter_i)*flowPerCapitaPerDay;
    end

    % define pollutant concentrations
    %  emi, at wastewater source
    emi = zeros(size_m,size_i);

    for iter_m = 1:size_m
        for iter_i = 1:size_i
            emi(iter_m,iter_i) = (mInfo.minSrcConc(iter_m) + mInfo.maxSrcConc(iter_m))/2;
        end
    end

    %  eminj, influent of cw
    eminj = zeros(size_m,size_j);

    for iter_m = 1:size_m
        for iter_j = 1:size_j
            eminj(iter_m,iter_j) = (mInfo.minSrcConc(iter_m) + mInfo.maxSrcConc(iter_m))/2;
        end
    end

    % kma, coeff 
    % tmj, treatment targets
    % cmj, background concentrations
    kma = zeros(1, size_m);
    tmj = zeros(size_m, size_j);
    cmj = zeros(size_m, size_j);

    for iter_m = 1:size_m
        kma(iter_m) = mInfo.kA(iter_m);
        for iter_j = 1:size_j
            tmj(iter_m,iter_j) = mInfo.treatmentTarget(iter_m);
            cmj(iter_m,iter_j) = mInfo.backgroundConc(iter_m);
        end
    end

    % capacity and area of cw j given design option k
    Qjk = zeros(size_j, size_k);
    Ajk = zeros(size_j, size_k);
    cjk = zeros(size_j, size_k);

    for iter_j = 1:size_j
        for iter_k = 1:size_k
            Qjk(iter_j,iter_k) = kInfo.Qk(iter_k);
            Ajk(iter_j,iter_k) = kInfo.Ak(iter_k);
            cjk(iter_j,iter_k) = kInfo.ck(iter_k);
        end
    end

%% variable types and boundaries
    start_y = (size_i*size_j+1);
    numVars = (size_i*size_j)+(size_i*size_j)+(size_j*size_k)+(size_i*size_j*size_k); 
    intcon = start_y:(numVars-(size_i*size_j*size_k)); % y and z are binary
    lb = zeros(numVars, 1);
    ub = ones(numVars, 1);
    
%% objective function 
    f = init(size_i,size_j,size_k);
    [x, y, z, s] = reset(size_i,size_j,size_k);
    y = pipeCostPerMetre*dij;
    z = cjk;
    f = consolidateVars(x,y,z,s);
    
%% constraints    
    %initialisation
    [Aeq, beq] = init(size_i, size_j, size_k);
    [A, b] = init(size_i, size_j, size_k);
    
    % pollutant treatment constraint
    %  pollutant decay factor
    wjkm = zeros(size_j, size_k, size_m);
    exponjk = -1.0*(Ajk./Qjk);
    for iter_m = 1:size_m
        wjkm(:,:,iter_m) = exp(exponjk.*kma(iter_m));
    end
    wjkm(:,1,:) = zeros(size_j,1,size_m); %overwrite /0 in option 0 for k
    
    %  pollutant initial concentration
%    constmj = bsxfun(@minus,eminj,cmj);
    %  pollutant treatment target
    netTargetmj = bsxfun(@minus,tmj,cmj);
    
    temp_i = 7;
    temp_m = 1;
    temp_j = 10;
    temp_k = 5;
    % coefficient calculation and assignment
    for iter_m = 1:temp_m
        for iter_j = 1:temp_j
            [Aeq_x, Aeq_y, Aeq_z, Aeq_s, beq_new] = reset(size_i, size_j, size_k);
            Aeq_z(iter_j,:) = cmj(iter_m,iter_j)*wjkm(iter_j,:,iter_m);
            temp = emi(iter_m,:)'*wjkm(iter_j,:,iter_m);
            Aeq_s(:,iter_j,:) = permute(temp, [1 3 2]);
            beq_new = netTargetmj(iter_m,iter_j);
            [Aeq, beq] = addNewConstraint(Aeq_x, Aeq_y, Aeq_z, Aeq_s, beq_new, Aeq, beq);
        end
    end
    
    % constraint for introduced variable s_ijk
    for iter_i = 1:temp_i
        for iter_j = 1:temp_j
            for iter_k = 1:temp_k
                [A_x, A_y, A_z, A_s, b_new] = reset(size_i, size_j, size_k);                
                A_s(iter_i, iter_j, iter_k) = -1;
                A_x(iter_i, iter_j) = 1;
                A_z(iter_j, iter_k) = 1;
                b_new = 1;
                [A, b] = addNewConstraint(A_x, A_y, A_z, A_s, b_new, A, b);
            end
        end
    end
    
    for iter_i = 1:temp_i
        for iter_j = 1:temp_j
            for iter_k = 1:temp_k
                [A_x, A_y, A_z, A_s, b_new] = reset(size_i, size_j, size_k);                
                A_s(iter_i, iter_j, iter_k) = 1;
                A_z(iter_j, iter_k) = -1;
                b_new = 0;
                [A, b] = addNewConstraint(A_x, A_y, A_z, A_s, b_new, A, b);
            end
        end
    end
    
    for iter_i = 1:temp_i
        for iter_j = 1:temp_j
            for iter_k = 1:temp_k
                [A_x, A_y, A_z, A_s, b_new] = reset(size_i, size_j, size_k);                
                A_s(iter_i, iter_j, iter_k) = 1;
                A_x(iter_i, iter_j) = -1;
                b_new = 0;
                [A, b] = addNewConstraint(A_x, A_y, A_z, A_s, b_new, A, b);
            end
        end
    end
    
    % pollutant amount conservation constraint
    
    for iter_m = 1:size_m
        [Aeq_x, Aeq_y, Aeq_z, Aeq_s, beq_new] = reset(size_i, size_j, size_k);
        for iter_j = 1:size_j
            Aeq_x(:,iter_j) = flowRate.*emi(iter_m,:);
        end
        beq_new = sum(sum(flowRate.*emi(iter_m,:)));
        [Aeq, beq] = addNewConstraint(Aeq_x, Aeq_y, Aeq_z, Aeq_s, beq_new, Aeq, beq);
    end
    
    % CW minimum capacity constraint
    
    for iter_j = 1:size_j
        [A_x, A_y, A_z, A_s, b_new] = reset(size_i, size_j, size_k);
        A_x(:,iter_j) = flowRate;
        A_z(iter_j,:) = -1.0*Qjk(iter_j,:);
        b_new = 0;
        [A, b] = addNewConstraint(A_x, A_y, A_z, A_s, b_new, A, b);
    end
    
    % logic constraints
    % xij
    for iter_i = 1:size_i
        [Aeq_x, Aeq_y, Aeq_z, Aeq_s, beq_new] = reset(size_i, size_j, size_k);
        Aeq_x(iter_i,:) = ones(1,size_j);
        beq_new = 1;
        [Aeq, beq] = addNewConstraint(Aeq_x, Aeq_y, Aeq_z, Aeq_s, beq_new, Aeq, beq);
    end
    % zjk
    for iter_j = 1:size_j
        [Aeq_x, Aeq_y, Aeq_z, Aeq_s, beq_new] = reset(size_i, size_j, size_k);
        Aeq_z(iter_j,:) = ones(1,size_k);
        beq_new = 1;
        [Aeq, beq] = addNewConstraint(Aeq_x, Aeq_y, Aeq_z, Aeq_s, beq_new, Aeq, beq);
    end
    %yij
    for iter_j = 1:size_j
        for iter_i = 1:size_i
            [A_x, A_y, A_z, A_s,b_new] = reset(size_i, size_j, size_k);
            A_x(iter_i,iter_j) = 1;
            A_y(iter_i,iter_j) = -1;
            b_new = 0;
            [A, b] = addNewConstraint(A_x, A_y, A_z, A_s, b_new, A, b);
        end
    end

%% run intlinprog
    options = optimoptions('intlinprog','Display','iter', 'MaxNodes', 10^16, 'LPMaxIterations', Inf);
    [res,~,exitflag] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options)
%     x = reshape(res(1:size_i*size_j)', size_i,size_j)
%     for iter_v = 1:length(vars)
%         fprintf('%12.2f \t%s\n',res(iter_v),vars{iter_v})    
%     end
end

function [A, b] = init(size_i, size_j, size_k)
    num_x = size_i*size_j;
    num_y = size_i*size_j;
    num_z = size_j*size_k;
    num_s = size_i*size_j*size_k;
    A = double.empty(0,num_x+num_y+num_z+num_s);
    b = double.empty(0,1);
end

function [A_x, A_y, A_z, A_s, b] = reset(size_i, size_j, size_k)
    A_x = zeros(size_i,size_j);
    A_y = zeros(size_i,size_j);
    A_z = zeros(size_j,size_k);
    A_s = zeros(size_i,size_j,size_k);
    b = 0;
end

function [A_new, b_new] = addNewConstraint(A_x, A_y, A_z, A_s, b, A_curr, b_curr)
    A = cat(2,reshape(A_x',1,[]),reshape(A_y',1,[]),reshape(A_z',1,[]),reshape(permute(A_s,[2 1 3]),1,[]));
    A_new = cat(1, A_curr, A);
    b_new = cat(1, b_curr, b);
end

function f = consolidateVars(x,y,z,s)
    f = addNewConstraint(x,y,z,s,0,[],[]);
end