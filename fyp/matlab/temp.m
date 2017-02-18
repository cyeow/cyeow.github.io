function [optval, x, y, z] = temp() 
    filename = 'importdata.xlsx';

    iInfo = readtable(filename, 'Sheet', 1, 'ReadRowNames', true);
    jInfo = readtable(filename, 'Sheet', 2, 'ReadRowNames', true);
    kInfo = readtable(filename, 'Sheet', 3, 'ReadRowNames', true);
    mInfo = readtable(filename, 'Sheet', 4, 'ReadRowNames', true);
    miscInfo = readtable(filename, 'Sheet', 5);

    pipeCostPerMetre = miscInfo.pipeCostPerMetre(1);
    flowPerCapitaPerDay = miscInfo.flowPerCapitaPerDay(1);

    size_i = height(iInfo);
    size_j = height(jInfo);
    size_k = height(kInfo);
	size_m = height(mInfo);

    % determine dij: calculate distance between wastewater sources and potential cw sites
    dij = zeros(size_i, size_j);

    for iter_i = 1:size_i
        for iter_j = 1:size_j
            dij(iter_i,iter_j) = distance(iInfo.latitude(iter_i), iInfo.longitude(iter_i),...
                                          jInfo.latitude(iter_j), jInfo.longitude(iter_j));
        end
    end

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

    %constraints initialisation
    [Aeq, beq, A, b] = init(size_i, size_j, size_k);
    
    % pollutant treatment constraint
    
    % pollutant decay factor
    wjkm = zeros(size_j, size_k, size_m);
    tempjk = exp(-1.0*(Ajk./Qjk));
    tempjk(:,1) = zeros(size_j,1); %overwrite /0 in option 0 for k
    for iter_m = 1:size_m
        wjkm(:,:,iter_m) = tempjk.*kma(iter_m);
    end
    
    % pollutant initial concentration
    constmj = bsxfun(@minus,eminj,cmj);
    % pollutant treatment target
    netTargetmj = bsxfun(@minus,tmj,cmj);
    
    % coefficient calculation and assignment
    for iter_m = 1:size_m
        for iter_j = 1:size_j
            [Aeq_x, Aeq_y, Aeq_z, beq, A_x, A_y, A_z, b] = reset(size_i, size_j, size_k);
            Aeq_z(iter_j,:) = constmj(iter_m,iter_j)*wjkm(iter_j,:,iter_m);
            beq = netTargetmj(iter_m,iter_j);
            Aeq = addNewConstraint(Aeq_x, Aeq_y, Aeq_z, beq, Aeq);
        end
    end
    % mass conservation
    
    
end

function [Aeq, beq, A, b] = init(size_i, size_j, size_k)
    Aeq = zeros(1,size_i*size_j+size_i*size_j+size_j*size_k);
    beq = 0;
    A = zeros(1,size_i*size_j+size_i*size_j+size_j*size_k);
    b = 0;
end

function [Aeq_x, Aeq_y, Aeq_z, beq, A_x, A_y, A_z, b] = reset(size_i, size_j, size_k)
    Aeq_x = zeros(size_i,size_j);
    Aeq_y = zeros(size_i,size_j);
    Aeq_z = zeros(size_j,size_k);
    A_x = zeros(size_i,size_j);
    A_y = zeros(size_i,size_j);
    A_z = zeros(size_j,size_k);
    beq = 0;
    b = 0;
end

function A_new = addNewConstraint(A_x, A_y, A_z, b, A_curr)
    A = cat(2,reshape(A_x',1,[]),reshape(A_y',1,[]),reshape(A_z',1,[]));
    A_new = cat(1, A_curr, A);
end