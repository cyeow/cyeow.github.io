    filename = 'importdata.xlsx';

    iInfo = readtable(filename, 'Sheet', 1, 'ReadRowNames', true);
    jInfo = readtable(filename, 'Sheet', 2, 'ReadRowNames', true);
    kInfo = readtable(filename, 'Sheet', 3, 'ReadRowNames', true);
    mInfo = readtable(filename, 'Sheet', 4, 'ReadRowNames', true);
    miscInfo = readtable(filename, 'Sheet', 5);

    pipeCostPerMetre = miscInfo.pipeCostPerMetre(1);
    flowPerCapitaPerDay = miscInfo.flowPerCapitaPerDay(1);

    i_size = height(iInfo);
    j_size = height(jInfo);
    k_size = height(kInfo);
    m_size = height(mInfo);

    % determine dij: calculate distance between wastewater sources and potential cw sites
    dij = zeros(height(iInfo), height(jInfo));

    for iter_i = 1:height(iInfo)
        for iter_j = 1:height(jInfo)
            dij(iter_i,iter_j) = distance(iInfo.latitude(iter_i), iInfo.longitude(iter_i),...
                                          jInfo.latitude(iter_j), jInfo.longitude(iter_j));
        end
    end

    % define fi: flow rate of wastewater sources
    flowRate = zeros(1,height(iInfo));

    for iter_i = 1:height(iInfo)
        flowRate(iter_i) = iInfo.population(iter_i)*flowPerCapitaPerDay;
    end

    % define pollutant concentrations
    %  emi, at wastewater source
    emi = zeros(height(mInfo),height(iInfo));

    for iter_m = 1:height(mInfo)
        for iter_i = 1:height(iInfo)
            emi(iter_m,iter_i) = (mInfo.minSrcConc(iter_m) + mInfo.maxSrcConc(iter_m))/2;
        end
    end

    %  eminj, influent of cw
    eminj = zeros(height(mInfo),height(jInfo));

    for iter_m = 1:height(mInfo)
        for iter_j = 1:height(jInfo)
            eminj(iter_m,iter_j) = (mInfo.minSrcConc(iter_m) + mInfo.maxSrcConc(iter_m))/2;
        end
    end

    % kma, coeff 
    % tmj, treatment targets
    % cmj, background concentrations
    kma = zeros(1, height(mInfo));
    tmj = zeros(height(mInfo), height(jInfo));
    cmj = zeros(height(mInfo), height(jInfo));

    for iter_m = 1:height(mInfo)
        kma(iter_m) = mInfo.kA(iter_m);
        for iter_j = 1:height(iInfo)
            tmj(iter_m,iter_j) = mInfo.treatmentTarget(iter_m);
            cmj(iter_m,iter_j) = mInfo.backgroundConc(iter_m);
        end
    end

    % capacity and area of cw j given design option k
    Qjk = zeros(height(jInfo), height(kInfo));
    Ajk = zeros(height(jInfo), height(kInfo));
    cjk = zeros(height(jInfo), height(kInfo));

    for iter_j = 1:height(jInfo)
        for iter_k = 1:height(kInfo)
            Qjk(iter_j,iter_k) = kInfo.Qk(iter_k);
            Ajk(iter_j,iter_k) = kInfo.Ak(iter_k);
            cjk(iter_j,iter_k) = kInfo.ck(iter_k);
        end
    end

    wjkm = zeros(height(jInfo), height(kInfo), height(mInfo));

    tempjk = exp(-1.0*(Ajk./Qjk));
    tempjk(:,1) = zeros(height(jInfo),1); %overwrite /0 in option 0 for k
    
    for iter_m = 1:height(mInfo)
        wjkm(:,:,iter_m) = tempjk.*kma(iter_m);
    end
    
    