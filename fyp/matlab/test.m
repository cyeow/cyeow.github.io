%setup
% importing tables
filename = 'importdata.xlsx';

iInfo = readtable(filename,...
    'Sheet', 1,... 
    'ReadRowNames', true);
jInfo = readtable(filename,...
    'Sheet', 2,...
    'ReadRowNames', true);
kInfo = readtable(filename,...
    'Sheet', 3,...
    'ReadRowNames', true);
mInfo = readtable(filename,...
    'Sheet', 4,...
    'ReadRowNames', true);
pipeCostPerMetre = readtable(filename,...
    'Sheet', 5);
flowPerCapitaPerDay = importdata('flow-per-capita-per-day.csv');

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

%  emoutj, effluent of cw
emoutj = zeros(height(mInfo),height(jInfo));

% kma, coeff 
kma = ones(1, height(mInfo));

% tmj, treatment targets
% cmj, background concentrations
tmj = zeros(height(mInfo), height(jInfo));
cmj = zeros(height(mInfo), height(jInfo));

for iter_m = 1:height(mInfo)
    for iter_j = 1:height(iInfo)
        tmj(iter_m,iter_j) = mInfo.treatmentTarget(iter_m);
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


% defining constraints

% numCols remain the same for the constraint arrays.
% from 1 to i*j represent xij, from i*j+1 to i*j+1+i*j represent yij,
% from i*j+1+i*j+1 to i*j+1+i*j+1+j*k represent zjk
% numCols = height(iInfo)*height(jInfo)+height(iInfo)*height(jInfo)+height(jInfo)*height(kInfo);

%  inequalities
% rows are determined by the number of inequality constraints
numRows = height(jInfo)+height(jInfo)+(height(jInfo)*height(iInfo));

% A = zeros(numRows, numCols);
% 
% for iter = 1:numRows
%     temp = zeros(1,numCols);
%     
%     if  iter <= height(jInfo)
%         for iter_j = 1:height(jInfo)
%             for iter_m = 1:height(mInfo)
%                 sum = 0;
%                 for iter_i = 1:height(iInfo)
%                     sum = sum + flowRate(iter_i)*emi(iter_m,iter_i) - eminj(iter_m,iter_j)*flowRate(iter_i);
%                 end
%                 temp(height(iInfo)*height(jInfo)+iter_i+iter_j) = sum;
%                 A(iter,:) = temp;
%             end
%         end
%     elseif iter <= height(jInfo)+height(jInfo)
%         
%     elseif iter <= height(jInfo)+height(jInfo)+(height(jInfo)*height(iInfo))
%         
%     end
%     
%     %add line to array A
% end

% i is row, j is col
% yij = zeros(1,height(jInfo));
% 
% for iter_m = 1:height(mInfo)
%         sum = 0;
%     for iter_j = 1:height(jInfo)
%         for iter_i = 1:height(iInfo)
%             sum = sum + flowRate(iter_i)*emi(iter_m,iter_i) - eminj(iter_m,iter_j)*flowRate(iter_i);
%         end
%         yij(1,iter_j) = sum;
%     end
%     
% end
%        A(iter,:) = temp;

resetA(height(iInfo), height(jInfo), height(kInfo));
% pollutant mass equality constraints
for iter_m = 1:height(mInfo)
    %sum all over 
    for iter_i = 1:height(iInfo)
        coeff = 0;
        for iter_j = 1:height(jInfo)
            newAeq_x(iter_i,iter_j)  = coeff + emi(iter_m, iter_i)*flowRate(iter_i);
        end
         
        newbeq = emi(iter_m, iter_i)*flowRate(iter_i);
        %concat to Aeq and beq
    end
end

resetA(height(iInfo), height(jInfo), height(kInfo));
% treatment target constraints
for iter_m = 1:height(mInfo)
    for iter_j = 1:height(jInfo)
        newb = tmj(iter_m, iter_j) + cmj(iter_m,iter_j);
        
        for iter_k = 1:height(kInfo)
            newA_z(iter_j,iter_k) = exp(-1.0*Ajk(iter_j,iter_k)/Qjk(iter_j,iter_k)*kma(iter_m));            
        end
        %concat A and  b
    end
end
        
% capacity of constructed wetlands constraints
resetA(height(iInfo), height(jInfo), height(kInfo));
for iter_j = 1:height(jInfo)
    for iter_i = 1:height(iInfo)
        newA_x(iter_i, iter_j) = flowRate(iter_i);
    end 
    for iter_k = 1:height(kInfo)
        newA_z(iter_j,iter_k) = -1.0*Qjk(iter_j,iter_k);
    end
    newb = 0;
    %concat
end

function resetA(height_i,height_j, height_k)
    newA_x = zeros(height_i, height_j);
    newA_y = zeros(height_i, height_j);
    newA_z = zeros(height_j, height_k);
    newAeq_x = zeros(height_i, height_j);
    newAeq_y = zeros(height_i, height_j);
    newAeq_z = zeros(height_j, height_k);
end