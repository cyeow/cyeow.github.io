filename = 'importdata.xlsx';

iInfo = readtable(filename, 'Sheet', 1, 'ReadRowNames', true);
jInfo = readtable(filename, 'Sheet', 2, 'ReadRowNames', true);
kInfo = readtable(filename, 'Sheet', 3, 'ReadRowNames', true);
mInfo = readtable(filename, 'Sheet', 4, 'ReadRowNames', true);
miscInfo = readtable(filename, 'Sheet', 5);