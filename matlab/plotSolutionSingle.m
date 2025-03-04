%% load data
if ~exist('data.dat', 'file')
    error('No data file, please compile and run the code first!');
end

fp = fopen('data.dat', 'r');

% read header
dimensions = fread(fp, [1,2], 'uint32')

solution = fread(fp, dimensions, 'float32');

fclose(fp);

%% plot solution
contourf(solution);