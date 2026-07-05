function data = case1_load_data(rootDir)
%case1_load_data Read Case 1 train/testA csv files into MATLAB structs.
%   The csv format follows the assignment statement:
%   odd rows  -> temperature values
%   even rows -> voltage values

if nargin < 1
    rootDir = pwd;
end

trainFile = find_first_file(rootDir, 'dataform_train*.csv');
testAFile = find_first_file(rootDir, 'dataform_testA*.csv');

data.train = read_one_dataset(trainFile, 'train');
data.testA = read_one_dataset(testAFile, 'testA');
end

function filePath = find_first_file(rootDir, pattern)
fileList = dir(fullfile(rootDir, pattern));
if isempty(fileList)
    error('Cannot find file with pattern: %s', pattern);
end

[~, order] = sort([fileList.datenum], 'descend');
filePath = fullfile(fileList(order(1)).folder, fileList(order(1)).name);
end

function dataset = read_one_dataset(filePath, tag)
raw = readmatrix(filePath);

if mod(size(raw, 1), 2) ~= 0
    error('Dataset %s has odd row count, which violates the csv format.', filePath);
end

temperatureRows = raw(1:2:end, :);
voltageRows = raw(2:2:end, :);

referenceTemperature = temperatureRows(1, :);
maxTemperatureDrift = max(abs(temperatureRows - referenceTemperature), [], 'all');
if maxTemperatureDrift > 1e-9
    error('Dataset %s has inconsistent temperature rows.', filePath);
end

if ~all(all(diff(voltageRows, 1, 2) < 0))
    error('Dataset %s contains non-monotone voltage rows.', filePath);
end

dataset.tag = tag;
dataset.filePath = filePath;
dataset.temperatureGrid = referenceTemperature;
dataset.voltageMatrix = voltageRows;
dataset.sampleCount = size(voltageRows, 1);
dataset.pointCount = size(voltageRows, 2);
dataset.minVoltage = min(voltageRows, [], 'all');
dataset.maxVoltage = max(voltageRows, [], 'all');
end

