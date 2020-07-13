function d = readTetFimBinary(binFile)

fid = fopen(binFile);
numPoints = fread(fid, 1, 'int32');
data = single(fread(fid, inf, 'single'));
n = floor(numel(data)/numPoints);
d = reshape(data(1:n*numPoints), numPoints, n);

end