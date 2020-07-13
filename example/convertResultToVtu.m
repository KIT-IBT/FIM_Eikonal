addpath('../../vtkToolbox/MATLAB')

box = vtkRead('box.vtu');
d = readTetFimBinary('result.bin');
box.pointData.seed1 = d(:,1);
box.pointData.seed2 = d(:,2);
box.pointData.seed3 = d(:,3);
box.pointData.seed4 = d(:,4);
vtkWrite(box, 'result.vtu');