%32 bit window
lapacklib = fullfile(matlabroot, ...
  'extern', 'lib', 'win32', 'microsoft', 'libmwlapack.lib');
blaslib = fullfile(matlabroot, ...
  'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');
mex('-v', '-largeArrayDims', 'ndfun.c', blaslib, lapacklib);

%64 bit window
lapacklib = fullfile(matlabroot, ...
  'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
blaslib = fullfile(matlabroot, ...
  'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
mex('-v', '-largeArrayDims', 'ndfun.c', blaslib, lapacklib);

%64 bit linux
mex -v -largeArrayDims ndfun.c -lmwlapack -lmwblas