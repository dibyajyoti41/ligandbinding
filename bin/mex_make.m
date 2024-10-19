% build file

% compile the cpp file
cc = mex.getCompilerConfigurations('C++');
devel = cc.Manufacturer;
files = 'OVERLAPPING_DTW_MEX.cpp';
flags = 'CXXFLAGS="-O3 -fPIC -march=native -funroll-loops -fwrapv -ffp-contract=fast"';

 mex('-v', '-R2018a', files, flags);
 