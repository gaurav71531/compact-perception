%% generate test data
nx = 40;
ny = 30;
nz = 20;

% genTestMatrices(nx, ny, nz);
%% get rank of matrices

% lambda-gamma plot

beta = 100;
betaGrid = beta;
lambdaGrid = 0.01:0.1:2;
gammaGrid = 0.01:0.01:45;
% gammaGrid = 0.01:0.1:20;
gammaGrid = fliplr(gammaGrid);
[gridB, gridL, gridG] = meshgrid(betaGrid, lambdaGrid, gammaGrid);

fName = sprintf('varUseComplete_%d_%d_%d.mat', nx, ny, nz);
method = 'CP';


param = struct('paramPair', [gridB(:), gridL(:), gridG(:)]);
out = infCPMain('fName', fName, 'method', method, 'param', param);


