function out = infCPMain(varargin)


% matObj = matfile('varUseNewComplete_100_70_30.mat');

param = [];
out = struct();
for i = 1:2:nargin
    switch varargin{i}
        case 'fName'
            fName = varargin{i+1};
            
        case 'method'
            method = varargin{i+1};
            
        case 'param'
            param = varargin{i+1};
    end
end
matObj = matfile(fName);
epsilon = 1e-4;
Sigma_X = matObj.Sigma_X;
Axy = matObj.Axy;
Ayz = matObj.Ayz;
Sigma_ey = matObj.Sigma_ey;
Sigma_ez = matObj.Sigma_ez;


Sigma_Y = Axy * Sigma_X * Axy' + Sigma_ey;
Sigma_Z = Ayz * Sigma_Y * Ayz' + Sigma_ez;
Sigma_ZgX = Ayz * Sigma_ey * Ayz' + Sigma_ez;

Sigma_XgY = Sigma_X - Sigma_X * Axy' /Sigma_Y * Axy * Sigma_X;
Sigma_YgZ = Sigma_Y - Sigma_Y * Ayz' /Sigma_Z * Ayz * Sigma_Y;

Sigma_ZX = Ayz * Axy * Sigma_X;
Sigma_XgZ = Sigma_X - Sigma_ZX' / Sigma_Z * Sigma_ZX;

niter = 100;

switch method       
    case 'CP'

        if isempty(param)      
            beta = 300;
            betaGrid = beta;
            lambdaGrid = 0.01:0.1:2;
            gammaGrid = 0.01:0.1:20;

            gammaGrid = fliplr(gammaGrid);
            [gridB, gridL, gridG] = meshgrid(betaGrid, lambdaGrid, gammaGrid);
            paramPair = [gridB(:), gridL(:), gridG(:)];
        else
            paramPair = param.paramPair;
        end
        
            

        MI_XB2_grid = zeros(size(paramPair,1),2+2+1+3);   
        rank_grid = zeros(size(paramPair,1),2);
        functional_grid = cell(size(paramPair,1),1);
        parfor pInd = 1:size(paramPair,1)
                betaUse = paramPair(pInd,1);
                gammaUse = paramPair(pInd,3);
                lambdaUse = paramPair(pInd,2);
                [phiMat, deltaMat, Sigma_si1, Sigma_si2, functionalIter, Isfeasible] = getDeltaPhiMat(Sigma_X, Axy, Ayz, Sigma_ey, Sigma_ez, betaUse, gammaUse, lambdaUse, niter, epsilon);
                functional_grid{pInd} = functionalIter;
                if Isfeasible
                    Sigma_B1 = phiMat * Sigma_X * phiMat' + Sigma_si1;
                    Sigma_B2 = deltaMat *Sigma_B1 * deltaMat' + Sigma_si2;

                    Sigma_B1gZ = phiMat * Sigma_XgZ * phiMat' + Sigma_si1;
                    Sigma_B1gY = phiMat * Sigma_XgY * phiMat' + Sigma_si1;
                    Sigma_B2gZ = deltaMat * Sigma_B1gZ * deltaMat' + Sigma_si2;
                    Sigma_B2gY = deltaMat * Sigma_B1gY * deltaMat' + Sigma_si2;
                    Sigma_B2gX = deltaMat * Sigma_si1 * deltaMat' + Sigma_si2;

                    MI_XB2Temp = getMutInf(Sigma_B2, Sigma_B2gX);
                    MI_B2ZTemp = getMutInf(Sigma_B2, Sigma_B2gZ);
                    MI_YB2Temp = getMutInf(Sigma_B2, Sigma_B2gY);
                    MI_XB2_grid(pInd,:) = [MI_B2ZTemp, MI_XB2Temp, MI_YB2Temp, rank(phiMat), rank(deltaMat), paramPair(pInd,:)];
                    rank_grid(pInd,:) = [rank(phiMat), rank(deltaMat)];
                end
        end
%         assignin('base', 'MI_XB2_grid', MI_XB2_grid);
%         assignin('base', 'rank_grid', rank_grid);
%         assignin('base', 'functional_grid', functional_grid);
        out.MI_XB2_grid = MI_XB2_grid;
        out.rank_grid = rank_grid;
        out.functional_grid = functional_grid;
end



function [phiMatFin, deltaMatFin, Sigma_si1Fin, Sigma_si2Fin,...
    funcIter, Isfeasible] = getDeltaPhiMat(Sigma_X, Axy, Ayz, Sigma_ey,...
                            Sigma_ez, beta, gamma, lambda, niter, epsilon)

nx = size(Sigma_X,1);

Sigma_Y = Axy * Sigma_X * Axy' + Sigma_ey;
Sigma_Z = Ayz * Sigma_Y * Ayz' + Sigma_ez;

Sigma_XgY = Sigma_X - Sigma_X * (Axy' /Sigma_Y) * Axy * Sigma_X;

Sigma_ZX = Ayz * Axy * Sigma_X;
Sigma_XgZ = Sigma_X - (Sigma_ZX' / Sigma_Z) * Sigma_ZX;

Sigma_si1 = cell(1,niter);
Sigma_si2 = cell(1,niter);

Sigma_B1 = cell(1,niter);
Sigma_B2 = cell(1,niter);

Sigma_B1gY = cell(1,niter);
Sigma_B1gZ = cell(1,niter);
Sigma_B2gZ = cell(1,niter);

deltaMat = cell(1,niter);
phiMat = cell(1,niter);

phiMat{1} = eye(nx);
deltaMat{1} = eye(nx);

Sigma_si1{1} = eye(nx);
Sigma_si2{1} = eye(nx);

Sigma_B1{1} = phiMat{1} *Sigma_X * phiMat{1}' + Sigma_si1{1};
Sigma_B2{1} = deltaMat{1} *Sigma_B1{1} * deltaMat{1}' + Sigma_si2{1};

Sigma_B1gY{1} = phiMat{1} * Sigma_XgY * phiMat{1}' + Sigma_si1{1};
Sigma_B1gZ{1} = phiMat{1} * Sigma_XgZ * phiMat{1}' + Sigma_si1{1};
Sigma_B2gZ{1} = deltaMat{1} * Sigma_B1gZ{1} * deltaMat{1}' + Sigma_si2{1};


constMat = eye(nx) - Sigma_XgY/Sigma_X;
constMat1 = (Sigma_ZX' / Sigma_Z) * (Sigma_ZX / Sigma_X);

Isfeasible = 1;
funcIter = zeros(niter,1);

for iterInd = 2:niter
    
    Sigma_B1gYInv = getSVDBasedInverse(Sigma_B1gY{iterInd-1}, epsilon);
    Sigma_B1Inv = getSVDBasedInverse(Sigma_B1{iterInd-1}, epsilon);
    Sigma_si2Inv = getSVDBasedInverse(Sigma_si2{iterInd-1}, epsilon);
    
    Sigma_B2Inv = getSVDBasedInverse(Sigma_B2{iterInd-1}, epsilon);
    
    matTemp = beta * Sigma_B1gYInv - (beta-1) * Sigma_B1Inv ...
        + lambda * deltaMat{iterInd-1}' * Sigma_B2Inv * deltaMat{iterInd-1};
    [~,E] = eig(matTemp);
    E = diag(E);
    Eneg = E<0;
    EnegInd = find(abs(E(Eneg)) > 1e-5);
    
    if length(EnegInd)>0
%         disp('neg start');
        Isfeasible = 0;
        break;
    end

    Sigma_si1{iterInd} = getSVDBasedInverse(matTemp, epsilon);
    
    Sigma_B2gZInv = getSVDBasedInverse(Sigma_B2gZ{iterInd-1}, epsilon);
    matTemp1 = beta * Sigma_B1gYInv * phiMat{iterInd-1} * constMat ...
       + lambda * gamma * deltaMat{iterInd-1}' * Sigma_B2gZInv * deltaMat{iterInd-1} * phiMat{iterInd-1} * constMat1;
    phiMat{iterInd} = Sigma_si1{iterInd} * matTemp1;
    
    Sigma_B1{iterInd} = phiMat{iterInd} *Sigma_X * phiMat{iterInd}' + Sigma_si1{iterInd};
    Sigma_B1gY{iterInd} = phiMat{iterInd} * Sigma_XgY * phiMat{iterInd}' + Sigma_si1{iterInd};
    Sigma_B1gZ{iterInd} = phiMat{iterInd} * Sigma_XgZ * phiMat{iterInd}' + Sigma_si1{iterInd};
    
    Sigma_B2gZ{iterInd} = deltaMat{iterInd-1} * Sigma_B1gZ{iterInd} * deltaMat{iterInd-1}' + Sigma_si2{iterInd-1};
    Sigma_B2{iterInd} = deltaMat{iterInd-1} *Sigma_B1{iterInd} * deltaMat{iterInd-1}' + Sigma_si2{iterInd-1};
    
    Sigma_B1Inv = getSVDBasedInverse(Sigma_B1{iterInd}, epsilon);
    Sigma_B2Inv = getSVDBasedInverse(Sigma_B2{iterInd}, epsilon);
    Sigma_B2gZInv = getSVDBasedInverse(Sigma_B2gZ{iterInd}, epsilon);
    
    matTemp2 = gamma * Sigma_B2gZInv - (gamma-1) * Sigma_B2Inv;
    Sigma_si2{iterInd} = getSVDBasedInverse(matTemp2, epsilon);
    
    matTemp3 = eye(nx) - Sigma_B1gZ{iterInd} * Sigma_B1Inv;
    deltaMat{iterInd} = gamma * Sigma_si2{iterInd} * Sigma_B2gZInv * deltaMat{iterInd-1} * matTemp3;  
    
    Sigma_B2{iterInd} = deltaMat{iterInd} *Sigma_B1{iterInd} * deltaMat{iterInd}' + Sigma_si2{iterInd};   
    
    Sigma_B2gZ{iterInd} = deltaMat{iterInd} * Sigma_B1gZ{iterInd} * deltaMat{iterInd}' + Sigma_si2{iterInd};
    
    funcIter(iterInd) = getMutInf(Sigma_B1{iterInd}, Sigma_si1{iterInd}) ...
        -beta*getMutInf(Sigma_B1{iterInd}, Sigma_B1gY{iterInd}) ...
        +lambda*getMutInf(Sigma_B2{iterInd}, Sigma_si2{iterInd}) ...
        -lambda*gamma*getMutInf(Sigma_B2{iterInd}, Sigma_B2gZ{iterInd});
     
end

phiMatFin = phiMat{end};
deltaMatFin = deltaMat{end};
Sigma_si1Fin = Sigma_si1{end};
Sigma_si2Fin = Sigma_si2{end};


function out = getSVDBasedInverse(X, epsilon)

n = size(X,1);
% if issymmetric(X)
%     [U,E] = eig(X);
%     E = diag(E);
%     Eneg = E<0;
%     EnegInd = find(abs(E(Eneg)) > 1e-5);
% %     EnegInd = Eneg(abs(E(Eneg))> 1e-5);
%     if sum(EnegInd)>0
%         disp('gg');
%     end
%     E(EnegInd) = 0;
%     EposInd = E>0;
%     EInv = E;
%     EInv(EposInd) = 1./EInv(EposInd);
%     out = U * diag(EInv) * U';
%         
%     
% else
    try
% %         [U,S,V] = svd(X);
% %         SUse = diag(S);
% % %     rUse = sum(SUse > epsilon);
% %         rUse = n;
% %         SUseInv = [1./SUse(1:rUse);zeros(n-rUse,1)];
% %         out = V * diag(SUseInv) * U';
        out = inv(X);
    catch
        disp('gg')
        out = 0;
    end
    
% end
% out = inv(X);



function out = getMutInf(Sigma_X, Sigma_XgY)

% out = log(det(Sigma_X)) - log(det(Sigma_XgY));
ee1 = eig(Sigma_X);
ee2 = eig(Sigma_XgY);
out = sum(log(ee1)) - sum(log(ee2));
if out < 0
%     fprintf('negative mut inf\n');
    out = Inf;
end