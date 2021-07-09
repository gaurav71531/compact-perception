function genTestMatrices(nx, ny, nz)
bigX = randn(nx+ny+nz,nx+ny+nz);
bigX = bigX*bigX';

Sigma_X = bigX(1:nx,1:nx);
Sigma_Y = bigX(nx+1:nx+ny,nx+1:nx+ny);
Sigma_Z = bigX(nx+ny+1:nx+ny+nz,nx+ny+1:nx+ny+nz);

Sigma_YX = bigX(nx+1:nx+ny,1:nx);
Sigma_YZ = bigX(nx+1:nx+ny,nx+ny+1:nx+ny+nz);
Sigma_ZX = bigX(nx+ny+1:nx+ny+nz,1:nx);

Axy = Sigma_YX/Sigma_X;
Ayz = Sigma_YZ' / Sigma_Y;

Sigma_ey = Sigma_Y - Sigma_YX / Sigma_X * Sigma_YX';
Sigma_ez = Sigma_Z - Sigma_YZ' / Sigma_Y * Sigma_YZ;

save(sprintf('varUseComplete_%d_%d_%d.mat', nx, ny, nz), 'Sigma_X', 'Sigma_Y', 'Sigma_Z', ...
    'Sigma_YX', 'Sigma_YZ', 'Sigma_ZX', 'Axy', 'Ayz', 'Sigma_ey', 'Sigma_ez');

