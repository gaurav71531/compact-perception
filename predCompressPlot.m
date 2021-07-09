clear
%%
nx_ny_nz_tuple = [50, 40, 30;...
                 100, 70, 30;...
                  100, 40, 100;...
                  150, 100, 50;...
                  150, 150, 150;...
                  300, 100, 50;...
                  ];
%%
for i = 1:size(nx_ny_nz_tuple, 1)
    nx = nx_ny_nz_tuple(i, 1);ny = nx_ny_nz_tuple(i, 2);nz = nx_ny_nz_tuple(i, 3);
    genTestMatrices(nx, ny, nz);
    fName = sprintf('varUseComplete_%d_%d_%d.mat', nx, ny, nz);
    fprintf('starting CP case for (nx, ny, nz) = (%d, %d, %d)\n', nx, ny, nz);
    method = 'CP';
    infCPMain('fName', fName, 'method', method);
    save(sprintf('CP_%d_%d_%d.mat',nx, ny, nz), 'MI_XB2_grid');
end

%% plot figures

nx_ny_nz_tuple = [50, 40, 30;...
                 100, 70, 30;...
                  100, 40, 100;...
                  150, 100, 50;...
                  150, 150, 150;...
                  300, 100, 50;...
                  ];             
for i = 1:size(nx_ny_nz_tuple,1)
    nx = nx_ny_nz_tuple(i, 1);ny = nx_ny_nz_tuple(i, 2);nz = nx_ny_nz_tuple(i, 3);
    load(sprintf('CP_%d_%d_%d.mat',nx, ny, nz));

    MI_Max = MI_XB2_grid(1,1);

    minThres = 1e-3;
    MI_MaxInd_XB2 = 1;
    MI_MinInd_XB2 = find(MI_XB2_grid(:,1)<1e-3,1,'first')-1; 

    figure;
    legendName = {};
    plot(MI_XB2_grid(MI_MaxInd_XB2:MI_MinInd_XB2,1),MI_XB2_grid(MI_MaxInd_XB2:MI_MinInd_XB2,2), 'b');
    set(0,'defaultTextInterpreter','latex');
    set(0, 'DefaultFigureColor', [1 1 1]);
    set(gca, 'fontsize', 18);
    xlabel('$I({\bf B}_{1};{\bf X}_{2})$');
    ylabel('$I({\bf X}_{0};{\bf B}_{0})$');
    title(sprintf('(%d, %d, %d)', nx, ny, nz));
    grid;
%     saveFileName = sprintf('CP_%d_%d_%d', nx, ny, nz);
%     export_fig(saveFileName, '-pdf');
end