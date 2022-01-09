function full_f_save = session_output_simulation(fol, t_offset, t_max, seed_start,varargin)
[wildcard,min_bytes,unstable_a_thresh,full,simFolder,unstable_timeout, stage2] = process_options(varargin,'wildcard','CVCrit','min_bytes',150,'unstable_a_thresh',1,'full',0,'simFolder','','unstable_timeout',[],'stage2',0);
%fol = 'fit18';
%t_offset = 25e4;
t_range = [1 t_max]+t_offset;
%seed_start = 99;

%% Load and sort all results
saveSubFolder = fullfile('exampleData','exampleResult',fol);
[R, fileNames] = loadAllFittingResults(saveSubFolder,'wildcard',wildcard,'verbose',1,'min_bytes',min_bytes);
R = [R.fit];
for i=1:length(R)
    n_ySp(i) = R(i).yCell.nSpikes;
    yLabels{i} =  R(i).yCell.label;
end
[~, sorted_ySp_inds] = sort(n_ySp,'descend');
R = R(sorted_ySp_inds);
yLabels=yLabels(sorted_ySp_inds);
fileNames=fileNames(sorted_ySp_inds);
nF = length(R);

%% Run all simulations
y_sim = nan(t_range(2)-t_range(1)+1,nF);
y_real = nan(t_range(2)-t_range(1)+1,nF);
for i=1:nF
    tstart = tic;
    switch wildcard
        case 'CVCrit'
            [y_sim_temp, y_unstable_temp, a_temp, u_temp, Pf_temp,n_temp] = R(i).simulateOutput(t_range,seed_start+i,'unstable_a_thresh',unstable_a_thresh,'full',full,'unstable_timeout',unstable_timeout,'stage2',stage2);
        otherwise
            [y_sim_temp, y_unstable_temp, a_temp, u_temp, Pf_temp] = R(i).simulateOutput_different_basis(t_range,seed_start+i,'unstable_a_thresh',unstable_a_thresh,'full',full,'unstable_timeout',unstable_timeout);
    end
    y_sim(:,i) = y_sim_temp;
    a(:,i) = a_temp;
    u(:,i) = u_temp;
    Pf(:,i) = Pf_temp;
    n(:,i) = n_temp;
    y_unstable(i)=y_unstable_temp;
    [x_real,y_real_temp]= R(i).getXY('sess_t',1);
    y_real(:,i) = y_real_temp(t_range(1):t_range(2));
    x_real=x_real(t_range(1):t_range(2),:);
    tend = toc(tstart);
    disp(['Done ' num2str(i) ' out of ' num2str(nF) ' in ' num2str(tend) 's'])
end

%% Save simulation results
if stage2
    stage2str='st2';
else
    stage2str='';
end
f_save = ['y_sim_' fol '_t' num2str(t_range(1)) 'thru' num2str(t_range(2)) '_seedStart' num2str(seed_start) '_nY' num2str(nF) '_' wildcard '_full' num2str(full) stage2str '_d' num2str(now) '.mat'];

if ~isempty(simFolder)
    saveSubFolder=[saveSubFolder filesep simFolder];
end
if ~exist(saveSubFolder,'dir')
    mkdir(saveSubFolder)
end
    full_f_save = [saveSubFolder filesep f_save];
save(full_f_save, 'y_sim', 'y_real','t_range','seed_start','nF','y_unstable','sorted_ySp_inds','yLabels','x_real','a','u','Pf','n')






