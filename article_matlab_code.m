% After extracting timeseries from brainstorm.

%% Averaged functional connectivity in beta band
srate=1000;
fmin=13;
fmax=30;
FrequencyBand=[fmin fmax];
filterorder=4/FrequencyBand(1);
v=reshape(Value,68,57,40000); % 68: ROIs nbre 
                              % 57: nbre of participants and epochs (19*3)
                              % 40000: epoch duration 40s since sampling frequency = 1000Hz
mean_v=squeeze(mean(v,2)); % averaged timeserie over all participants
nchannels = size(mean_v,1);
b1 = fir1(floor(filterorder*srate),[FrequencyBand(1) FrequencyBand(2)]/(srate/2)); % filtering with FIR
for j = 1:nchannels
    FilteredSignal(j,:) = filtfilt(b1,1,double(mean_v(j,:)));
end
% PLV
plv=smoothed_PLV_v2(FilteredSignal,srate,40); % 40: Window of analysis (Duration in time)
plv_beta=squeeze(plv)+squeeze(plv)';
thresh_plv_beta=ThreshMat(plv,10); % thresholding the matrix
thresh_plv_beta=thresh_plv_beta+thresh_plv_beta';
% AEC
[env,~] = osl_envelope(FilteredSignal); % can be downloaded from:https://github.com/OHBA-analysis/osl-core
[AEC_beta, ~] =corr(env','type','Pearson');
matrix=triu(AEC_beta); % thresholding the matrix
new_matrix=matrix.*~eye(size(matrix));
thresh_AEC_beta=ThreshMat(new_matrix,10);
thresh_AEC_beta=thresh_AEC_beta+thresh_AEC_beta';
% PLI
pli_beta=Phase_lag_index(FilteredSignal);
thresh_pli_beta=ThreshMat(triu(pli_beta),10); % thresholding the matrix
thresh_pli_beta=thresh_pli_beta+thresh_pli_beta';
% PLV after applying orthogonalisation correction
nodeData = ROInets.remove_source_leakage(FilteredSignal, 'symmetric'); % orthogonalisation can be downloaded from: https://github.com/OHBA-analysis/MEG-ROI-nets
plv=smoothed_PLV_v2(nodeData,srate,40);
plv_orth_beta=squeeze(plv)+squeeze(plv)';
thresh_orth_plv_beta=ThreshMat(plv,10); % thresholding the matrix
thresh_orth_plv_beta=thresh_orth_plv_beta+thresh_orth_plv_beta';
% AEC after applying orthogonalisation correction
[env,~] = osl_envelope(nodeData);
[AEC_orth_beta, ~] =corr(env','type','Pearson');
matrix=triu(AEC_orth_beta); % thresholding the matrix
new_matrix=matrix.*~eye(size(matrix));
thresh_orth_AEC_beta=ThreshMat(new_matrix,10);
thresh_orth_AEC_beta=thresh_orth_AEC_beta+thresh_orth_AEC_beta';


%% Functional connectivity by epoch
plv_beta_epochs=zeros(57,68,68);
thresh_plv_beta_epochs=zeros(57,68,68);
AEC_beta_epochs=zeros(57,68,68);
thresh_AEC_beta_epochs=zeros(57,68,68);
pli_beta_epochs=zeros(57,68,68);
thresh_pli_beta_epochs=zeros(57,68,68);
plv_orth_beta_epochs=zeros(57,68,68);
thresh_orth_plv_beta_epochs=zeros(57,68,68);
AEC_orth_beta_epochs=zeros(57,68,68);
thresh_orth_AEC_beta_epochs=zeros(57,68,68);

for i=1:57
    for j = 1:nchannels
        FilteredSignal(j,:) = filtfilt(b1,1,double(squeeze(v(j,i,:))));
    end
    % PLV
    plv=smoothed_PLV_v2(FilteredSignal,srate,40);
    plv_beta_epochs(i,:,:)=squeeze(plv)+squeeze(plv)';
    thresh_plv=ThreshMat(squeeze(plv),10);
    thresh_plv_beta_epochs(i,:,:)=thresh_plv+thresh_plv';
    % AEC
    [env,~] = osl_envelope(FilteredSignal);
    [AEC_beta_epochs(i,:,:), ~] =corr(env','type','Pearson');
    matrix=triu(squeeze(AEC_beta(i,:,:))); % thresholding the matrix
    new_matrix=matrix.*~eye(size(matrix));
    thresh_AEC=ThreshMat(new_matrix,10);
    thresh_AEC_beta_epochs(i,:,:)=thresh_AEC+thresh_AEC';
    % PLI
    pli_beta_epochs(i,:,:)=Phase_lag_index(FilteredSignal);
    thresh_pli=ThreshMat(triu(squeeze(pli_beta_epochs(i,:,:))),10); % thresholding the matrix
    thresh_pli_beta_epochs(i,:,:)=thresh_pli+thresh_pli';
    % PLV after applying orthogonalisation correction
    nodeData = ROInets.remove_source_leakage(FilteredSignal, 'symmetric'); % orthogonalisation
    plv=smoothed_PLV_v2(nodeData,srate,40);
    plv_orth_beta_epochs(i,:,:)=squeeze(plv)+squeeze(plv)';
    thresh_orth_plv=ThreshMat(plv,10); % thresholding the matrix
    thresh_orth_plv_beta_epochs(i,:,:)=thresh_orth_plv+thresh_orth_plv';
    % AEC after applying orthogonalisation correction
    [env,~] = osl_envelope(nodeData);
    [AEC_orth_beta(i,:,:), ~] =corr(env','type','Pearson');
    matrix=triu(squeeze(AEC_orth_beta(i,:,:))); % thresholding the matrix
    new_matrix=matrix.*~eye(size(matrix));
    thresh_orth_AEC=ThreshMat(new_matrix,10);
    thresh_orth_AEC_beta_epochs(i,:,:)=thresh_orth_AEC+thresh_orth_AEC';   
end

%% reordering connectivity matrices (in this example I chose PLV)
ind=[27 28 43 44 7 8 23 24 9 10 35 36 65 66 13 14 61 62 17 18 31 32 67 68 1 2 63 64 59 60 15 16 51 52 57 58 55 56 5 6 37 38 41 42 39 40 25 26 29 30 11 12 49 50 33 34 45 46 53 54 3 4 47 48 21 22 19 20];
% ind is the new order of ROIs used in the article
new_mat=thresh_plv_beta(ind,ind); 
C=ones(68,1); % dividing the brain regions in groups according to which brain lobe they are part of
C(9:26)=C(9:26)*2;
C(27:34)=C(27:34)*3;
C(35:52)=C(35:52)*4;
C(53:58)=C(53:58)*5;
C(59:68)=C(59:68)*6;
[X,Y,INDSORT] = grid_communities(C);
imagesc(new_mat(INDSORT,INDSORT));
colormap(hot)
hold on;
plot(X,Y,'r','linewidth',6);

%% Participant correlation values between thresholded fMRI and EEG connectivity methods (in this example I chose PLV)
k=1;
for i=1:3:55 % since each participant has 3 epochs
mat(1,:,:)=squeeze(thresh_plv_beta_epochs(i,:,:)); % thresholded EEG functional connectivity matrix
mat(2,:,:)=squeeze(thresh_plv_beta_epochs(i+1,:,:)); 
mat(3,:,:)=squeeze(thresh_plv_beta_epochs(i+2,:,:)); 
corr_plv(k,1)=corr2(thresh_FC,squeeze(mean(mat,1))); % thresh_FC is the thresholded fMRI connectivity matrix
k=k+1;
end
% the box plots presented in the article are obtained using plotly 

%% EEG edges correlation with fMRI (in this example I chose PLV)
k=1;
for i=1:67
    for j=i+1:68
        edge_fMRI(k,1)=FC_fMRI(i,j); % FC_fMRI is the fMRI connectivity matrix without threshold
        edge_PLV(k,1)=plv_beta(i,j); % EEG functional connectivity matrix without threshold
        k=k+1;
    end
end
[ro, pval] =corr(edge_PLV,edge_fMRI,'type','Spearman');
figure,
plot (edge_PLV,edge_fMRI,'bo','MarkerFaceColor','b','MarkerSize',10),
grid ON
h=lsline;
set(h,'LineWidth',4,'color','r')