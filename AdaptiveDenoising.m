function [Data_denoised,k] = AdaptiveDenoising(Data,Segment,CriterionOrNumber,verbosity)
% ADAPTIVEDENOISING. Denoising of CEST data using principal component
% analysis and a data-driven determination of the number of components
%
% INPUT:
%                Data: Noisy Z-spectra (x,y,z,dw)
%             Segment: Segmentation of the 'useable' voxels, i.e. not
%                      outside the brain, not affected by fat. By default
%                      the entire dataset is used.
%   CriterionOrNumber: Fixes number of components or method name as string,
%                      i.e. 'Malinowski', 'Nelson', 'Median'. By default
%                      the 'Median' criterion is used.
%           verbosity: Set to false to supress output
%
% OUTPUT:
%       Data_denoised: Denoised Z-spectra (x,y,z,dw)
%                   k: Number of components used
%
% Reference:
%   Breitling J, Deshmane A, Goerke S, Korzowski A, Herz K, Ladd ME,
%   Scheffler K, Bachert P, Zaiss M. Adaptive denoising for chemical
%   exchange saturation transfer MR imaging. NMR Biomed. (2019);32:e4133.


    %%
    if nargin < 2
        Segment = true(size(Data,1:3));
    end
    if nargin < 3
        CriterionOrNumber = 'Median';
    end
    if nargin < 4
        verbosity = true;
    end

    %% Prepare data
    % Vectorize data such that X = [Z-spectra,Voxel]
    if isnumeric(Segment); Segment = logical(Segment); end
    segment_idx = find(Segment); % vectorization indices
    Z1 = permute(Data,[4 1 2 3]); % permutation before vectorization
    X = Z1(:,segment_idx)';  % vectorize only useful data
    
    % centering of data, 'substract mean Z-spectrum'
    meanX = mean(X,1);
    data =  bsxfun(@minus,X,meanX);
    
    % calculate covariance matrix
    CovMat = cov(data);
    
    % calculate eigenvalues and eigenvectors of covariance matrix
    % -> the principle components (PC)
    [PC,D] = eigs(CovMat,size(X,2));
    
    % sort eigenvalues in size (seems to be necessary for older MatLab)
    [d,ix] = sort(diag(D),'descend');
    D = diag(d);
    PC = PC(:,ix);
 
    %% Malinowskis factor (empirical) indicator function
    % Malinowski, E. R. Determination of the number of factors and the
    % experimental error in a data matrix. Anal. Chem. 49, 612–617 (1977).
  
    IND = nan(size(X,2)-1,1);
    for ii_IE = 1: size(X,2)-1
       IND(ii_IE) = (sum(diag(D(ii_IE:end, ii_IE:end))./(size(X,2)-ii_IE)))^0.5 / (size(X,2)-ii_IE)^2;
    end
    [~,PCind] = min(IND);
    
    %% Nelson coefficient of determination
    % "Nelson, L.R. (2005). Some observations on the scree test, and
    % on coefficient alpha. Thai Journal of Educational Research and
    % Measurement (ISSN 1685-6740): 3(1), 1-17."
    
    R2 = nan(size(X,2)-2,1);
    for ii_CNG = 1:size(X,2)-2
        x = ii_CNG:size(X,2);
        x_new = [ones(numel(x),1) x'];
        Y = diag(D(ii_CNG:end,ii_CNG:end));
        b = x_new\Y;
        R2(ii_CNG) = (corr2(x_new*b,Y))^2;
    end
    
    % choose R^2 > 0.80 . Not clearly stated (but 0.80 mentioned)
    Scree_R2_ind = find(R2>0.80,1,'first');
    
    %% Median Noise Estimation
    % "Manjón, J.V. et al. (2015). MRI noise estimation and denoising using
    % non-local PCA. Medical Image Analysis 22 (2015) 35–47"
    
    MEDind = find(sqrt(diag(D)) >= 1.29*sqrt(median(d(sqrt(diag(D)) < 2*median(sqrt(diag(D)))))),1,'last');
    
    %% Output for recommended number of factors
    if verbosity
        fprintf('>> AdaptiveDenoising: Malinowskis factor (empirical) indicator function  proposes %d eigenvalues.\n', PCind)
        fprintf('>>                  : Nelson coefficient of determination  proposes %d eigenvalues.\n', Scree_R2_ind)
        fprintf('>>                  : Median method  proposes %d eigenvalues.\n', MEDind)
    end
    
    %% Set number of components to use
    
    if isnumeric(CriterionOrNumber)
        k = CriterionOrNumber;
        verbosity_string = ['>>                  : Use the ' num2str(k) ' eigenvalues.\n'];
    elseif strcmpi(CriterionOrNumber,'Nelson')
        verbosity_string = ['>>                  : Use the proposed ' num2str(Scree_R2_ind) ' eigenvalues of Nelson.\n'];
        k = Scree_R2_ind;
    elseif strcmpi(CriterionOrNumber,'Median')
        verbosity_string = ['>>                  : Use the proposed ' num2str(MEDind) ' eigenvalues of Median method.\n'];
        k = MEDind;
    elseif strcmpi(CriterionOrNumber,'None')
        verbosity_string =  '>>                  : No denoising applied.\n';
        k = numel(d);
    elseif strcmpi(CriterionOrNumber,'Malinowski')
        verbosity_string = ['>>                  : Use the proposed ' num2str(PCind) ' eigenvalues of Malinowski.\n'];
        k = PCind;
    end
    
    if verbosity
        fprintf(verbosity_string)
    end
    
    %% Projection onto the first k PCs
    Xd= (data*PC(:,1:k))*PC(:,1:k)';
    
    % add back the mean
    Xd = bsxfun(@plus,Xd,meanX);
    Xd = Xd';
    
    %% Sort data back
    Data_denoised = nan(size(Z1)); % allocate
    Data_denoised(:,segment_idx)=Xd;     % unvevtorize
    Data_denoised = permute(Data_denoised,[2 3 4 1]); % backpermute
    
end