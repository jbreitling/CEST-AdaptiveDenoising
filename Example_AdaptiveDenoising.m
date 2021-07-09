%% Tutorial for adaptive denoising
% Small script to show the denoising of CEST data using different criteria
% to determine the optimal number of components to retain.
%
% Reference:
%   Breitling J, Deshmane A, Goerke S, Korzowski A, Herz K, Ladd ME,
%   Scheffler K, Bachert P, Zaiss M. Adaptive denoising for chemical
%   exchange saturation transfer MR imaging. NMR Biomed. (2019);32:e4133.

%% Load example dataset:
%   Z_noise: Motion corrected, normalized and B0-corrected dataset acquired
%            at 7T in a healthy volunteer
%        dw: Frequency offsets for the Z-spectra
%   Segment: Manual segmentation of the parenchym

cd(fileparts(which(mfilename)))
load('ExampleDataset.mat')


%% Perform adaptive denoising
%   [Data_denoised,k] = PCA_denoising(Data,Segment,CriterionOrNumber,verbosity)
%                Data: Noisy Z-spectra (x,y,z,dw)
%             Segment: Segmentation of the 'useable' voxels, i.e. not
%                      outside the brain, not affected by fat
%   CriterionOrNumber: Fixes number of components or method name as string,
%                      i.e. 'Malinowski', 'Nelson', 'Median'
%           verbosity: Set to false to supress output
%
%       Data_denoised: Denoised Z-spectra (x,y,z,dw)
%                   k: Number of components used

Z_denoised = AdaptiveDenoising(Z_noise,Segment,'Median',true);

%% Plot example Z-spectra

voxel = [40,30,8]; % Set voxel location


figure, hold on
plot(dw,squeeze(Z_noise(voxel(1),voxel(2),voxel(3),:)),'.-','color','r');
plot(dw,squeeze(Z_denoised(voxel(1),voxel(2),voxel(3),:)),'.-','color','k');
xlim([-10 10]);
ylim([0 1]);
box on;
set(gca,'xdir','reverse');
xlabel('\Delta\omega [ppm]')
ylabel('Z')
legend({'noisy','denoised'})
hold off