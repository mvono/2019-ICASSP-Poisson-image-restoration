%-------------------------------------------------------------------------%
%                    SPLIT-AND-AUGMENTED GIBBS SAMPLER (SPA)              %
%                    APPLIED TO POISSONIAN IMAGE RESTORATION              %
%                      WITH TOTAL VARIATION REGULARIZATION                %
%-------------------------------------------------------------------------%
% File: SPA_analysis.m
% Author: M. VONO
% Created on: 12/05/2019
% Last modified : 12/05/2019

clearvars;
close all;
set(0,'DefaultFigureWindowStyle','docked'); % to dock figures
addpath('../utils/'); % to use special functions
addpath('../src/'); % to launch SPA

%-------------------------------------------------------------------------%
% REF.                                                                    %
% M. VONO et al.,                                                         %
% "Bayesian image restoration under Poisson noise and log-concave prior", %
% in Proc. IEEE Int. Conf. Acoust., Speech, and Signal Processing         %
% (ICASSP), Brighton, U.K., May 2019                                      %
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Pre-processing.

    % Load image
i = 8;
type = 'micro';
load(['image_' type '_',num2str(i)]);
im = double(eval(['image_' type '_',num2str(i)]));
N = size(im,1);
M = 30;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Launch SPA MCMC algorithm.

    % Parameters.
beta = 1; % regul. param. on TV prior (between 0.1 and 1)
rho = 1; % regul. param. associated to the splitting scheme
N_MC = 1e5; % total nb. of MCMC samples
N_bi = 5e4; % nb. of burn-in iter.
n_rep = 1;
    
    % MCMC sampling
for i = 1:n_rep

    rng(i);

    for m = 1:length(M)
        
        im = M(m) * im ./ max(im(:)); % scaling of the noisy/blurred image
        B = fspecial('gaussian',8,1); % variance = 1
        [FB,FBC,F2B,Bx] = HXconv(im,B,'Hx');
        obs = poissrnd(Bx);
        obs(isnan(obs)) = 0;

        [X_MC,time_SPA] = ...
        SPA_analysis(rho,beta,N,k,N_MC,N_bi,FBC,F2B,FB,obs,M(m))
        
        MAE_MMSE_SPA(i,m) = sum(sum(abs(im-mean(X_MC,3))))/N^2; 
        disp(['MAE_MMSE = ' num2str(MAE_MMSE_SPA(i,m))]);
        MAE_MMSE_norm_SPA(i,m) = MAE_MMSE_SPA(i,m)/M(m);
        disp(['MAE_norm_MMSE = ' num2str(MAE_MMSE_norm_SPA(i,m))]);
        X_MMSE_SPA(i,m,:,:) = mean(X_MC,3);
        figure(); imagesc(mean(X_MC,3)); colormap(gray); colorbar;
        disp(['Iteration n° ' num2str(i) ' for M = ' num2str(M(m)) ' finished.']);

    end
end