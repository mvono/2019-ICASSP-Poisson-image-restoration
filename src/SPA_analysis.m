%-------------------------------------------------------------------------%
%                  SPLIT-AND-AUGMENTED GIBBS SAMPLER                      %
%-------------------------------------------------------------------------%

function [X_MC,time_SPA] = SPA_analysis(rho,beta,N,k,N_MC,N_bi,FBC,F2B,FB,obs,W,WT,M)

% Parameters
alpha = rho; % regul. param. associated to the data aug. scheme
lambdaPMYULA = rho^2; % param. P-MYULA
gammaPMYULA = rho^2 / 4; % param. P-MYULA

% Initialization
z1 = randi([0 M],N,N);
z2 = randi([0 M],N,k);
z3 = randi([0 M],N,N);
u1 = randi([0 M],N,N);
u2 = randi([0 M],N,k);
u3 = randi([0 M],N,N);
        
h = waitbar(0,'SPA sampling in progress...');
tic;

for t = 1:N_MC

	% Sampling the variable of interest x
    cov = (rho^2) ./ (F2B  + 2);
    moy = (1 / rho^2) * cov .* (FBC .* fft2(z1-u1) ...
                  + fft2(z2 + z3 - u2 - u3)); 
    eps = sqrt(0.5) * (randn(N,N) + sqrt(-1)*randn(N,N));
    X = real(ifft2(moy + eps .* sqrt(cov)));
    if t > N_bi
        X_MC(:,:,t-N_bi) = X;
    end

    % Sampling the auxiliary variable z1
    u = randn(N,N);
    gradH1 = (1 / rho^2) * real(ifft2(fft2(z1-u1) - FB .* fft2(X)));
    proxH1 = 0.5 * (z1 - lambdaPMYULA ...
            + sqrt((z1 - lambdaPMYULA).^2 + 4 * lambdaPMYULA * obs));
    z1 = (1 - gammaPMYULA / lambdaPMYULA) * z1 ...
            - gammaPMYULA * gradH1 ...
            + (gammaPMYULA / lambdaPMYULA) * proxH1 ...
            + sqrt(2 * gammaPMYULA) * u;

    % Sampling the auxiliary variable z2
    u = randn(N,N);
    gradH2 = (1 / rho^2) * (z2 - X - u2);
    [proxH2,~] = chambolle_prox_TV_stop(z2, ...
                 'lambda', beta*lambdaPMYULA, ...
                 'maxiter', 10);
    z2 = (1 - gammaPMYULA / lambdaPMYULA) * z2 ...
            - gammaPMYULA * gradH2 ...
            + (gammaPMYULA / lambdaPMYULA) * proxH2 ...
            + sqrt(2 * gammaPMYULA) * u;

    % Sampling the auxiliary variable z3
    u = randn(N,N);
    gradH3 = (1 / rho^2) * (z3 - X - u3);
    proxH3 = max(z3,0);
    z3 = (1 - gammaPMYULA / lambdaPMYULA) * z3 ...
            - gammaPMYULA * gradH3 ...
            + (gammaPMYULA / lambdaPMYULA) * proxH3 ...
            + sqrt(2 * gammaPMYULA) * u;

    % Sampling the auxiliary variable u1
    u1 = (alpha^2 / (alpha^2 + rho^2)) * ...
            (z1 - real(ifft2(FB .* fft2(X)))) ...
            + alpha * rho /sqrt(alpha^2 + rho^2) * randn(N,N);

    % Sampling the auxiliary variable u2
    u2 = (alpha^2 / (alpha^2 + rho^2)) * (z2 - X) ...
          + alpha * rho /sqrt(alpha^2 + rho^2) * randn(N,N);

    % Sampling the auxiliary variable u3
    u3 = (alpha^2 / (alpha^2 + rho^2)) * (z3 - X) ...
          + alpha * rho /sqrt(alpha^2 + rho^2) * randn(N,N);

    % Update of the waitbar
    waitbar(t/N_MC);

end

time_SPA = toc;

disp(['SPA sampling finished. Total run time = ' num2str(time_SPA) '.']);