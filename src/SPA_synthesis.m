%-------------------------------------------------------------------------%
%                  SPLIT-AND-AUGMENTED GIBBS SAMPLER                      %
%-------------------------------------------------------------------------%

function [X_MC,time_SPA] = SPA_synthesis(rho,beta,N,k,N_MC,N_bi,FBC,F2B,FB,obs,W,WT,M)

% Parameters
alpha = 0.1 * rho; % regul. param. associated to the data aug. scheme
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
    eta1 = z1 - u1 + rho * randn(N,N);
    eta2 = z2 - u2 + rho * randn(N,k);
    eta3 = z3 - u3 + rho * randn(N,N);
    eta = WT(real(ifft2(FBC .* fft2(eta1)))) + eta2 + WT(eta3);
    ratio = 1 ./ (1 + 1 ./ (F2B + 1));
    X = eta - WT(real(ifft2(ratio .* fft2(W(eta)))));
    if t > N_bi
        X_MC(:,:,t-N_bi) = W(X);
    end

    % Sampling the auxiliary variable z1
    u = randn(N,N);
    gradH1 = (1 / rho^2) * real(ifft2(fft2(z1-u1) - FB .* fft2(W(X))));
    proxH1 = 0.5 * (z1 - lambdaPMYULA ...
            + sqrt((z1 - lambdaPMYULA).^2 + 4 * lambdaPMYULA * obs));
    z1 = (1 - gammaPMYULA / lambdaPMYULA) * z1 ...
            - gammaPMYULA * gradH1 ...
            + (gammaPMYULA / lambdaPMYULA) * proxH1 ...
            + sqrt(2 * gammaPMYULA) * u;

    % Sampling the auxiliary variable z2
    u = randn(N,k);
    gradH2 = (1 / rho^2) * (z2 - X - u2);
    proxH2 = sign(z2) .* max(abs(z2) - beta * lambdaPMYULA, zeros(N,k));
    z2 = (1 - gammaPMYULA / lambdaPMYULA) * z2 ...
            - gammaPMYULA * gradH2 ...
            + (gammaPMYULA / lambdaPMYULA) * proxH2 ...
            + sqrt(2 * gammaPMYULA) * u;

    % Sampling the auxiliary variable z3
    u = randn(N,N);
    gradH3 = (1 / rho^2) * (z3 - W(X) - u3);
    proxH3 = max(z3,0);
    z3 = (1 - gammaPMYULA / lambdaPMYULA) * z3 ...
            - gammaPMYULA * gradH3 ...
            + (gammaPMYULA / lambdaPMYULA) * proxH3 ...
            + sqrt(2 * gammaPMYULA) * u;

    % Sampling the auxiliary variable u1
    u1 = (alpha^2 / (alpha^2 + rho^2)) * ...
            (z1 - real(ifft2(FB .* fft2(W(X))))) ...
            + alpha * rho /sqrt(alpha^2 + rho^2) * randn(N,N);

    % Sampling the auxiliary variable u2
    u2 = (alpha^2 / (alpha^2 + rho^2)) * (z2 - X) ...
            + alpha * rho /sqrt(alpha^2 + rho^2) * randn(N,k);

    % Sampling the auxiliary variable u3
    u3 = (alpha^2 / (alpha^2 + rho^2)) * (z3 - W(X)) ...
            + alpha * rho /sqrt(alpha^2 + rho^2) * randn(N,N);

    % Update of the waitbar
    waitbar(t/N_MC);

end

close(h);
time_SPA = toc;

disp(['SPA sampling finished. Total run time = ' num2str(time_SPA) '.']);