function [estIF, estIA, estMode] = AVNCMD(g, fs, iniIF, beta, tol)
%
% This code implements the AVNCMD algorithm 
%
% This code is based on the VNCMD code available from http://cn.mathworks.com/matlabcentral/fileexchange/64292 from the following paper
% [1] Chen S, Dong X, Peng Z, et al. Nonlinear chirp mode decomposition: A variational method. IEEE Transactions on Signal Processing, 2017.
% Please check the accompanying license and the license of [1] before using. 
%
%
% Inputs:
%    g:   sampled time-domain signal
%    fs:  sample frequency
%    iniIF:   initial IFs
%    beta:    filter parameter 
%    tol:     tolerance of convergence criterion
% Outputs:
%    estIF:     estimated IFs at each iteration 
%    estIA:     estimated group delays at each iteration 
%    estMode:   estimated modes at each iteration 
%
% Authors: Hao Liang (haoliang@stu.xmu.edu.cn) and Xiaotong Tu (xttu@xmu.edu.cn)
% Last modified by: 21/10/03
%

% Parameter setting
[K, N] = size(iniIF);  % K denotes the number of the modes, N the number of the samples
t = (0:N-1)/fs;   % time variables

% Construct the second-order difference matrix H
e = ones(N,1); e2 = -2*e;
H = spdiags([e e2 e], 0:2, N-2, N);  
HtH = H'*H;

% Form the 2K block second-order difference matrix D
D1 = blkdiag(H, H); tempm = repmat('D1,', 1, K);
D = eval(sprintf('blkdiag(%s)',tempm(1:end-1)));  

% Initialization
sinm = zeros(K,N); cosm = zeros(K,N);
uk = zeros(K,N); vk = zeros(K,N);   % the two demodulated quadrature signals
iternum = 300;  % the maximum allowable iterations
IFsetiter = zeros(K,N,iternum+1); IFsetiter(:,:,1) = iniIF;     % IF record 
Modeset_iter = zeros(K,N,iternum+1); Modeset_iter(:,:,1) = 0;   % Mode record

% Iniitialize the dictionary matrix A
A = zeros(N, 2*K*N);
for i = 1:K
    sinm(i,:) = sin(2*pi*(cumtrapz(t,iniIF(i,:))));
    cosm(i,:) = cos(2*pi*(cumtrapz(t,iniIF(i,:))));
    Sk = spdiags(sinm(i,:)', 0, N, N);
    Ck = spdiags(cosm(i,:)', 0, N, N);
    Ak = [Ck, Sk];
    A(1:N, 2*(i-1)*N+1:2*i*N) = Ak;
end

% Start iterations
iter = 1;        % iteration counter
sDif = tol + 1;  % tolerance counter

while ( sDif > tol &&  iter <= iternum )
    
    % Gradually increase the filter parameter during the iterations
    beta_thin = 10^(iter/36-10); 
    if beta_thin > beta
        beta_thin = beta;
    end
    
    % Estimating the Nonlinear Chirp Signal
    x = Estimate_NCS(A, D, g(:), K, N);  x = x';
    
    % This loop implements the Section 3.3 "Data-driven Implementation"
    for i = 1:K
        
        % Extract the two demodulated signals, uk and vk
        temp_uk = x(1,2*(i-1)*N+1:(2*i-1)*N); temp_vk = x(1,(2*i-1)*N+1:2*i*N);
        uk(i,:) = temp_uk; vk(i,:) = temp_vk;
        
        % Employing the arctangent demodulation technique to update the IF increment
        ukdif = Differ(vk(i,:),1/fs); vkdif = Differ(uk(i,:),1/fs);
        deltaIF = (uk(i,:).*ukdif - vk(i,:).*vkdif)./(uk(i,:).^2 + vk(i,:).^2)/2/pi;
        
        % IF increment may be corrected by a low-pass fliter
        deltaIF = (2/beta_thin*HtH + speye(N))\deltaIF';
        
        % Update the IF
        iniIF(i,:) = abs(iniIF(i,:) - deltaIF');
        
        % Update cos and sin functions 
        sinm(i,:) = sin(2*pi*(cumtrapz(t,iniIF(i,:))));
        cosm(i,:) = cos(2*pi*(cumtrapz(t,iniIF(i,:))));
        Sk = spdiags(sinm(i,:)', 0, N, N);
        Ck = spdiags(cosm(i,:)', 0, N, N);
        Ak = [Ck, Sk]; xk = [uk(i,:) vk(i,:)].';
        
        % Estimating the k-th mode
        Modeset_iter(i,:,iter+1) = Ak*xk;
        
        % Construct the dictionary for the next iteration
        A(1:N, 2*(i-1)*N+1:2*i*N) = Ak;
        
    end
    
    % Record the IFs
    IFsetiter(:,:,iter+1) = iniIF;
    
    % Stopping criteria
    sDif = 0;
    for i = 1:K
        sDif = sDif + (norm(Modeset_iter(i,:,iter+1) - Modeset_iter(i,:,iter))/norm(Modeset_iter(i,:,iter))).^2;
    end
    iter
    sDif
    iter = iter + 1;
    
end

% Demodulated results
estIF = IFsetiter(:,:,1:iter);         % Estimated IF
estMode = Modeset_iter(:,:,1:iter);    % Estimated Modes
estIA = sqrt(uk.^2 + vk.^2);           % Estimated IAs

end

