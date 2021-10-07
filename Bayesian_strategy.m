function w = Bayesian_strategy(Phi, y)
%
% This code implements the Section 3.2 "Bayesian Strategy" 
%
% The optimization problem may be expressed as
%    minimize   alpha*|| y - Phi w ||_2^2 + || w ||_2^2
%
% This code is based on the Fast_RVM code available from http://www.miketipping.com/sparsebayes.htm
% [1] Tipping, Michael E, Sparse Bayesian learning and the relevance vector machine, Journal of machine learning research, 2001.
% [2] Tipping, Michael E and Faul, Anita C, Fast marginal likelihood maximisation for sparse Bayesian models, International workshop on artificial intelligence and statistics, 2003
% Please check the accompanying license and the license of [1] and [2] before using. 
%
% Inputs:
%    Phi:  transformed dictionary
%    y:  transformed sampled signal
% Outputs:
%    x:  solution of the above optimization problem
%
% Author: Hao Liang
% Last modified by: 21/10/04
%

% Parameter setting
gamma_0 =  std(y)^2/1e2;
eta = 1e-8;
maxIter = 1000;

% Find initial gamma
[~,m] = size(Phi);
PHIt = Phi'*y;
PHI2 = sum(Phi.^2)';
ratio = (PHIt.^2)./PHI2;
[maxr,index] = max(ratio);
gamma = PHI2(index)/(maxr-gamma_0);

% Compute initial mu, Sig, S, Q
phi = Phi(:,index);  % phi_i
Hessian = gamma + phi'*phi/gamma_0;
Sig = 1/Hessian;
mu = Sig*PHIt(index)/gamma_0;
left = Phi'*phi/gamma_0;
S = PHI2/gamma_0-Sig*left.^2;
Q = PHIt/gamma_0-Sig*PHIt(index)/gamma_0*left;

for count = 1:maxIter
    
    % Calculate si and qi
    s = S; q = Q;
    s(index) = gamma.*S(index)./(gamma-S(index));
    q(index) = gamma.*Q(index)./(gamma-S(index));
    theta = q.^2-s;
    
    % Choice the next alpha that maximizes marginal likelihood
    ml = -inf*ones(1,m);
    ig0 = find(theta>0);
    
    % Index for re-estimate
    [ire,~,which] = intersect(ig0,index);
    if ~isempty(ire)  % If not empty
        Gamma = s(ire).^2./theta(ire);
        delta = (gamma(which)-Gamma)./(Gamma.*gamma(which));
        ml(ire) = Q(ire).^2.*delta./(S(ire).*delta+1)-log(1+S(ire).*delta);
    end
    
    % Index for adding
    iad = setdiff(ig0,ire);
    if ~isempty(iad)
        ml(iad) = (Q(iad).^2-S(iad))./S(iad)+log(S(iad)./(Q(iad).^2));
    end
    is0 = setdiff([1:m],ig0);
    
    % Index for deleting
    [ide,~,which] = intersect(is0,index);
    if ~isempty(ide)
        ml(ide) = Q(ide).^2./(S(ide)-gamma(which))-log(1-S(ide)./gamma(which));
    end

    [ML(count),idx] = max(ml);
    
    % Stopping criteria
    if count > 2 && abs(ML(count)-ML(count-1)) < abs(ML(count)-ML(1))*eta
        break;
    end

    % Update gamma
    which = find(index==idx);
    if theta(idx) > 0
        if ~isempty(which)   % Re-estimate
            Gamma = s(idx)^2/theta(idx);
            Sigii = Sig(which,which); mui = mu(which); Sigi = Sig(:,which);
            delta = Gamma-gamma(which);
            ki = delta/(1+Sigii*delta);
            mu = mu-ki*mui*Sigi;
            Sig = Sig-ki*Sigi*Sigi';
            comm = Phi'*(phi*Sigi)/gamma_0;
            S = S + ki*comm.^2;
            Q = Q + ki*mui*comm;
            
            gamma(which) = Gamma;
        else    % Adding
            Gamma = s(idx)^2/theta(idx);
            phii = Phi(:,idx); Sigii = 1/(Gamma+S(idx)); mui = Sigii*Q(idx);
            comm1 = Sig*(phi'*phii)/gamma_0;
            ei = phii-phi*comm1;
            off = -Sigii*comm1;
            Sig = [Sig+Sigii*comm1*comm1', off; off', Sigii];
            mu = [mu-mui*comm1; mui];
            comm2 = Phi'*ei/gamma_0;
            S = S - Sigii*comm2.^2;
            Q = Q - mui*comm2;
            
            index = [index;idx];
            gamma = [gamma;Gamma];
            phi = [phi,phii];
        end
    else
        if ~isempty(which)   % Deleting
            Sigii = Sig(which,which); mui = mu(which); Sigi = Sig(:,which);
            Sig = Sig-Sigi*Sigi'/Sigii; Sig(:,which) = []; Sig(which,:) = [];
            mu  = mu-mui/Sigii*Sigi; mu(which) = [];
            comm = Phi'*(phi*Sigi)/gamma_0;
            S = S + comm.^2/Sigii;
            Q = Q + mui/Sigii*comm;
            
            index(which) = [];
            gamma(which) = [];
            phi(:,which) = [];
        end
    end
end

% Estimation of w
weights	= mu; used = index;
w = zeros(m,1);  w(used) = weights;

end
