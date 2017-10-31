%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variational Inference for Gaussian Mixture Model
% Jie Yang 2015
%
% Note: All parameters named by My_** are
%       parameters to be iterated/optimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GMM_VI()

load('data.mat')

n_obs=size(X,2);
n_dim=size(X,1);

n_iter=100;
K=[2 4 10 25];

c=10;
a0=n_dim;

X_mean=mean(X,2);
X_cov=cov(X');

n_rng=2;
figure;
for i_K=1:length(K)
    
    My_phi = zeros(n_obs,K(i_K),n_iter+1);
    My_n = zeros(K(i_K),n_iter+1);
    My_mu = zeros(n_dim,K(i_K),n_iter+1);
    My_covA = zeros(n_dim,n_dim,K(i_K),n_iter+1);
    My_alpha = zeros(K(i_K),n_iter+1);
    My_a = zeros(K(i_K),n_iter+1);
    My_B = zeros(n_dim,n_dim,K(i_K),n_iter+1);
    My_t1 = zeros(K(i_K),n_iter+1);
    My_t2 = zeros(n_obs,K(i_K),n_iter+1);
    My_t3 = zeros(K(i_K),n_iter+1);
    My_t4 = zeros(K(i_K),n_iter+1);
    
    P_x=zeros(n_iter+1,1);
    P_c=zeros(n_iter+1,1);
    P_pi=zeros(n_iter+1,1);
    P_lambda=zeros(n_iter+1,1);
    P_mu=zeros(n_iter+1,1);
    P_pos=zeros(n_iter+1,1);
    Q=zeros(n_iter+1,1);
    PQ=zeros(n_iter+1,1);
    
    %-------------------------------------
    % initialize My_parameters in some way
    %-------------------------------------
    
    rng(n_rng); 
    rand_mean=round(abs(max(max(X)))/2);
    My_mu(:,:,1)=X_mean*ones(1,K(i_K))+rand_mean*(rand(n_dim,K(i_K))-0.5);
    for ii_K=1:K(i_K)
        My_covA(:,:,ii_K,1)=X_cov;
    end;
    My_alpha(:,1)=ones(K(i_K),1);
    My_a(:,1)=a0*ones(K(i_K),1);
    for ii_K=1:K(i_K)
        My_B(:,:,ii_K,1)=n_dim/10*X_cov;
    end;
    
    for i_iter=2:n_iter+1
        
        %-------------------------------------
        % iteration of My_parameters
        %-------------------------------------
        
        % My_phi
        for ii_K=1:K(i_K)
            My_t1(ii_K,i_iter)=sum(psi((1-(1:n_dim)+My_a(ii_K,i_iter-1))/2))...
                -2*sum(log(diag(chol(My_B(:,:,ii_K,i_iter-1)))));
            My_t2(:,ii_K,i_iter)=diag((X-My_mu(:,ii_K,i_iter-1)*ones(1,n_obs))'...
                *(My_a(ii_K,i_iter-1)*inv(My_B(:,:,ii_K,i_iter-1)))*(X-My_mu(:,ii_K,i_iter-1)*ones(1,n_obs)));
            My_t3(ii_K,i_iter)=trace(My_a(ii_K,i_iter-1)*inv(My_B(:,:,ii_K,i_iter-1))...
                *My_covA(:,:,ii_K,i_iter-1));
            My_t4(ii_K,i_iter)=psi(My_alpha(ii_K,i_iter-1))-psi(sum(My_alpha(:,i_iter-1)));
            My_phi(:,ii_K,i_iter)=exp(0.5*My_t1(ii_K,i_iter)...
                -0.5*My_t2(:,ii_K,i_iter)-0.5*My_t3(ii_K,i_iter)+My_t4(ii_K,i_iter));
        end;
        My_phi(:,:,i_iter)=My_phi(:,:,i_iter)./(sum(My_phi(:,:,i_iter),2)*ones(1,K(i_K)));
        
        % My_n
        My_n(:,i_iter)=squeeze(sum(My_phi(:,:,i_iter)));
        
        % My_covA
        for ii_K=1:K(i_K)
            My_covA(:,:,ii_K,i_iter)=inv(1/c*eye(n_dim)+My_n(ii_K,i_iter)*...
                My_a(ii_K,i_iter-1)*inv(My_B(:,:,ii_K,i_iter-1)));
        end;
        
        % My_mu
        for ii_K=1:K(i_K)
            My_mu(:,ii_K,i_iter)=My_covA(:,:,ii_K,i_iter)*(My_a(ii_K,i_iter-1)*...
                inv(My_B(:,:,ii_K,i_iter-1))*sum(ones(n_dim,1)*My_phi(:,ii_K,i_iter)'.*X,2));
        end;
        
        % My_alpha
        My_alpha(:,i_iter)=My_alpha(:,1)+My_n(:,i_iter);
        
        % My_a
        My_a(:,i_iter)=My_a(:,1)+My_n(:,i_iter);
        
        % My_B
        for ii_K=1:K(i_K)
            My_B(:,:,ii_K,i_iter)=((X-My_mu(:,ii_K,i_iter)*ones(1,n_obs)).*(ones(n_dim,1)*...
                sqrt(My_phi(:,ii_K,i_iter))'))*((X'-ones(n_obs,1)*My_mu(:,ii_K,i_iter)').*...
                (sqrt(My_phi(:,ii_K,i_iter))*ones(1,n_dim)))...
                +sum(My_phi(:,ii_K,i_iter))*My_covA(:,:,ii_K,i_iter)+My_B(:,:,ii_K,1);
        end;
        
        %-------------------------------------
        % objective function (log)
        %-------------------------------------
        
        % log P(x|c,mu,lambda)
        P_x_temp1=zeros(n_obs,1);
        P_x_temp2=zeros(K(i_K),1);
        Det_B=zeros(K(i_K),1);
        for ii_K=1:K(i_K)
            Det_B(ii_K)=2*sum(log(diag(chol(My_B(:,:,ii_K,i_iter)))));
            P_x_temp1=P_x_temp1-...
                0.5*My_phi(:,ii_K,i_iter).*...
                (diag((X'-ones(n_obs,1)*My_mu(:,ii_K,i_iter)')*(My_a(ii_K,i_iter)*...
                inv(My_B(:,:,ii_K,i_iter)))*(X-My_mu(:,ii_K,i_iter)*ones(1,n_obs)))+...
                trace((My_a(ii_K,i_iter)*inv(My_B(:,:,ii_K,i_iter)))*My_covA(:,:,ii_K,i_iter)));
            P_x_temp2 (ii_K)=0.5*sum(My_n(ii_K,i_iter)*(n_dim*(n_dim-1)/4*log(pi)+n_dim*log(2)...
                +sum(psi(My_a(ii_K,i_iter)/2+(1-(1:n_dim))/2))-Det_B(ii_K)));
        end;
        P_x(i_iter)=-n_dim/2*n_obs*log(2*pi)+sum(P_x_temp1)+sum(P_x_temp2);
        
        % log P(c|pi)
        P_c_temp=zeros(n_obs,1);
        for ii_K=1:K(i_K)
            P_c_temp=P_c_temp+My_phi(:,ii_K,i_iter)*...
                (psi(My_alpha(ii_K,i_iter))-psi(sum(My_alpha(:,i_iter))));
        end;
        P_c(i_iter)=sum(P_c_temp);
        
        % log P(pi)
        P_pi(i_iter)=gammaln(sum(My_alpha(:,i_iter)))-sum(gammaln(My_alpha(:,i_iter)))+...
            sum((My_alpha(:,i_iter)-1).*(psi(My_alpha(:,i_iter))-psi(sum(My_alpha(:,i_iter)))));
        
        % log P(mu)
        P_mu_temp=zeros(K(i_K),1);
        for ii_K=1:K(i_K)
            P_mu_temp(ii_K)=sum(diag(My_covA(:,:,ii_K,i_iter)...
                +My_mu(:,ii_K,i_iter)*My_mu(:,ii_K,i_iter)'));
        end;
        P_mu(i_iter)=-K(i_K)*(log(2*pi)+0.5*log(c))-1/2/c*sum(P_mu_temp);
        
        % log P(lambda)
        psi_d=zeros(K(i_K),1);
        tr_B=zeros(K(i_K),1);
        for ii_K=1:K(i_K)
            psi_d(ii_K)=n_dim*(n_dim-1)/4*log(pi)+sum(psi(My_a(ii_K,i_iter)/2+(1-(1:n_dim))/2));
            tr_B(ii_K)=sum(trace(n_dim/10*X_cov*My_a(ii_K,i_iter)*inv(My_B(:,:,ii_K,i_iter))));
        end;
        gamma_d0=n_dim*(n_dim-1)/4*log(pi)+sum(gammaln(a0/2+(1-(1:n_dim))/2));
        P_lambda(i_iter)=(a0-n_dim-1)/2*sum(psi_d+n_dim*log(2)-Det_B)-0.5*sum(tr_B)-...
            a0/2*(n_dim*log(2)+2*sum(log(diag(chol(n_dim/10*X_cov)))))*K(i_K)-gamma_d0*K(i_K);
        
        % log P_pos
        P_pos(i_iter)=P_x(i_iter)+P_c(i_iter)+P_pi(i_iter)+P_mu(i_iter)+P_lambda(i_iter);
        
        % log q
        Gamma_d=zeros(K(i_K),1);
        for ii_K=1:K(i_K)
            Gamma_d(ii_K)=n_dim*(n_dim-1)/4*log(pi)+...
                sum(gammaln(My_a(ii_K,i_iter)/2+(1-(1:n_dim))/2));
            Q(i_iter)=Q(i_iter)-(n_dim+1)/2*Det_B(ii_K)+n_dim*(n_dim+1)/2*log(2)...
                +Gamma_d(ii_K)-(My_a(ii_K,i_iter)-n_dim-1)/2*psi_d(ii_K)+My_a(ii_K,i_iter)*n_dim/2;
        end;
        Det_A=zeros(K(i_K),1);
        for ii_K=1:K(i_K)
            Det_A(ii_K)=2*sum(log(diag(chol(My_covA(:,:,ii_K,i_iter)))));
            Q(i_iter)=Q(i_iter)+n_dim/2*(1+log(2*pi))+0.5*Det_A(ii_K);
        end;
        Q(i_iter)=Q(i_iter)+sum(-sum(My_phi(:,:,i_iter).*log(My_phi(:,:,i_iter)),2));
        Q(i_iter)=Q(i_iter)-gammaln(sum(My_alpha(:,i_iter)))+sum(gammaln(My_alpha(:,i_iter)))+...
            (sum(My_alpha(:,i_iter))-K(i_K))*psi(sum(My_alpha(:,i_iter)))...
            -sum((My_alpha(:,i_iter)-1).*psi(My_alpha(:,i_iter)));
        
        % pq
        PQ(i_iter)=P_pos(i_iter)+Q(i_iter);
        
    end;
    
    %-------------------------------------
    % Plot
    %-------------------------------------
    
    subplot(2,4,i_K)
    plot(1:n_iter,PQ(2:end))
    xlim([1 n_iter]);
    ylim([-1500 -1000]);
    xlabel('iteration')
    ylabel('objective function')
    title(['K=',num2str(K(i_K))]);
    set(gca,'FontSize',12,'FontWeight','b');
    
    subplot(2,4,i_K+4)
    [~,My_phi_M] = max(My_phi(:,:,n_iter),[],2);
    lk=1;
    for ii_K = 1:K(i_K)
        hold on;
        if ~isempty(find(My_phi_M==ii_K))
            scatter(X(1,find(My_phi_M==ii_K)),X(2,find(My_phi_M==ii_K)),'*')
            lgd{lk}=['Cluster ' num2str(lk)];
            lk=lk+1;
        end
    end;
    xlabel('Dimension 1')
    ylabel('Dimension 2')
    title(['K=',num2str(K(i_K))]);
    legend(lgd')
    set(gca,'FontSize',12,'FontWeight','b');
    axis square
end