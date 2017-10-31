%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expectation Maximization for Gaussian Mixture Model
% Jie Yang 2015
%
% Note: All parameters named by My_** are
%       parameters to be iterated/optimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GMM_EM()

load('data.mat')
n_obs=size(X,2);
n_dim=size(X,1);

n_iter=100;
K=[2 4 8 10];

X_mean=mean(X,2);
X_cov=cov(X');

n_rng=0;
figure;

for i_K=1:length(K)
    
    My_pi=zeros(K(i_K),n_iter+1);
    My_mu=zeros(n_dim,K(i_K),n_iter+1);
    My_cov=zeros(n_dim,n_dim,K(i_K),n_iter+1);
    My_phi=zeros(n_obs,K(i_K),n_iter+1);
    My_n=zeros(K(i_K),n_iter+1);
    My_obj_ft=zeros(n_iter+1,1);
    
    %--------------------------------------
    % initialize My_parameters in some way
    %--------------------------------------
    
    rng(n_rng);
    My_pi(:,1)=1/K(i_K)*ones(1,K(i_K));
    rand_mean=round(abs(max(max(X)))/2);
    My_mu(:,:,1)=X_mean*ones(1,K(i_K))+rand_mean*(rand(n_dim,K(i_K))-0.5);
    for ii_K=1:K(i_K)
        My_cov(:,:,ii_K,1)=X_cov;
    end;
    My_obj_ft_temp=zeros(n_obs,1);
    for ii_K=1:K(i_K)
        My_obj_ft_temp=My_obj_ft_temp + My_pi(ii_K,1)*mvnpdf(X',My_mu(:,ii_K,1)',My_cov(:,:,ii_K,1));
    end;
    My_obj_ft(1)=sum(log(My_obj_ft_temp));
    
    for i_iter=2:n_iter+1
        
        %--------------------------------------
        % E-step
        %--------------------------------------
        
        for ii_K=1:K(i_K)
            My_phi(:,ii_K,i_iter)=My_pi(ii_K,i_iter-1)*mvnpdf(X',My_mu(:,ii_K,i_iter-1)',My_cov(:,:,ii_K,i_iter-1));
        end;
        My_phi(:,:,i_iter)=My_phi(:,:,i_iter)./(sum(My_phi(:,:,i_iter),2)*ones(1,K(i_K)));
        
        %--------------------------------------
        % M-step
        %--------------------------------------
        
        My_n(:,i_iter)=squeeze(sum(My_phi(:,:,i_iter)));
        My_pi(:,i_iter) =My_n(:,i_iter)/n_obs;
        for ii_K=1:K(i_K)
            My_mu(:,ii_K,i_iter)=sum(ones(n_dim,1)*My_phi(:,ii_K,i_iter)'.*X,2)/My_n(ii_K,i_iter);
            My_cov(:,:,ii_K,i_iter)=((X-My_mu(:,ii_K,i_iter)*ones(1,n_obs)).*(ones(2,1)...
                *sqrt(My_phi(:,ii_K,i_iter))'))*((X'-ones(n_obs,1)*My_mu(:,ii_K,i_iter)')...
                .*(sqrt(My_phi(:,ii_K,i_iter))*ones(1,2)))/My_n(ii_K,i_iter);
        end;
        My_obj_ft_temp=zeros(n_obs,1);
        for ii_K=1:K(i_K)
            My_obj_ft_temp=My_obj_ft_temp + My_pi(ii_K,i_iter)*mvnpdf(X',My_mu(:,ii_K,i_iter)',My_cov(:,:,ii_K,i_iter));
        end;
        My_obj_ft(i_iter)=sum(log(My_obj_ft_temp));
        
    end;
    
    %--------------------------------------
    % Plot
    %--------------------------------------
    
    subplot(2,4,i_K)
    plot(1:n_iter,My_obj_ft(2:end))
    xlim([1 n_iter]);
    ylim([-1500 -1000]);
    xlabel('iteration')
    ylabel('log likelihood')
    title(['K=',num2str(K(i_K))]);
    set(gca,'FontSize',12,'FontWeight','b');
    
    subplot(2,4,i_K+4)
    [~,My_phi_M] = max(My_phi(:,:,n_iter),[],2);
    lgd=cell(K(i_K),1);
    for ii_K = 1:K(i_K)
        hold on;
        scatter(X(1,find(My_phi_M==ii_K)),X(2,find(My_phi_M==ii_K)),'*')
        lgd{ii_K}=['Cluster ' num2str(ii_K)];
    end;
    xlabel('Dimension 1')
    ylabel('Dimension 2')
    title(['K=',num2str(K(i_K))]);
    legend(lgd)
    set(gca,'FontSize',12,'FontWeight','b');
    axis square
    
end