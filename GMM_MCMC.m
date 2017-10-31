%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gibbs Sampling for Gaussian Mixture Model
% Jie Yang 2015
%
% Note: All parameters named by My_** are
%       parameters to be iterated/optimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GMM_MCMC()

load('data.mat')
n_obs=size(X,2);
n_dim=size(X,1);

n_iter=500;

X_mean=mean(X,2);
X_cov=cov(X');

c0=0.1;
a=n_dim;
A=X_cov;
B=c0*n_dim*A;

alpha=1;

My_phi=zeros(n_obs,n_obs+1,n_iter+1);
My_lambda=zeros(n_dim,n_dim,n_obs+1,n_iter+1);
My_mu=zeros(n_dim,n_obs+1,n_iter+1);
My_inds=zeros(n_obs,n_iter+1);
My_nK=zeros(n_iter+1,1);
ind_cnt=zeros(6,n_iter+1);

%--------------------------------------
% initialize My_parameters in some way
%--------------------------------------

n_rng=0;
rng(n_rng);
My_lambda(:,:,1,1)=wishrnd(inv(B),a);
My_mu(:,1,1)=mvnrnd(X_mean',inv(c0*My_lambda(:,:,1,1)));
My_inds(:,1)=1; % index of cluster per observation
My_nK(1)=1; % number of cluster

for i_iter=2:n_iter+1
    
    My_mu(:,:,i_iter)=My_mu(:,:,i_iter-1);
    My_lambda(:,:,:,i_iter)=My_lambda(:,:,:,i_iter-1);
    My_inds(:,i_iter)=My_inds(:,i_iter-1);
    
    for i_obs=1:n_obs
        
        My_nK_temp=max(My_inds(:,i_iter));
        for i_K=1:My_nK_temp
            if sum(My_inds(setdiff(1:n_obs,i_obs),i_iter)==i_K)
                My_phi(i_obs,i_K,i_iter)=mvnpdf(X(:,i_obs)',...
                    My_mu(:,i_K,i_iter)',...
                    inv(My_lambda(:,:,i_K,i_iter)))*sum(My_inds(setdiff(1:n_obs,i_obs),i_iter)==i_K)/(alpha+n_obs-1);
            end;
        end;
        i_K=My_nK_temp+1;
        
        %--------------------------------------
        % iteration of My_parameters
        %--------------------------------------
        
        My_phi(i_obs,i_K,i_iter)=alpha/(alpha+n_obs-1)*(c0/(pi*(1+c0)))^(n_dim/2)/(det(B))^(-a/2)*...
            det(B+c0/(1+c0)*(X(:,i_obs)-X_mean)*(X(:,i_obs)-X_mean)')^(-(a+1)/2)*...
            exp(sum(gammaln((a+1)/2+(1-(1:n_dim))/2)-gammaln(a/2+(1-(1:n_dim))/2)));
        My_phi(i_obs,:,i_iter)=My_phi(i_obs,:,i_iter)/sum(My_phi(i_obs,:,i_iter),2);
        My_inds(i_obs,i_iter)=sum(rand(1)>cumsum(My_phi(i_obs,:,i_iter)))+1;
        
        if (My_inds(i_obs,i_iter)==i_K)
            s=1;c1=c0+s;a1=a+s;
            m1=c0/(c0+s)*X_mean+1/(s+c0)*X(:,i_obs);
            B1=B+s/(a*s+1)*(X(:,i_obs)-X_mean)*(X(:,i_obs)-X_mean)';
            My_lambda(:,:,i_K,i_iter)=wishrnd(inv(B1),a1);
            My_mu(:,i_K,i_iter)=mvnrnd(m1,inv(c1*My_lambda(:,:,i_K,i_iter)));
        end;
    end;
    
    My_nK(i_iter)=max(My_inds(:,i_iter));
    ind_delete=[];
    for i_K=1:My_nK(i_iter)
        s=sum(My_inds(:,i_iter)==i_K);
        if ~s
            ind_delete=[ind_delete i_K];
        else
            c1=c0+s;a1=a+s;
            m1=c0/(c0+s)*X_mean+1/(s+c0)*sum(X(:,My_inds(:,i_iter)==i_K),2);
            Xmean=mean(X(:,My_inds(:,i_iter)==i_K),2);
            B1=B+s/(a*s+1)*(Xmean-X_mean)*(Xmean-X_mean)'+(X(:,My_inds(:,i_iter)==i_K)...
                -Xmean*ones(1,s))*(X(:,My_inds(:,i_iter)==i_K)-Xmean*ones(1,s))';
            My_lambda(:,:,i_K,i_iter)=wishrnd(inv(B1),a1);
            My_mu(:,i_K,i_iter)=mvnrnd(m1,inv(c1*My_lambda(:,:,i_K,i_iter)));
        end;
    end;
    
    % new cluster inds
    if ~isempty(ind_delete)
        ind_save=setdiff(1:My_nK(i_iter),ind_delete);
        mu_temp=My_mu(:,ind_save,i_iter);
        lambda_temp=My_lambda(:,:,ind_save,i_iter);
        My_nK(i_iter)=My_nK(i_iter) - length(ind_delete);
        My_mu(:,1:My_nK(i_iter),i_iter)=mu_temp;
        My_lambda(:,:,1:My_nK(i_iter),i_iter)=lambda_temp;
        
        My_inds_temp=My_inds(:,i_iter);
        for k=1:length(ind_save)
            My_inds_temp(My_inds_temp == ind_save(k))=k;
        end;
        My_inds(:,i_iter)=My_inds_temp;
    end;
    
    ind_cnt_temp=histc(My_inds(:,i_iter),unique(My_inds(:,i_iter)));
    if length(ind_cnt_temp)>6
        ind_cnt_temp=ind_cnt_temp(1:6);
    end;
    ind_cnt(1:length(ind_cnt_temp),i_iter)=ind_cnt_temp;
    
    disp(['iteration time = ' num2str(i_iter)])
end;
ind_cnt=sort(ind_cnt)';
ind_cnt=flip(ind_cnt,2);

%--------------------------------------
% Plot
%--------------------------------------

figure;
subplot(1,3,1);
plot(1:n_iter,ind_cnt(2:end,:))
legend('cluster 1','cluster 2','cluster 3','cluster 4','cluster 5','cluster 6')
xlabel('iteration');
ylabel('number of observations per cluster')
set(gca,'FontSize',12,'FontWeight','b');

subplot(1,3,2);
plot(1:n_iter,My_nK(2:end))
xlim([1 n_iter])
ylim([1 max(My_nK)])
xlabel('iteration');
ylabel('number of clusters')
set(gca,'FontSize',12,'FontWeight','b');

subplot(1,3,3);
lgd=cell(My_nK(n_iter+1),1);
for ii_K = 1:My_nK(n_iter+1)
    hold on;
    dnX=find(My_inds(:,n_iter+1)==ii_K);
    scatter(X(1,dnX),X(2,dnX),'*')
    lgd{ii_K}=['Cluster ' num2str(ii_K)];
end;
xlabel('Dimension 1')
ylabel('Dimension 2')
legend(lgd);
set(gca,'FontSize',12,'FontWeight','b');
axis square


