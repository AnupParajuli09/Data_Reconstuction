%%%% ============= Node Wise Dissimilar Load Profile: PCA Spatial Data Reconstruction ============= %%%%

clear all
close all
clc

file='Data.xlsx';

% =================  Bus power data of IEEE 123 bus ================= %
data1=xlsread(file,1); % Bus power data
bus=data1(:,1);
phA_kw=data1(:,2);
phA_kVAR=data1(:,3);
phB_kw=data1(:,4);
phB_kVAR=data1(:,5);
phC_kw=data1(:,6);
phC_kVAR=data1(:,7);
T_kw=data1(:,8);
T_kVAR=data1(:,9);
Nbus=max(bus);

% ============ Sample load profile generation: 1-min ============== %
% T=96; 
T=24*60;

% % ==== Dissimilar Load Profile ====
data2=xlsread(file,5);
trainpu=data2(2:1441,2:126);
% trainpu=data2(2:97,2:126);
% plot(data2(2:1441,1),trainpu)
% plot(data2(2:97,1),trainpu)
% L=zeros(T,Nbus);
% for n=1:Nbus
%     for t=1:T
%         L(t,n)=(1.3-0.8)*rand+0.8;
%     end
% end
% plot(L)

% % Test load profile 15-min
data3=xlsread(file,6);
testpu=data2(2:97,2:126);
% figure
% plot(data3(2:97,1),testpu)
T1=24*60/15;
% for n=1:Nbus
% for t=1:T1
%     Lt(t,n)=(1.3-0.8)*rand+0.8;
% end
% end


% % ======== Offline calculation of PCA parameters: Bus wise ========== %
P_L_demo=zeros(Nbus,T);
for n=1:Nbus
    for t=1:T
        P_L_demo(n,t)=trainpu(t,n)*T_kw(n,1);
    end
end
size_P=size(P_L_demo);

% Mean vector calculation
P_avg=zeros(Nbus,1);
for n=1:Nbus
    S=0;
    for t=1:T
        S=S+P_L_demo(n,t);
    end
    P_avg(n,1)=(1/T)*S;
end

% Correlation Matrix Calculation
Corr_P_L=zeros(Nbus,Nbus);
for m=1:Nbus
    for n=1:Nbus
        S1=0; S2=0; S3=0;
        for t=1:T
            S1=S1+(P_L_demo(m,t)-P_avg(m,1))*(P_L_demo(n,t)-P_avg(n,1));
            S2=S2+power((P_L_demo(m,t)-P_avg(m,1)),2);
            S3=S3+power((P_L_demo(n,t)-P_avg(n,1)),2);
        end
        Corr_P_L(m,n)=S1/sqrt(S2*S3);
    end
end



% % Covariance matrix calculation
% P_L_norm=P_L_demo-P_avg;
% % Cov_P_L=cov(P_L_norm)
% Cov_P_L=zeros(Nbus,Nbus);
% for m=1:Nbus
%     for n=1:Nbus
%         S=0;
%         for t=1:T
%             S=S+(P_L_demo(m,t)-P_avg(m,1))*(P_L_demo(n,t)-P_avg(n,1));
%         end
%         Cov_P_L(m,n)=(1/T)*S;
%     end
% end
% size_cov=size(Cov_P_L);
% Cov_P_L_pu=Cov_P_L/1000;
% 
% % Eigen values (D) and vectors (V) of covariance matrix
% [V,D]=eig(Cov_P_L);
% 
% % Sort eigen vectors
% non_sortV=V;
% colSums = sum(V);
% [sumsorted , sortedIndices] = sort(colSums, 'descend');  % Sort indices in descending order
% sortedV = V(:, sortedIndices);  % Rearrange columns based on sorted indices
% 
% %%%%%%%%%%%%%%% Online data collection and data recovery %%%%%%%%%%%%%
% itrmax=50;
% for itr=1:itrmax
% Ttest=5; % Testing time interval
% L=50; % Total number of data collecting nodes
% 
% 
% for n=1:Nbus
%     for t=1:T1
%         S_L_test(n,t)=testpu(t,n)*T_kw(n,1);
%     end
% end
% % figure
% % plot(S_L_test(:,1))
% 
% 
% % UL=sortedV(:,1:L);
% UL=sortedV(:,2);
% S_L=zeros(Nbus,1); % Original power at collected nodes
% % n_L=(randperm(Nbus, L))'; % Index of selected node 
% 
% 
% % Identify nodes having loads
% k=1;
% for n=1:Nbus
%     if T_kw(n,1)~=0
%         n_load(k,1)=n;
%         k=k+1;
%     end
% end
% 
% i_L=(randperm(length(n_load), L))'; % Index of selected item from n_load
% 
% % Selected Node and their actual measurement
% n_L=zeros(L,1);
% for l=1:L
%     n_L(l,1)=n_load(i_L(l,1),1);
% end
% 
% for l=1:L
%     S_L(n_L(l,1),1,1)=testpu(Ttest,n_L(l,1))*T_kw(n_L(l,1),1)-P_avg(n_L(l,1),1);
% end
% 
% % Recovered data
% S_L_rec=P_avg+UL*UL'*S_L;
% 
% % Integrated normalized absolute error (INAE) and Data Plot
% S_test_tot=zeros(1,T1);
% for t=1:T1
%     for n=1:Nbus
%         S_test_tot(1,t)=S_test_tot(1,t)+S_L_test(n,t);
%     end
% end
% 
% 
% S_rec_tot=0;
% for n=1:Nbus
%     S_rec_tot=S_rec_tot+abs(S_L_test(n,Ttest)-S_L_rec(n,1));
% end
% INAE=S_rec_tot/S_test_tot(1,Ttest)
% 
% Error(itr,1)=INAE*100;
% end
% 
% % figure
% % T=1:Nbus;
% % p=plot(T,[S_L_test(:,Ttest),S_L_rec(:,1)])
% % legend(["Original" "Reconstruct"])
% % p(1).Marker = 'o';
% % p(2).Marker = '*';
% % title('Reconstructed power data')
% % ylabel('Power (kW)')
% % xlabel('Bus number')
% 
% 
