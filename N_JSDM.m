% EMAIL:songyc@njupt.edu.cn
clear all
close all
clc
derad = pi/180;
radeg = 180/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('co64_10.mat');AS=10;
K=32;                    % number of users
kelma =64;               % number of antennas
NEA=10;                  %NAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lm=300;
Tc=100;                   %coherence time
twpi = 2*pi;
GT=kelma/K;
d1=0:1/2:(kelma-1)/2;     %
%%%%%%%%%%%%%%azimuth group

%%%%%%%%%%%%%%%%%%%%%%%%
for SNR=1:8
    snr=5*(SNR-3);
    rho=1/10^(snr/10);
    C1(SNR)=0;
    C2(SNR)=0;
    C3(SNR)=0;
    for pp=1:Lm
        pp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%system model
        Num=0;
        for ii=1:K
            THA(ii)=rand*240;
            THa(ii)=(floor(THA(ii)))/2-60;
        end
        THa=sort(THa);THA=sort(THA);
        for ii=1:K
            RA{ii}=Cod{floor(THA(ii))+1};
        end
        
        %%%%%%%%%%%%%%%%%%% prebeamformer
        WWW=[];eeee=zeros(K,1);eeee1=zeros(K,1);
        for ii=1:K
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%optimal N-JSDM
            WT2=zeros(kelma,0);
            wa{ii}=zeros(kelma,0);
            OM{ii}=ones(K,1);
            RETempa=eye(kelma)*0;
            RET{ii}=eye(kelma)*0;
            for kk4=1:K
                if abs(THa(kk4)-THa(ii))>NEA
                    OM{ii}(kk4)=0;
                    OM2{ii}(kk4,:)=0;
                    RETempa=RETempa+RA{kk4};
                else
                    if kk4<ii
                        WT2=[WT2 wa{kk4}];
                    end
                    RET{ii}=RET{ii}+RA{kk4};
                end
            end
            CA=[];        WT3=null(WT2');
            [A2 B2]=eig(WT3'*RETempa*WT3);BB2=diag(B2);
            for iii=1:size(BB2,1)
                if real(BB2(iii))<0.01;
                    CA=[CA,A2(:,iii)];
                end
            end
            if size(CA,2)==0
            else
                EEE=CA;
                CA1=WT3*CA;
                EEE=CA1;
                TTA=EEE'*RA{ii}*EEE;
                [VT1 DD]=eig(TTA);
                ad=diag(abs(DD));
                V1=zeros(size(VT1,1),0);
                mmax=max(ad);
                while max(ad)>=mmax/10
                    [b aw]=max(ad);
                    ad(aw)=-Inf;
                    V1=[V1,VT1(:,aw)];
                end
                wa{ii}=[wa{ii} EEE*V1];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%constrained N-JSDM
            WT2=zeros(kelma,0);tt=size(wa{ii},2);
            wa1{ii}=zeros(kelma,0);
            RETempa=eye(kelma)*0;
            nnum=0;
            for kk4=1:K
                if abs(THa(kk4)-THa(ii))>NEA
                    RETempa=RETempa+RA{kk4};
                else
                    if kk4<ii
                        nnum=nnum+1;
                        WT2=[WT2 wa1{kk4}];
                    end
                end
            end
            CA=[];        WT3=null(WT2');
            [A2 B2]=eig(WT3'*RETempa*WT3);BB2=diag(B2);
            for iii=1:size(BB2,1)
                if real(BB2(iii))<0.01
                    CA=[CA,A2(:,iii)];
                end
            end
            G=floor((nnum+1)*GT)-size(WT2,2);
            
            CA2=[];
            
            [AE2 BE2]=eig(RA{ii});BB2=diag(BE2);
            tttt=max(real(BB2))/10;
            for iii=1:size(BB2,1)
                if real(BB2(iii))>=tttt
                    CA2=[CA2,AE2(:,iii)];
                end
            end
            
            RATE=CA2*CA2';
            if size(CA,2)==0
            elseif size(CA,2)==1
                G1=1;
                EEE=CA;
                CA1=WT3*CA;
                EEE=CA1;
                TTA=EEE'*RATE*EEE;
                [VT1 DD]=eig(TTA);numt=min(tt,round(G1));
                aw=zeros(numt,1);
                ad=diag(abs(DD));
                for iii=1:numt;
                    [b aw(iii)]=max(ad);
                    ad(aw(iii))=-Inf;
                end
                V1=VT1(:,aw);
                wa1{ii}=EEE*V1;
                wa1{ii}=[wa1{ii},zeros(kelma,G-G1)];
            else 
                EEE=CA;numt=min(tt,round(G));
                CA1=WT3*CA;
                EEE=CA1;
                TTA=EEE'*RATE*EEE;
                [VT1 DD]=eig(TTA);
                aw=zeros(numt,1);
                ad=diag(abs(DD));
                for iii=1:numt;
                    [b aw(iii)]=max(ad);
                    ad(aw(iii))=-Inf;
                end
                V1=VT1(:,aw);
                wa1{ii}=EEE*V1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% second stage precoding
        
        HT=[];
        for ii=1:K
            [AA BB CC]=svd(RA{ii});
            RA1{ii}=AA*sqrt(BB)*CC';
            HTe{ii}=0.5*(randn(1,kelma)+randn(1,kelma)*j)*RA1{ii};
            HT=[HT;HTe{ii}];
        end
        for ii=1:K
            W{ii}=wa{ii};
            PT{ii}=(HT*W{ii}).*OM{ii};
            W11{ii}=wa1{ii};
            PT1{ii}=(HT*W11{ii}).*OM{ii};
        end
        WTT=[];HTT=[];WTT1=[];HTT1=[];
        for ii=1:K
            WTT=[WTT wa{ii}];
            HTT=[HTT PT{ii}];
            WTT1=[WTT1 wa1{ii}];
            HTT1=[HTT1 PT1{ii}];
        end
        % HTT
        for kk=1:size(HTT,1)
            ttem=HTT(kk,:);
            ttem1=HTT1(kk,:);
            n(kk)=numel(ttem(ttem~=0));
            n1(kk)=numel(ttem1(ttem1~=0));
        end
        num=max(n);                           %DTL  optimal N-JSDM
        num1=max(n1);                         %DTL  constrained N-JSDM
        P=pinv(HTT);
        P1=pinv(HTT1);vall=size(HTT1,2);
        
        for ii=1:size(P,2)
            P(:,ii)=P(:,ii)/norm(P(:,ii));
        end
        for ii=1:size(P1,2)
            P1(:,ii)=P1(:,ii)/norm(P1(:,ii));
        end
        WT=WTT*P;WT1=WTT1*P1;
        H=HT;
        W1=WT;
        W2=pinv(HT);
        for ii=1:size(W2,2)
            W2(:,ii)=W2(:,ii)/norm(W2(:,ii));
        end
        W3=WT1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%effective spectral efficiency
        Cs=0;
        Cs2=0;
        Cs3=0;
        for kk=1:K
            Su=0;
            Su2=0;
            Su3=0;
            for ii=1:K
                Su=abs(H(kk,:)*W1(:,ii))^2+Su;
                Su2=abs(H(kk,:)*W2(:,ii))^2+Su2;
                Su3=abs(H(kk,:)*W3(:,ii))^2+Su3;
            end
            Ct(kk)=log2(1+abs((H(kk,:)*W1(:,kk)))^2/(1/10^(snr/10)+Su-abs(H(kk,:)*W1(:,kk))^2)); %optimal N-JSDM
            Ct2(kk)=log2(1+abs((H(kk,:)*W2(:,kk)))^2/(1/10^(snr/10)+Su2-abs(H(kk,:)*W2(:,kk))^2));%full precoding
            Ct3(kk)=log2(1+abs((H(kk,:)*W3(:,kk)))^2/(1/10^(snr/10)+Su3-abs(H(kk,:)*W3(:,kk))^2));%constrained N-JSDM
            Cs=Cs+Ct(kk)*(1-num/Tc);
            Cs2=Cs2+Ct2(kk)*(1-kelma/Tc);
            Cs3=Cs3+Ct3(kk)*(1-num1/Tc);
        end
        C1(SNR)=C1(SNR)+Cs;
        C2(SNR)=C2(SNR)+Cs2;
        C3(SNR)=C3(SNR)+Cs3;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
m1=C1/Lm; %optimal N-JSDM
m2=C2/Lm;%full precoding
m3=C3/Lm;%constrained N-JSDM
x=1:8;
plot(x,m1)
hold on
plot(x,m2,'r-*')
plot(x,m3,'k-*')