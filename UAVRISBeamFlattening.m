clear all
close all
clc
global P sigma2 beta0 M N H w dx dy lambda
H=100; %RIS altitude
sigma2dBm=-110;
sigma2=10^(sigma2dBm/10)/1000;
PdBm=20;
P=10^(PdBm/10)/1000;
beta0dB=-40;
beta0=10^(beta0dB/10);
fc=2.4e9;
My=8;
Mz=8;
M=My*Mz; %number of BS antennas
c=3e8;
lambda=c/fc;
dx=lambda/10;
dy=lambda/10;
Nx=400;
Ny=1;
N=Nx*Ny; %number of RIS elements
q=[100,0].'; %RIS horizontal location 
w=[1000,0].'; %user location
phiR=asin(norm(q)/norm([q;H])); %OK
etaR=atan(q(2)/q(1)); %OK
PhiR=sin(phiR)*cos(etaR);
OmegaR=sin(phiR)*sin(etaR);
%receive array response of RIS:
aR=kron([exp(-1i*2*pi*(0:Nx-1)*dx/lambda*sin(phiR)*cos(etaR))].',[exp(-1i*2*pi*(0:Ny-1)*dy/lambda*sin(phiR)*sin(etaR))].'); %(5)
phiR2=pi/2-asin(norm(q)/norm([q;H]));
etaR2=atan(q(2)/q(1)); 
aTs=kron([exp(-1i*2*pi*(0:My-1)*dx/lambda*sin(phiR2)*cos(etaR2))].',[exp(-1i*2*pi*(0:Mz-1)*dy/lambda*sin(phiR2)*sin(etaR2))].'); %(???)
betaG=beta0/(H^2+norm(q)^2); %(1)
dG=norm([q;H]);
G=sqrt(betaG)*exp(-1i*2*pi*dG/lambda)*aR*aTs'; %(6)
phiT=asin(norm(q-w)/norm([q;H]-[w;0])); %OK
etaT=acos((w(1)-q(1))/norm(q-w)); %OK
PhiT=sin(phiT)*cos(etaT);
OmegaT=sin(phiT)*sin(etaT);
aT=kron([exp(-1i*2*pi*(0:Nx-1)*dx/lambda*sin(phiT)*cos(etaT))].',[exp(-1i*2*pi*(0:Ny-1)*dy/lambda*sin(phiT)*sin(etaT))].'); %(7)
betah=beta0/(H^2+norm(q-w)^2); %(2)
dh=sqrt(H^2+norm(w-q)^2);
hH=sqrt(betah)*exp(-1i*2*pi*dh/lambda)*aT'; %(8)
%theta=2*pi*rand(1,N);%random phase shifts 
theta=zeros(1,N);%random phase shifts 
for ny=1:Ny
    for nx=1:Nx
        %fprintf('nx %d, ny %d, (nx-1)*Ny+nx %d \n',nx,ny, (nx-1)*Ny+nx)
        theta((nx-1)*Ny+ny)=0-2*pi*(nx-1)*dx/lambda*(PhiT-PhiR)-2*pi*(ny-1)*dy/lambda*(OmegaT-OmegaR); %(13)
    end
end
Theta=diag(exp(1i*theta));
%v is the BS precoding vector (it is the eigenvector corresponding to the largest eigenvalue of the matrix htildeH*htildeH')
v=aTs/sqrt(M); %proposition 1
gamma=P/sigma2*abs(hH*Theta*G*v)^2; %(10) 
gamma1=P/sigma2*beta0^2*M*N^2/(H^2+norm(q-w)^2)/(H^2+norm(q)^2); %(14) 
gamma1dB=pow2db(gamma1);
rho=norm(w)/H; %distance over height ratio
ksi=1/2-sqrt(1/4-1/rho^2);
qstar=ksi*norm(w); % (15) optimal location of RIS

fprintf('gamma is %d, gamma1 is %d, rho = %d, qstar is %d  meters\n', gamma, gamma1, rho, qstar)
drawfig8()

Ns=20
Delta=1/(Ns*dx/lambda)



% PhiR=q(1)/norm([q;H]);
% OmegaR=q(2)/norm([q;H]);
% PhiT=(w(1)-q(1))/norm([q;H]-[w;0]);
% OmegaT=(w(2)-q(2))/norm([q;H]-[w;0]);


%Deltaspan=abs((CalculatePhiT(q,[xl;0])-CalculatePhiR(q,[xl;0]))-(CalculatePhiT(q,[xu;0])-CalculatePhiR(q,[xu;0])))
xl=155;
xu=325;
drawfig9and10(256, xl,xu)

drawfigure11(256, 500, 1500, -300,+300)

drawfigure12(20, 500, 1500, -300,+300)

function PhiR=CalculatePhiR(q,w)
    global H
    PhiR=q(1)/norm([q;H]);
end

function OmegaR=CalculateOmegaR(q,w)
    global H
    OmegaR=q(2)/norm([q;H]);
end

function PhiT=CalculatePhiT(q,w)
    global H    
    PhiT=(w(1)-q(1))/norm([q;H]-[w;0]);
end

function OmegaT=CalculateOmegaT(q,w)
    global H
    OmegaT=(w(2)-q(2))/norm([q;H]-[w;0]);
end

function drawfig8()
    global P sigma2 beta0 M H w
    figure;
    hold on;
    Narray=200:100:800;
    q=[10.1,0];
    gammaarray1=pow2db(P/sigma2*beta0^2*M*Narray.^2/(H^2+norm(q-w)^2)/(H^2+norm(q)^2));
    q=[500,0];
    gammaarray2=pow2db(P/sigma2*beta0^2*M*Narray.^2/(H^2+norm(q-w)^2)/(H^2+norm(q)^2))  ;
    plot(Narray, gammaarray1)
    plot(Narray, gammaarray2)
    grid on
    hold off
end

function drawfig9and10(N,xl,xu)
    global P sigma2 beta0 M N H dx dy lambda
    w=[1000;0];
    qxarray=-150:10:(xl+xu)/2;
    Lstararray=[];
    f2array=[];
    f1array=[];
    for i=1:numel(qxarray)
        q=[qxarray(i);0];
        Deltaspan=abs((CalculatePhiT(q,[xl;0])-CalculatePhiR(q,[xl;0]))-(CalculatePhiT(q,[xu;0])-CalculatePhiR(q,[xu;0])));
        Lstararray=[Lstararray,ceil(sqrt(Deltaspan*N*dx/lambda))];
        f2array=[f2array,max((H^2+norm(q-[xl;0])^2)*(H^2+norm(q)^2),(H^2+norm(q-[xu;0])^2)*(H^2+norm(q)^2))];
        f1array=[f1array,4/pi^2*M*N^2/ceil(sqrt(Deltaspan*N*dx/lambda))^2]; %(39)
    end
    figure
    hold  on
    title('Figure 9')
    yyaxis left
    ylabel('The required number of sub arrays')
    stairs(qxarray,Lstararray)
    yyaxis right
    ylabel('The worst case concatenated pat loss (dB)')
    plot(qxarray,80+pow2db(f2array))
    grid on
    xlabel('AIRS placement along the x axis q_x (m)')
    hold off
    figure
    plot(qxarray,pow2db(P*beta0^2*f1array./f2array/sigma2))
    title('Figure 10')
    ylabel('The worst case SNR (dB)')
    xlabel('AIRS placement along the x axis q_x (m)')
    grid on
end


function drawfigure11(N,xl,xu, yl, yu)
    global sigma2 beta0 M H dx dy lambda 
    PdBmarray=20:30;
    q=[1000,0].'; %RIS horizontal location 
    Nx=N; Ny=1;
    f1array=[];
    f2array=[];
    Deltaspanx=abs((CalculatePhiT(q,[xl;yl])-CalculatePhiR(q,[xl;yl]))-(CalculatePhiT(q,[xu;yu])-CalculatePhiR(q,[xu;yu])));
    Deltaspany=abs(CalculateOmegaT(q,[xl;yl])-CalculateOmegaT(q,[xu;yu]));
    for i=1:numel(PdBmarray)
        P=10^(PdBmarray(i)/10)/1000;
        f1array=[f1array,16/pi^4*Nx^2/ceil(sqrt(Deltaspanx*Nx*dx/lambda))^2*Ny^2/ceil(sqrt(Deltaspany*Ny*dy/lambda))^2];
        f2array=[f2array,max((H^2+norm(q-[xl;yl])^2)*(H^2+norm(q)^2),(H^2+norm(q-[xu;yu])^2)*(H^2+norm(q)^2))];
    end
    figure
    hold
    plot(PdBmarray,pow2db(10.^(PdBmarray./10)./1000.*beta0^2.*M.*f1array./f2array/sigma2), marker='*') %3dB shifted for some reason
    title('Figure 11')
    ylabel('The worst case SNR (dB)')
    xlabel('Number of AIRS reflecting elements') 
    grid on
end

function drawfigure12(Ny,xl,xu, yl, yu)
    global P sigma2 beta0 M H dx dy lambda 
    Narray=400:100:1000;
    q=[100,0].'; %RIS horizontal location 
    %Ny=20;
    f1array=[];
    f2array=[];
    f3array=[];
    for i=1:numel(Narray)
        Nx=Narray(i)/Ny;
        
        Deltaspanx=abs((CalculatePhiT(q,[xl;yl])-CalculatePhiR(q,[xl;yl]))-(CalculatePhiT(q,[xu;yu])-CalculatePhiR(q,[xu;yu])));
        Deltaspany=abs(CalculateOmegaT(q,[xl;yl])-CalculateOmegaT(q,[xu;yu]));
        16/pi^4*Nx^2/ceil(sqrt(Deltaspanx*Nx*dx/lambda))^2*Ny^2/ceil(sqrt(Deltaspany*Ny*dy/lambda))^2;
        f1array=[f1array,16/pi^4*Nx^2/ceil(sqrt(Deltaspanx*Nx*dx/lambda))^2*Ny^2/ceil(sqrt(Deltaspany*Ny*dy/lambda))^2];
        f2array=[f2array,max((H^2+norm(q-[xl;yl])^2)*(H^2+norm(q)^2),(H^2+norm(q-[xu;yu])^2)*(H^2+norm(q)^2))];
        f3array=[f3array,4/pi^2*Nx^2*abs(sin(pi*Ny*dy/lambda*(CalculateOmegaT(q,[xu;yu])-CalculateOmegaR(q,[xu,yu])))/sin(pi*dy/lambda*(CalculateOmegaT(q,[xu;yu])-CalculateOmegaR(q,[xu,yu]))))^2/ceil(sqrt(Deltaspanx*Nx*dx/lambda))^2];
    end
    figure
    hold
    plot(Narray,pow2db(P*beta0^2*M*f1array./f2array/sigma2), marker='*')
    plot(Narray,pow2db(P*beta0^2*M*f3array./f2array/sigma2), marker='o')%f3array is wrong
    title('Figure 12')
    ylabel('The worst case SNR (dB)')
    xlabel('Number of AIRS reflecting elements')
    grid on
end