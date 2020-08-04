% EMAIL:songyc@njupt.edu.cn   ¿¼ÂÇÔëÉù
clear all
close all
clc
derad = pi/180;
radeg = 180/pi;
twpi = 2*pi;
kelma =32;               % azimuth   
d1=0:1/2:(kelma-1)/2;     % 
%%%%%%%%%%%%%%azimuth group
ASa=8;
%%%%%%%%%%%%%%%%%%%%%%%%
N=240;
THa=zeros(N,1);
lowg=-60;upg=60;
for ii=1:N+1
    ii
    THa(ii)=lowg+(upg-lowg)/N*(ii-1);
    Cod{ii}=cova(THa(ii),ASa,kelma);
end