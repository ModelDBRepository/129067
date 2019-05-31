%This script reproduces the grid and place cell interference patterns of
%Fig 7A from Welday et al.

figure(1);
v=25; % speed in cm/sec 
Freq=7; % theta base freq in Hz
    W=2*pi*Freq; % angular freq
d = [23.1579    28.5766    37.2030    53.3025    71.0747    79.9804   91.4141   106.6508]*3.8e-3;    %slopes of speed modulation for theta cell VCOs                  
Npt=1001; % number of data points on the time axis
t=linspace(-2,2,Npt);  % time axis: -2 sec to 2 sec

clf
% plot oscillators
for i=1:8
subplot(12,1,i)
    plot(t,cos((W+d(i)*v)*t))
    set(gca,'xtick',[]); axis tight;
end

% compute grid cells 
s=cos((W+d(2)*v)*t)+cos((W+d(8)*v)*t); % use cell 1 and cell 7
a=cos(d(2)*v*t)+cos(d(8)*v*t);
b=sin(d(2)*v*t)+sin(d(8)*v*t);
en=sqrt(a.^2+b.^2); % envelope
ph=unwrap(angle(a+sqrt(-1)*b));  % phase
% unwrap phase
for i=2:length(ph)
    if ph(i)<ph(i-1)-pi/2
        ph(i:end)=ph(i:end)+pi;
    elseif ph(i)>ph(i-1)+pi/2
        ph(i:end)=ph(i:end)-pi;
    end
end
ca=cos(W*t+ph); % carrier

% plot grid cells
subplot(12,1,9)
    plot(t, s)
    hold on
    plot(t,en,'r-')
    set(gca,'xtick',[]); axis tight;
subplot(12,1,10)
    plot(t,ca, 'g-')
    set(gca,'xtick',[]); axis tight;
    
% compute place cells 
S=0;A=0;B=0;
slp=[];
wv=[];
for i=1:8 % add all oscillators
    S=S+cos((W+d(i)*v)*t);
    wv=[wv; cos((W+d(i)*v)*t)];
    slp(i)=mean(diff(((W+d(i)*v)*t)'));
    A=A+cos(d(i)*v*t);
    B=B+sin(d(i)*v*t);
end
EN=sqrt(A.^2+B.^2); %envelope
% PH=atan(B./A); % phase
PH=unwrap(angle(A+sqrt(-1)*B));  % phase
% unwrap phase
for i=2:length(PH)
    if PH(i)<PH(i-1)-pi/2
        PH(i:end)=PH(i:end)+pi;        
    elseif PH(i)>PH(i-1)+pi/2
        PH(i:end)=PH(i:end)-pi;
    end     
end
CA=cos(W*t+PH); % carrier
%CA=cos(W*t+PH+pi); % carrier (phase + pi)
    
% plot place cells
subplot(12,1,11)
    plot(t, S)
    hold on
    plot(t,EN,'r-')
    set(gca,'xtick',[]); axis tight;
subplot(12,1,12)
    plot(t,CA, 'g-')
    set(gca,'xtick',[]); axis tight;

    figure(1);