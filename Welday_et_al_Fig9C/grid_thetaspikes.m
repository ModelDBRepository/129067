%%%CREATES THREE SIMUMLATED THETA CELL SPIKE TRAINS FROM COSINE VCOs
%%%WHICH WILL INTERFERE WITH ONE ANOTHER TO FORM THE GRID CELL IN FIG. 7F

load 'trackingdata';
pixpercm = 5; %pixels per cm in the tracking data

dX=diff(rsX); 
dY=diff(rsY); 
dX=[dX(1); dX]';
dY=[dY(1); dY]';
rsS=sqrt(dX.^2 + dY.^2)/(pixpercm*(rsTS(2)-rsTS(1))); %compute speed for modulating base frequency

%%%%% upsample tracking rate from 30 Hz to 500 Hz
tvec=interp1(rsTS, rsTS, [rsTS(1):.002:(rsTS(length(rsTS)))]); %rat path tracking timestamp vector 
pvec_x=interp1(rsTS, rsX, tvec); clear rsX; %rat path tracking X coordinate 
pvec_y=interp1(rsTS, rsY, tvec); clear rsY; %rat path tracking Y coordinate 
ratspeed=interp1(rsTS, rsS, tvec); clear rsS; %rat path tracking instantaneous speeds
clear rsTS; %delete 30 Hz timestamp vector


base_freq = 7; %VCO base frequency in Hz
speedslope = .025; %slope for linear dependence of base_freq upon running speed, in units of Hz/cm/s

%preferred vector lengths
rhovec = [11.2237   11.2237   11.2237];
%preferred vector orientations in radians
dirvec = [3.6652    5.7596    7.8540];
%VCO starting phases in radians    
startphases = [4.1888    5.2360    2.0944];
%shvec is a polar vector [r theta] specifying distance and direction by which to shift the envelope function ([0 0] for no shift)
shvec = [-4 60];
%adjust startphases to shift envelope as specified by shvec
startphases=startphases-[rhovec*shvec(1)*.0125].*cos([dirvec-deg2rad(shvec(2))]); 

%center the tracking data at the origin of the plane
center_x=(min(pvec_x)+max(pvec_x))/2;
center_y=(min(pvec_y)+max(pvec_y))/2;

%convert coordinate values from pixels to cm
pvec_x=(pvec_x-center_x)/pixpercm;
pvec_y=(pvec_y-center_y)/pixpercm;

%%%make sure position vectors are column vectors
if (size(pvec_x,1)==1)
    pvec_x=pvec_x';
end
if (size(pvec_y,1)==1)
    pvec_y=pvec_y';
end

dt=tvec(2)-tvec(1); %integration time step (assumed constant)

%%compute the speed-dependent base frequency (lfpfreq) and phase (lfpphase) at each time step
lfpfreq=ones(1,length(pvec_x)-1).*((base_freq+ratspeed(length(ratspeed)-1)*speedslope)*dt); 
lfpphase=cumsum([0 2*pi*lfpfreq]);

cullperc=.25; %theta cell spiking rate is determined by the percentage of "culled" spikes

for i=1:length(rhovec) %loop through VCOs
     pvec=[pvec_x pvec_y]*[[cos(dirvec(i)) sin(dirvec(i))];[-sin(dirvec(i)) cos(dirvec(i))]]; %pvec is rat's position along the VCO's preferred vector
     pvec=pvec(:,1); %drop y coordinates (x is position along the preferred vector)
     fvvec=diff(pvec)*dt; %speed along preferred vector
     foff=rhovec(i)*fvvec; %frequency of each VCO is offset from base_freq by a velocity-dependent value
     thfreq=lfpfreq'+foff;

    VCO(i).thphase=cumsum([startphases(i) [2*pi*(lfpfreq'+foff)]']);
    VCO(i).spiketrain=(((cos(VCO(i).thphase)+1).*(rand(size(VCO(i).thphase,1),size(VCO(i).thphase,2))))>.5); 
    ggg=find(VCO(i).spiketrain>0); ggh=rand(1,length(ggg)); ggg=ggg(find(ggh<cullperc)); 
    VCO(i).spiketrain(ggg)=0;
    VCO(i).spiketrain=tvec(find(VCO(i).spiketrain>0))-tvec(1);
    VCO(i).spiketrain=VCO(i).spiketrain(find(VCO(i).spiketrain>0));
    ggh=rand(1,length(VCO(i).spiketrain)); ggg=find(ggh<.5);
    VCO(i).spiketrain=VCO(i).spiketrain(ggg); %knock down to 50 Hz firing rate for theta cell
    if i==1
        mempot=cos(VCO(i).thphase);
    else
        mempot=mempot+cos(VCO(i).thphase);
    end
end


%save theta cell spike train time stamps (in units of milliseconds for NEURON)
 if (1)
     for i=1:length(rhovec)
         spkts0_r1=VCO(i).spiketrain'*1000; 
         save(['gridtheta_' num2str(i) '.dat'], 'spkts0_r1', '-ascii');
     end
 end

%plot graph approximating what NEURON simulation output should look like
threshlevel = 0.75;
spos=[];
 tout=10;
 thresh=max(mempot)*threshlevel;
 for t=1:length(tvec)
     if (mempot(t)>thresh)
         if (tout==10)
             spos=[spos; [pvec_x(t) pvec_y(t)]];
             tout=1;
         else
             tout=tout+1;
         end
     else
         tout=10;
     end
 end

figure(1);
hold off;
plot(pvec_x,pvec_y,'-k');
hold on;
scatter(spos(:,1),spos(:,2),'.r')
axis tight; axis square;
