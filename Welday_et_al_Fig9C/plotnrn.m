%This file generates a path plot for NEURON simulation results stored in:
%SPIKEOUT.dat

%Position tracking data is read in from the file:
%trackingdata.mat

spikes=load('-ascii','SPIKEOUT.dat');   %spike times from NEURON simulation
load 'trackingdata.mat';                %position data for path plot
spikes=spikes/1000;
sdex=[];
for i=1:length(spikes)
sdex=[sdex; find(rsTS>spikes(i),1,'first')];
end
figure(2); clf; hold off;
plot(rsX,rsY,'k'); axis square;
hold on;
scatter(rsX(sdex),rsY(sdex),'.r'); axis square;
set(gca,'XLim',[200 600],'YLim',[25 425]);