clear;

phz=NaN(12,6); %VCO phases are defined by a 12x6 matrix representing the output from a matrix of theta CPGs in Fig. 7C
rotation_angle = 0; %orientation of the tuning function is zero by default
q=0.65; % spike threshold is 0.65 by default

%%  Note that above, the phz matrix is assigned to all NaNs, meaning that the target neuron
%%  receives no input from any ring oscillator. Assigning a value of N=1-12 to an entry of the 
%%  phz matrix causes the target neuron to receive input from the theta cell at the Nth position
%%  in the ring oscillator at the corresponding matrix position. For example, phz(8,4)=10 causes 
%%  the target neuron to receive input from the 10th theta cell in the ring at row 8, column 4.
%%  Note that this scheme only permits the target neuron to receive input from one theta cell 
%%  per ring, and that by default, the weights of all inputs are uniformly equal to one.

%%%1) border cell against bottom edge in square box
% phz(10,:)=[3 4 6 8 12 10];

%%%2) large spacing entorhinal grid cell in square box 
%phz(:,4)=[10 NaN NaN NaN 10 NaN NaN NaN 4 NaN NaN NaN];

%%%3) small spacing entorhinal grid cell in square box 
phz(:,5)=[NaN 4 NaN NaN NaN 10 NaN NaN NaN 4 NaN NaN];

%%%4) single-field CA3 place cell in square box
% rotation_angle = 8*pi/6; %orientation of the tuning function is zero by default
% phz(:,3)=   [12    12     1     1     1     1     1     1    12    11    11    11];
% phz(:,4)=   [11    12     1     1     1     1     1     1    11    11    11    11];
% phz(:,5)=   [11    12     1     2     3     2     1     1    11    10     9    10];

%%%5) curved edge border cell in cylinder 
%   phz(12,:)=[2 2 3 5 8 2];
%   phz(1,:)=[2 2 3 5 8 2];
%   phz(2,:)=[2 2 3 5 8 2];
%   phz(3,:)=[2 2 3 5 8 2];
%   phz(4,:)=[2 2 3 5 8 2];
%   phz(5,:)=[2 2 3 5 8 2];
%   q=.75

%%%6) lumpy border cell against right edge in square box
% phz(1,:)=[2 3 4 6 9 4];
% phz(3,3:4)=[7 1]

%%%7) multi-field dentate place cell in square box 
%    phz(:,3)=   [ 9     1     2     2     8     5     3     8     9     9     5     3];
%    phz(:,4)=   [ 4     1     4     1     6     1     2     5     6    10    10    12];
%    phz(:,5)=   [ 1     3     3    11     8     1     1     6     7     3     7     5];

 
cells_per_ring = 12; %number of theta cells assumed to reside in each ring oscillator CPG
cell_phases=(0:(cells_per_ring-1))*2*pi/cells_per_ring; %phases for each cell in a VCO ring
phz=(phz./12)*2*pi; %convert phase matrix from cell numbers to VCO phases

figure(3);
clf
nrows=6; %number of rows in the CPG matrix
ncols=12; %number of columns in the CPG matrix
        Nmesh=50;
        MAX=5;
        ss=linspace(-MAX,MAX,Nmesh);
        [xx,yy]=meshgrid(ss,ss);
        
minrho=0.14; %smallest preferred vector length in the VCO matrix

cosnum=1; Esum=[];

for col=1:ncols
    theta(col,1:nrows)=pi+rotation_angle+(((col-1))/12)*2*pi; %assign the orientations of preferred vectors by column
    for row=1:nrows
         rho(col,row)=minrho*(sqrt(3)^(row-1)); 
        if isnan(phz(col,row)); %assign which cell in this ring projects to the target
            weight(col,row)=0; %weighting coefficient is 1 for the projection cell, implicitly 0 for all others
        else
            weight(col,row)=1; %weighting coefficient is 1 for the projection cell, implicitly 0 for all others
            px=cos(-theta(col,row));
            py=sin(-theta(col,row));
            %analytic signal
            s(cosnum,:,:)=exp(i*(rho(col,row)*px*xx + rho(col,row)*py*yy+phz(col,row)));
         if isempty(Esum);
             Esum=squeeze(s(cosnum,:,:));
         else
             Esum=Esum+squeeze(s(cosnum,:,:));
         end        
            cosnum=cosnum+1;
        end
    end
end

  
        Esum=abs(Esum);
  
        subplot(1,2,1);
        imagesc(ss,ss,Esum)
        axis equal
        axis tight
        colormap jet;
        set(gca,'xtick',[],'ytick',[]);
%colorbar;

        subplot(1,2,2);
        Emax=max(Esum(:));
        EE=max(0, Esum-q*Emax);
        EEmax=max(EE(:));

        imagesc(ss,ss,EE/EEmax)
        axis equal
        axis tight
        colormap jet;
        set(gca,'xtick',[],'ytick',[]);
