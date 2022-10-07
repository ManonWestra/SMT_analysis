% Please, cite the article
% "APM_GUI: analyzing particle movement on the cell membrane and determining confinement"
% S.A. Menchon, M.G. Martin and C.G. Dotti,

% This file defines the confined regions using the confinement index defined in 'Lforonepath'.

% Some modifications by Manon Westra, 2020
% 
% OUPUT:
% temp  = Results: 
%           1. radius 2. x center 3. y center 4. dwell time (s)
%           5. start frame 6. end frame
% temp1 = Total dwell time confined & percentage of total time confined
% temp2 = Total time & critical L for this track


function [temp, temp1, temp2] = RfromL(A,L,alfa,Lcm,Tc,Lmax,Lcmax)

temp1(1,1)=0;

S=size(L);

temp2(1)=(L(S(1),1)-L(1,1));
Res=[];
%%%%%It defines the threshold, Lc
Lc=0;
nLc=0;
% calculate average L of the track
for i=1:S(1)
    Lc=Lc+L(i,2);
    nLc=nLc+1;
end
% multiply average L by factor alfa
Lc=alfa*Lc/nLc;
% when this Lc is lower than minimum Lc, Lc will be Lcm
if Lc<Lcm
  Lc=Lcm;
end

% extra option, when L is very high, Lcmax will be Lc
if max(L(:,2))>Lmax
  Lc=Lcmax;
end

temp2(2)=Lc;
%%%%%%%%%%%
% It defines the intervals where the confinement index is greater than the threshold
k1=0;
k2=1;
Tint(1,1)=0;
for i=1:S(1)
  if L(i,2)>Lc
    if k1==0
      % k2 is confinement zone number, this is the start frame
      Tint(k2,1)=i;
      k1=1;
    end
  else
    if k1==1
      % define the last frame of this k2 confinement zone
      Tint(k2,2)=i-1;
      k1=0;
      k2=k2+1;
    end
  end
end

if k1==1
  % when last frame was still confined, this is the end of confinement
  Tint(k2,2)=S(1);
end

% number of confinement zones
STint=size(Tint);

%It defines confined regions if the particle has a confinement index greater than or equal to the threshold in an interval of times longer than Tc (Tc is parameters(6) in 'explotting.m')
k1=0;
tDwell=0;
if Tint(1,1)~=0
  for j=1:STint(1)  % loop over confinement zones
    if (L(Tint(j,2),1)-L(Tint(j,1),1))>=Tc   % is time longer or equal to Tc? it's time NOT frames!
      tDwell=tDwell+L(Tint(j,2),1)-L(Tint(j,1),1);
      k1=k1+1;      % k1 is confinement zone number
      kT=0;         % ??? not used
      rcmx=0;
      rcmy=0;
      ncm=0;
      for i=Tint(j,1)+1:Tint(j,2)-1     % loop over all the confined frames
        ncm=ncm+1;
        rcmx=rcmx+A(i,1);               % sum all these x coordinates
        rcmy=rcmy+A(i,2);               % sum all these y coordinates
      end
      rcmx=rcmx/ncm;                    % center x
      rcmy=rcmy/ncm;                    % center y
      Rct=0;
      for i=Tint(j,1)+1:Tint(j,2)-1     % loop over all the confined frames 
        % define which coordinates are furthest away from center to define
        % radius of confinement zone
        d=sqrt((rcmx-A(i,1))^2+(rcmy-A(i,2))^2);
        if Rct<d
            Rct=d;
        end
      end
      Res(k1,1)=Rct;
      Res(k1,2)=rcmx;
      Res(k1,3)=rcmy;
      Res(k1,4)=L(Tint(j,2),1)-L(Tint(j,1),1);
      Res(k1,5)=Tint(j,1);
      Res(k1,6)=Tint(j,2);
    end
  end
end

if k1~=0
  temp1(1,1)=tDwell;
  temp1(1,2)=100*tDwell/(L(S(1),1)-L(1,1));
end  

temp=Res;