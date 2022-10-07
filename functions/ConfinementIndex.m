% Please, cite the article
% "APM_GUI: analyzing particle movement on the cell membrane and determining confinement"
% S.A. Menchon, M.G. Martin and C.G. Dotti,

% This script defines the confinement index.
% The original model is described in 
% "Detection of Temporary Lateral Confinement of Membrane Proteins Using Single-Particle Tracking Analysis"
% R. Simson, E.D. Sheets and K. Jacobson
% Biophysical Journal, 69: 989-993 (1995).

% A few modifications were applied based on:
% "Detection of confinement and jumps in single-molecule membrane trajectories"
% N. Meilhac, L. Le Guyader, L. Salome and N. Desteinville
% Physical Review E, 73: 011915 (2006).

% Almost the same as script Lforonepath.m from APM_GUI article
% Slightly modified by Manon by addition of extra output D
% Fix bug that last point of track was not in the analysis
% Also add the segments in the last part of the track, to get better
% detection at the end of the track
% Calculate running D in middle of segment intead of forward

function [temp, temp1, D] = ConfinementIndex(A,Sm,delta,long,npointsMSDset,Smin,Dinstant)

S=size(A,1);
A(S+1,:) = [NaN NaN A(S,3)+1];    % to allow indexing after last frame

% % It calculates instantaneous diffusion coefficient.  OLD WAY, 10 points
% forward
% npointsMSD=min(long-1,10);
% for j=1:S
%   if j>S-npointsMSD 
%     Dins(j)=MSqDin(A(S-npointsMSD:S,:),delta,4);
%   else
%     Dins(j)=MSqDin(A(j:j+npointsMSD,:),delta,4);
%   end
%   if Dins(j)<0
%     Dins(j)=0;
%   end
% end

% It calculates instantaneous diffusion coefficient. Middle point in window
npointsMSD=min(long-1,npointsMSDset);
hp = round(0.5 * npointsMSD);
for j=1:S
    if j<=hp
        Dins(j)=MSqDin(A(1:1+npointsMSD,:),delta,4);
    elseif j>S-hp 
        Dins(j)=MSqDin(A(S-npointsMSD:S,:),delta,4);
    else
        Dins(j)=MSqDin(A(j-hp:j-hp+npointsMSD,:),delta,4);
    end
    if Dins(j)<0
        Dins(j)=0;
    end
end

D = max(Dins);
% D = Dinstant;
% D = 0.05;
imin = Smin-1;

%%%%%It defines the last point to be considered for the analysis
k=0;

for i=1:S
  L(i,1)=delta*A(i,3);
  L(i,2)=0;
  NL(i)=0;
  if k==0
    if A(S,3)-A(S-i,3)>=Sm-1
      k=1;
      nmax=S-i; 
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate logpsi vector with NaN, to prevent mix of zeros coming from
% calculated zero and not filled in
logpsi = NaN(S-imin,((Sm-1)-imin+1));

%%%%%It defines the confinement index for each point in the trajectory until the last point considered for the analysis
for n=1:nmax
  k=0;
%%%%%It defines imax and imin (times could be non-consecutive because of the blinking) 
  for i=0:Sm
    if k==0
      slong=A(n+i,3)-A(n,3);
      if slong>Sm-1
        k=1;
        imax=i-1;
      end
    end
  end
%%%%%%%%%%
  NN(n)=imax;
  NM(n)=imax-imin+1;
%%%%%for segments with length from imin to imax
  for j=imin:imax
%%%%%It defines Rmax for each segment 
    Rmax=0;
    for i=1:j
      d=sqrt((A(n+i,1)-A(n,1))^2+(A(n+i,2)-A(n,2))^2); 
      if d>Rmax
        Rmax=d;
      end
    end
%%%%%%%%%%
    Ts=delta*(A(n+j,3)-A(n,3)); %Time between the initial an final point in the segment
    logpsi(n,j-imin+1)=0.2048-2.5117*D*Ts/Rmax^2; %Saxton's model
    if logpsi(n,j-imin+1)<=log10(0.1)
      logpsi(n,j-imin+1)=(-1)*logpsi(n,j-imin+1)-1;
    else
      logpsi(n,j-imin+1)=0;
    end
  end
%%%%%%%%%%
end
%%%%%%%%%%

%%%%%It defines the confinement index for points further than nmax
for n=nmax+1:S-imin
%%%%%It defines imax and imin (times could be non-consecutive because of the blinking) 
  imax = S-n;
%%%%%%%%%%
  NN(n)=imax;
  NM(n)=imax-imin+1;
%%%%%for segments with length from imin to imax
  for j=imin:imax
%%%%%It defines Rmax for each segment 
    Rmax=0;
    for i=1:j
      d=sqrt((A(n+i,1)-A(n,1))^2+(A(n+i,2)-A(n,2))^2); 
      if d>Rmax
        Rmax=d;
      end
    end
%%%%%%%%%%
    Ts=delta*(A(n+j,3)-A(n,3)); %Time between the initial an final point in the segment
    logpsi(n,j-imin+1)=0.2048-2.5117*D*Ts/Rmax^2; %Saxton's model
    if logpsi(n,j-imin+1)<=log10(0.1)
      logpsi(n,j-imin+1)=(-1)*logpsi(n,j-imin+1)-1;
    else
      logpsi(n,j-imin+1)=0;
    end
  end
%%%%%%%%%%
end
%%%%%%%%%%

for n=1:S-imin
  for i=0:NN(n)
    if i<Smin
      iniL=1;
    else
      iniL=i-(Smin-2);
    end
    for j=iniL:NM(n)
      L(n+i,2)=L(n+i,2)+logpsi(n,j);
      NL(n+i)=NL(n+i)+1;
    end
  end
end

for i=1:S
  if NL(i)>0
    L(i,2)=L(i,2)/NL(i); % confinement index
  end
end

temp=L;
temp1=Dins;


