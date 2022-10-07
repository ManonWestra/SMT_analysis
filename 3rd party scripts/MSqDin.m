% Please, cite the article
% "APM_GUI: analyzing particle movement on the cell membrane and determining confinement"
% S.A. Menchon, M.G. Martin and C.G. Dotti,

% This script calculates the Mean Square Displacement (MSD), see
% "Direct observation of Brownian motion of lipids in a membrane"
% G.M. Lee, A. Ishihara and K.A. Jacobson

% In general we use 10-20 points of the paths to have the plot MSD vs. time but we only take 4 points to calculate the slope based on:
% "Single Particle Tracking: The distribution of Diffusion Coefficients"
% M.J. Saxton
% Biophysical Journal, 72: 1744-1753 (1997).

function temp1 = MSqDin(aA,adelta,anp)

ad=size(aA);

for ai=1:aA(ad(1),3)-aA(1,3)+2
  msd(ai)=0;
  v(ai)=0;
end;

for an=1:ad(1)-1
  for ai=1:ad(1)-an
    at=aA(ai+an,3)-aA(ai,3);
    msd(at+1)=msd(at+1)+((aA(ai+an,1)-aA(ai,1))^2+(aA(ai+an,2)-aA(ai,2))^2);
    v(at+1)=v(at+1)+1;
  end;
end;

for ai=1:aA(ad(1),3)-aA(1,3)+2
  if v(ai)~=0
  msd(ai)=msd(ai)/v(ai);
  end;
end;

at=0:adelta:(aA(ad(1),3)-aA(1,3)+1)*adelta;

ann=1;

msdnn(1)=0;
atnn(1)=0;

for ai=1:aA(ad(1),3)-aA(1,3)+1
  if v(ai+1)>0
      ann=ann+1;
      msdnn(ann)=msd(ai+1);
      atnn(ann)=at(ai+1);
  end;
end;

adnn=size(msdnn);

ab=anp;%5
if adnn(2)<ab
ab=adnn(2);
end;
p = polyfit(atnn(1:ab),msdnn(1:ab),1);
temp1=p(1)/4;


