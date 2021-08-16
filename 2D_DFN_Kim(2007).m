% This program generates 2D fractal fracture networks by means of
% Multiplicative cascade process.
% Refering the article "Connectivity properties of two-dimensional fracture
% networks with stochastic fractal correlation".
% Input Units: [m] for DFN
% [mm] for aperture
% Internal calculation units: All unit will be transformed to [m]
% Aperture unit is transformed to [m] after exponential is taken.
%======================================================================
% VARIABLES
%======================================================================
% lratio l of the equation notation.
% FDc Fractal dimension of fracture center distribution.
% FDl Fractal dimension of fracture trace length distribution.
% It is the same as (a=FDl+1)
% L Side length of generation domain.
% alpha Fracture density term.
% lmin Min. length of fracture.
% nSet Number of fracture sets. Max. values is 3.
% If nSet = 0, then it generates random orientations.
% fracAnglei Orientation. i = 1~3. Angle measured from the x-axis
% in counterclockwise.
% FisherKi Fisher constant for nSet=1~3.
% OrienProbi Probability of each orientation. i=1~3
% q Multifractal index
% Aden Aperture density (ea/m)
% Amplitude Amplitude of self affine fractal
% H Hurst exponent
% Maxi Maximum value of the measured data
% Mini Minimum value of the measrued data
% STD STD value of the measured data
% bFlag Boundary effect flag. 1: on, 2: off
% niter number of iterations of the process
% nFrac Number of generated fractures.
% p() Probability of Multiplicative cascade process.
% if lratio=2, then 4 values will be generated.
% if lratio=3, then 9 values will be generated.
% Criteria Integer to chech probability values are properly
generated
% ranPos() Random position of probability
% pPar() Probability of parent cells
% pOffs() Probability of offspring cells
% inum Integer indexing iteration of loops
% jnum Integer indexing iteration of loops
% knum Integer indexing iteration of loops
% pCellNum Parent Cell Number
% offIndex indices of a offspring cell number of a reference parent cell
% ncol Cell numbers of offspring cells which are the no. 1 in
% their parent cells
% CellCoord The origin coordinates of each cell.
% FracProb Probability of fracture.
% cumProb Cummulative probability of cells.
% fracLoc Fracture locations according to proper probability.
% AngDev() Angular Deviation.
% FracData1() FracData1(a, b),
% 'a' refers fracture number
% 'b' refers data type. (fracture center coordinates).
% 1: x coordinate, 2: y coordinate
% 3: length, 4: orientaion
% FracData2() FracData2(a, b), - for Extended domain
% 'a' refers fracture number
% 'b' refers data type.
% 1: x coordinate of x1, 2: y coordinate of y1
% 3: x coordinate of x2, 2: y coordinate of y2
% FracData3() FracData3(a, b), - for Desired domain
% 'a' refers fracture number
% 'b' refers data type.
% 1: x coordinate of x1, 2: y coordinate of y1
% 3: x coordinate of x2, 2: y coordinate of y2
% 5: Orientation
% RoseData() Data for drawing a Rose Daigram.
% fLeng Fracture length after trimming fractures. This values are
% used for aperture generation.
% Maxlevel maximal number of recursions
% Sigma initial standard deviation
% Xval() Generated random aperture values
% delta() Array holding standard deviations
% Eps Acceptable error
% Ubound Upper boundary of the generation
% Lbound Lower boundary of the generation
% Range Boolean variable deciding whether Xval are in range
% 1 : True
% 2 : False
% Outbound1 Locations of data of which data are bigger than Ubound
% Outbound2 Locations of data of which data are smaller than Lbound
% Awidth Width of one aperture value
% fVol Fracture Volume [m^3/m]
% fAper() Aperture Distributioins of all fractures
% fPoro Fracture porosity (fVol/L^2)
%======================================================================
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**********************************************************************
% DFN Data
%**********************************************************************
lratio=input('Input Scale ratio, lratio : ', 's');
lratio=str2double(lratio);
FDc=input('Input Fractal dimension of center, FDc : ', 's');
FDc=str2double(FDc);
FDl=input('Input Fractal dimension of length, FDl : ', 's');
FDl=str2double(FDl);
L=input('Input Side Length of the Domain, L : ', 's');
L=str2double(L);
bFlag=input('Input Boundary Effect Flag, bFlag : ', 's');
bFlag=str2double(bFlag);
alpha=input('Input Fracture Density, alpha : ', 's');
alpha=str2double(alpha);
lmin=input('Input Min. Fracture Length, lmin : ', 's');
lmin=str2double(lmin);
q=input('Input Multifractal Index, q : ', 's');
q=str2double(q);
% If q==1, Poisson process like results will happen.
nSet=input('Input No. of Fracture Sets, nSet : ', 's');
nSet=str2double(nSet);
switch nSet
case 1
% Unit of input angle is in 'degree'.
fracAngle1=input('Input Frac. orientation, fracAngle1 : ',
's');
fracAngle1=str2double(fracAngle1);
FisherK1=input('Input Fisher Cosnt., FisherK1 : ', 's');
FisherK1=str2double(FisherK1);
OrienPorb1=input('Input Porb. of Orien., OrienPorb1 : ', 's');
OrienPorb1=str2double(OrienPorb1);
if OrienProb1>1
OrienProb1=1;
end % end of "if, OrienProb1"
case 2
fracAngle1=input('Input Frac. orientation, fracAngle1 : ',
's');
fracAngle1=str2double(fracAngle1);
FisherK1=input('Input Fisher Cosnt., FisherK1 : ', 's');
FisherK1=str2double(FisherK1);
OrienProb1=input('Input Porb. of Orien., OrienProb1 : ', 's');
OrienProb1=str2double(OrienProb1);
fracAngle2=input('Input Frac. orientation, fracAngle2 : ',
's');
fracAngle2=str2double(fracAngle2);
FisherK2=input('Input Fisher Cosnt., FisherK2 : ', 's');
FisherK2=str2double(FisherK2);
OrienProb2=input('Input Porb. of Orien., OrienProb2 : ', 's');
OrienProb2=str2double(OrienProb2);
if (OrienProb1+OrienProb2)>1
error('Probability should be less than 1')
end % end of "if, OrienProb1+OrienPorb2"
case 3
fracAngle1=input('Input Frac. orientation, fracAngle1 : ',
's');
fracAngle1=str2double(fracAngle1);
FisherK1=input('Input Fisher Cosnt., FisherK1 : ', 's');
FisherK1=str2double(FisherK1);
OrienProb1=input('Input Porb. of Orien., OrienProb1 : ', 's');
OrienProb1=str2double(OrienProb1);
fracAngle2=input('Input Frac. orientation, fracAngle2 : ',
's');
fracAngle2=str2double(fracAngle2);
FisherK2=input('Input Fisher Cosnt., FisherK2 : ', 's');
FisherK2=str2double(FisherK2);
OrienProb2=input('Input Porb. of Orien., OrienProb2 : ', 's');
OrienProb2=str2double(OrienProb2);
fracAngle3=input('Input Frac. orientation, fracAngle3 : ',
's');
fracAngle3=str2double(fracAngle3);
FisherK3=input('Input Fisher Cosnt., FisherK3 : ', 's');
FisherK3=str2double(FisherK3);
OrienProb3=input('Input Porb. of Orien., OrienProb3 : ', 's');
OrienProb3=str2double(OrienProb3);
if (OrienProb1+OrienProb2+OrienProb3)>1
    error('Probability should be less than 1')
end % end of "if, OrienProb1+OrienPorb2+OrienPorb3"
end % end of "switch, nSet"
%**********************************************************************
% Aperture Data
%**********************************************************************
Aden=input('Input Aperture density (ea/m) : ', 's');
Aden=str2double(Aden);
Amplitude=input('Input Amplitude : ', 's');
Amplitude=str2double(Amplitude);
H=input('Input Hurst Exponent : ', 's');
H=str2double(H);
Maxi=input('Input Maximum value : ', 's');
Maxi=str2double(Maxi);
Mini=input('Input Minimum value : ', 's');
Mini=str2double(Mini);
STD=input('Input STD value : ', 's');
STD=str2double(STD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Read Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch bFlag
case 1
L=2*L;
end % end of "switch, bFlag"
niter=round(0.5+log(L/lmin)/log(lratio));
nFrac=round(0.5+(alpha/FDl)*(L^FDc)/(lmin^FDl));
Criteria=1;
clear p;
switch q
case 1 % Weak fractal clustering
while Criteria==1
p(1:(lratio^2-1))=rand(1,(lratio^2-1))./lratio;
const=-FDc*log(lratio)-sum(p.*log(p));
p(lratio^2)=real(const/lambertw(const));
if sum(p.*log(p))==-FDc*log(lratio)
Criteria=2;
clear const
end % end of "if, sum(p.*log(p))"
end % end of "while, Criteria"
otherwise % Strong fractal clustering
while Criteria==1
p(1:(lratio^2-1))=rand(1,(lratio^2-1))./lratio^2;
proba=(1/lratio)^((q-1)*FDc)-sum(p.^q);
if proba>0
Criteria=2;
p(lratio^2)=proba^(1/q);
clear proba
end % end of "if, proba"
end % end of "while, Criteria"
end % end of "switch, q"
clear Criteria;
% Randomly permute probability for the first iteration
% Local cell ordering is followed:
% 1 | 3
% ---------------
% 2 | 4
ranPos=randperm(lratio^2);
pPar=p(ranPos);
% Generation of offspring cells and Assign probability
for iteration=2:niter
pCellNum=1;
for inum=1:lratio:(lratio)^iteration
iline=(inum+(lratio-1))/lratio;
ncol=(1+(iline-1)*lratio^(iteration+1):lratio:...
lratio^iteration+(iline-1)*lratio^(iteration+1));
switch lratio
case 2
offindex=[ncol, ncol+1, ncol+lratio^iteration, ...
ncol+lratio^iteration+1];
case 3
offindex=[ncol, ncol+1, ncol+2,
ncol+lratio^iteration, ...
ncol+lratio^iteration+1, ncol+lratio^iteration+2, ...
ncol+2*lratio^iteration,
ncol+2*lratio^iteration+1, ...
ncol+2*lratio^iteration+2];
end % end of "switch, lratio"
offindex=reshape(offindex, lratio^(iteration-1), lratio^2);
% Assigning prob. in y-direction (j-direction).
% jline refers row number in 'offindex'.
jline=1;
for jnum=1:lratio:(lratio)^iteration
linepos=offindex(jline,:);
ranPos=randperm(lratio^2);
pOffs(linepos)=pPar(pCellNum)*p(ranPos);
pCellNum=pCellNum+1;
jline=jline+1;
end % end of for "jnum"
end % end of for "inum"
pPar=pOffs;
pOffs=[];
end % end of for "iteration"
clear pOffs; clear offindex; clear ncol; clear linepos
xCoord=(-L/2:L/lratio^niter:L/2*(1-lratio^(-niter)))';
yCoord=(-L/2:L/lratio^niter:L/2*(1-lratio^(-niter)))';
for inum=1:lratio^niter
ncol=(1+(inum-1)*lratio^niter:inum*lratio^niter);
% CellCoord(cell no., coordinate), coordinate :1-(x), 2-(y)
CellCoord(ncol, 1)=xCoord;
CellCoord(ncol, 2)=yCoord(inum);
end % end of "for, inum"
cumProb=cumsum(pPar);
if cumProb(length(cumProb))>1
FracProb=rand(1, nFrac).*cumProb(length(cumProb));
else
FracProb=rand(1, nFrac);
end % end of "if, cumProb"
clear pPar;
for inum=1:lratio^(2*niter)
% inum refers cell number.
switch inum
case 1
fracLoc=find(FracProb>0 & FracProb<=cumProb(inum));
FracData1(fracLoc,1)=CellCoord(inum,1)+rand(size(fracLoc))*...
(L/lratio^niter);
FracData1(fracLoc,2)=CellCoord(inum,2)+rand(size(fracLoc))*...
(L/lratio^niter);
FracData1(fracLoc,3)=((alpha/FDl)*(L^FDc./fracLoc)).^(1/FDl);
otherwise
fracLoc=find(FracProb>cumProb(inum-1) &...
FracProb<=cumProb(inum));
FracData1(fracLoc,1)=CellCoord(inum,1)+rand(size(fracLoc))*...
(L/lratio^niter);
FracData1(fracLoc,2)=CellCoord(inum,2)+rand(size(fracLoc))*...
(L/lratio^niter);
FracData1(fracLoc,3)=((alpha/FDl)*(L^FDc./fracLoc)).^(1/FDl);
end % end of "switch, inum"
end % end of "for, inum"
row1=find(FracData1(:,3)>0);
FracData1=FracData1(row1,:);
nFrac=length(FracData1);
if cumProb(length(cumProb))>1
FracProb=rand(1, nFrac).*cumProb(length(cumProb));
else
FracProb=rand(1, nFrac);
end % end of "if, cumProb"
% Assigning fracture orientations.
switch nSet
case 0
FracData1(:,4)=rand(nFrac,1).*pi;
case 1
if OrienProb1~=1
FracProb=rand(nFrac,1);
fracLoc=find(FracProb<=OrienProb1);
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(fracLoc,4)=FracData1(fracLoc,4).*(pi/180);
fracLoc=find(FracProb>OrienProb1 & FracProb<=1);
FracData1(fracLoc, 4)=rand(size(fracLoc)).*pi;
else
AngDev=rad2deg(acos(log(1.-rand(nFrac,1))./FisherK1+1));
AngProb=rand(nFrac,1);
AngLoc=find(AngProb<=0.5);
FracData1(AngLoc,4)=fracAngle1+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
FracData1(AngLoc,4)=fracAngle1-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(:,4)=FracData1(:,4).*(pi/180);
end % end of "if, OrienProb1"
case 2
if (OrienProb1+OrienProb2)~=1
FracProb=rand(1,nFrac);
fracLoc=find(FracProb<=OrienProb1);
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(fracLoc,4)=FracData1(fracLoc,4).*(pi/180);
fracLoc=find(FracProb>OrienProb1 &...
FracProb<=(OrienProb1+OrienProb2));
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle2+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle2-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(fracLoc,4)=FracData1(fracLoc,4).*(pi/180);
fracLoc=find(FracProb>(OrienProb1+OrienProb2) &
FracProb<=1);
FracData1(fracLoc, 4)=rand(size(fracLoc)).*pi;
else
FracProb=rand(1,nFrac);
fracLoc=find(FracProb<=OrienProb1);
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
fracLoc=find(FracProb>OrienProb1 &...
FracProb<=(OrienProb1+OrienProb2));
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle2+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle2-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(:,4)=FracData1(:,4).*(pi/180);
end % end of "if, (OrienProb1+OrienProb2)"
case 3
if (OrienProb1+OrienProb2+OrienProb3)~=1
FracProb=rand(1,nFrac);
fracLoc=find(FracProb<=OrienProb1);
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(fracLoc,4)=FracData1(fracLoc,4).*(pi/180);
fracLoc=find(FracProb>OrienProb1 &...
FracProb<=(OrienProb1+OrienProb2));
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle2+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle2-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(fracLoc,4)=FracData1(fracLoc,4).*(pi/180);
fracLoc=find(FracProb>(OrienProb1+OrienProb2) &...
FracProb<=(OrienProb1+OrienProb2+OrienProb3));
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle3+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle3-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(fracLoc,4)=FracData1(fracLoc,4).*(pi/180);
fracLoc=find(FracProb>(OrienProb1+OrienProb2+OrienProb3)
&...
FracProb<=1);
FracData1(fracLoc, 4)=rand(size(fracLoc)).*pi;
else
FracProb=rand(1,nFrac);
fracLoc=find(FracProb<=OrienProb1);
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle1-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
fracLoc=find(FracProb>OrienProb1 &...
FracProb<=(OrienProb1+OrienProb2));
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle2+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle2-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
fracLoc=find(FracProb>(OrienProb1+OrienProb2) &...
FracProb<=(OrienProb1+OrienProb2+OrienProb3));
AngDev=rad2deg(acos(log(1.-
rand(size(fracLoc)))./FisherK1+1));
AngProb=rand(size(fracLoc));
AngLoc=find(AngProb<=0.5);
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle3+AngDev(AngLoc);
AngLoc=find(AngProb>0.5);
if length(AngLoc)>0
AngLoc2=fracLoc(AngLoc);
FracData1(AngLoc2,4)=fracAngle3-AngDev(AngLoc);
end % end of if, "length(AngLoc)>0"
FracData1(:,4)=FracData1(:,4).*(pi/180);
end % end of "if, (OrienProb1+OrienProb2+OrienProb3)"
end % end of "switch, nSet"
clear AngDev; clear AngProb; clear AngLoc; clear AngLco2;
% Calculate end points of fractures.
FracData2(:,1)=FracData1(:,1)+FracData1(:,3)./2.*cos(FracData1(:,4));
FracData2(:,2)=FracData1(:,2)+FracData1(:,3)./2.*sin(FracData1(:,4));
FracData2(:,3)=FracData1(:,1)-FracData1(:,3)./2.*cos(FracData1(:,4));
FracData2(:,4)=FracData1(:,2)-FracData1(:,3)./2.*sin(FracData1(:,4));
switch bFlag
case 1
%**********************************************************************
% Added for considering the Boundary Effect
%**********************************************************************
% Cut the sampling domain from the generation domain
L=L/2;
FracData3=FracData2;
% Identify x1 & x2 are greater than L/2
TrimIndex=find(FracData3(:,1)>L/2 & FracData3(:,3)>L/2);
FracData3(TrimIndex,:)=100*L;
% Identify x1 & x2 are less than -L/2
TrimIndex=find(FracData3(:,1)<-L/2 & FracData3(:,3)<-L/2);
FracData3(TrimIndex,:)=100*L;
% Identify y1 & y2 are greater than L/2
TrimIndex=find(FracData3(:,2)>L/2 & FracData3(:,4)>L/2);
FracData3(TrimIndex,:)=100*L;
% Identify y1 & y2 are less than -L/2
TrimIndex=find(FracData3(:,2)<-L/2 & FracData3(:,2)<-L/2);
FracData3(TrimIndex,:)=100*L;
%**********************************************************************
otherwise
FracData3=FracData2;
clear FracData2;
end % end of "switch, bFlag"
% Trim fracture of which ends are out of the domain.
% x1 & x2 > L/2
TrimIndex=find(FracData3(:,1)>L/2 & FracData3(:,1)~=100*L);
FracData3(TrimIndex,1)=L/2;
FracData3(TrimIndex,2)=tan(FracData1(TrimIndex,4)).*L/2+...
(FracData1(TrimIndex,2)-tan(FracData1(TrimIndex,4)).*...
FracData1(TrimIndex,1));
TrimIndex=find(FracData3(:,3)>L/2 & FracData3(:,3)~=100*L);
FracData3(TrimIndex,3)=L/2;
FracData3(TrimIndex,4)=tan(FracData1(TrimIndex,4)).*L/2+...
(FracData1(TrimIndex,2)-tan(FracData1(TrimIndex,4)).*...
FracData1(TrimIndex,1));
% x1 & x2 < -L/2
TrimIndex=find(FracData3(:,1)<-L/2);
FracData3(TrimIndex,1)=-L/2;
FracData3(TrimIndex,2)=tan(FracData1(TrimIndex,4)).*(-L/2)+...
(FracData1(TrimIndex,2)-tan(FracData1(TrimIndex,4)).*...
FracData1(TrimIndex,1));
TrimIndex=find(FracData3(:,3)<-L/2);
FracData3(TrimIndex,3)=-L/2;
FracData3(TrimIndex,4)=tan(FracData1(TrimIndex,4)).*(-L/2)+...
(FracData1(TrimIndex,2)-tan(FracData1(TrimIndex,4)).*...
FracData1(TrimIndex,1));
% y1 & y2 > L/2
TrimIndex=find(FracData3(:,2)>L/2 & FracData3(:,2)~=100*L);
FracData3(TrimIndex,2)=L/2;
FracData3(TrimIndex,1)=(L/2-FracData1(TrimIndex,2)+...
tan(FracData1(TrimIndex,4)).*FracData1(TrimIndex,1))./...
tan(FracData1(TrimIndex,4));
TrimIndex=find(FracData3(:,4)>L/2 & FracData3(:,4)~=100*L);
FracData3(TrimIndex,4)=L/2;
FracData3(TrimIndex,3)=(L/2-FracData1(TrimIndex,2)+...
tan(FracData1(TrimIndex,4)).*FracData1(TrimIndex,1))./...
tan(FracData1(TrimIndex,4));
% y1 & y2 < -L/2
TrimIndex=find(FracData3(:,2)<-L/2);
FracData3(TrimIndex,2)=-L/2;
FracData3(TrimIndex,1)=(-L/2-FracData1(TrimIndex,2)+...
tan(FracData1(TrimIndex,4)).*FracData1(TrimIndex,1))./...
tan(FracData1(TrimIndex,4));
TrimIndex=find(FracData3(:,4)<-L/2);
FracData3(TrimIndex,4)=-L/2;
FracData3(TrimIndex,3)=(-L/2-FracData1(TrimIndex,2)+...
tan(FracData1(TrimIndex,4)).*FracData1(TrimIndex,1))./...
tan(FracData1(TrimIndex,4));
switch bFlag
case 1
% Identify x1 & x2 are greater than L/2
TrimIndex=find(FracData3(:,1)>L/2 & FracData3(:,3)>L/2);
FracData3(TrimIndex,:)=100*L;
% Identify x1 & x2 are less than -L/2
TrimIndex=find(FracData3(:,1)<-L/2 & FracData3(:,3)<-L/2);
FracData3(TrimIndex,:)=100*L;
% Identify y1 & y2 are greater than L/2
TrimIndex=find(FracData3(:,2)>L/2 & FracData3(:,4)>L/2);
FracData3(TrimIndex,:)=100*L;
% Identify y1 & y2 are less than -L/2
TrimIndex=find(FracData3(:,2)<-L/2 & FracData3(:,2)<-L/2);
FracData3(TrimIndex,:)=100*L;
% Remove fractures which are not in the domain
TrimIndex=find(FracData3(:,1)~= 100*L);
FracData3=FracData3(TrimIndex,:);
FracData3(:,5)=FracData1(TrimIndex,4);
[nFrac bbb]=size(FracData3);
clear bbb;
otherwise
FracData3(:,5)=FracData1(:,4);
end % end of "switch, bFlag"
% Recalculate fracture length within the domain.
xLeng=FracData3(:,1)-FracData3(:,3);
yLeng=FracData3(:,2)-FracData3(:,4);
fLeng=sqrt(xLeng.^2+yLeng.^2);
clear xLeng; clear yLeng;
% Draw Fracture traces using hold function.
figure(1);
Xcoords=[FracData3(1,1), FracData3(1,3)];
Ycoords=[FracData3(1,2), FracData3(1,4)];
plot(Xcoords, Ycoords);
hold on;
for inum=2:nFrac
Xcoords=[FracData3(inum,1), FracData3(inum,3)];
Ycoords=[FracData3(inum,2), FracData3(inum,4)];
plot(Xcoords, Ycoords);
end % end of "for, inum"
hold off;
title('Fractal Discrete Fracture Networks');
axis([-L/2 L/2 -L/2 L/2]);
xlabel('L = ');
ylabel('L = ');
RoseData=[FracData3(:,5); FracData3(:,5)+pi];
figure(2);
rose(RoseData);
title('Rose Diagram of FDFN');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aperture Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iAper=1:nFrac
Maxlevel=round(0.5+log(Aden*fLeng(iAper))/log(2));
if Maxlevel<1
Maxlevel=1;
end % end of "if, Maxlevel"
Sigma=log(Amplitude/(3.86*10^(-4)*(2-H)^10.42));
Eps=10^-6;
N=2^Maxlevel;
Awidth=fLeng(iAper)/(N+1);
Xval=zeros(N+1, 1);
delta=zeros(Maxlevel,1);
index=(1:Maxlevel)';
Ubound=log(Maxi+STD*3);
if (Mini-STD)<=0
Lbound=log(Mini);
elseif (Mini-STD*2)<=0
Lbound=log(Mini-STD);
elseif (Mini-STD*3)<=0
Lbound=log(Mini-STD*2);
else
Lbound=log(Mini-STD*3);
end % end of "if, (Mini-STD)"
delta=Sigma*0.5.^(index.*H)*sqrt(0.5)*sqrt(1-2^(2*H-2));
Range=2;
while Range==2
Xval(N+1)=Sigma*randn;
if Xval(N+1)<=Ubound & Xval(N+1)>= Lbound
Range=1;
end % end of "if, Xval(N+1)"
end % end of "while, Range"
D=N;
d=N/2;
for level=1:Maxlevel
index=(d:D:N-d)';
index=index+1;
indexL=index-d;
indexR=index+d;
Xval(index)=0.5.*(Xval(indexL)+Xval(indexR));
index=(0:d:N)';
index=index+1;
Xval(index)=Xval(index)+delta(level).*randn(size(index));
Range=2;
while Range==2
Outbound=find(Xval<Lbound | Xval>Ubound);
if length(Outbound)==0
Range=1;
else
Xvalnew=zeros(size(Xval))+100*Ubound;
Xvalnew(Outbound)=Xval(Outbound)...
+delta(level).*randn(size(Outbound));
Inbound=find(Xvalnew>Lbound & Xvalnew<Ubound);
Xval(Inbound)=Xvalnew(Inbound);
clear Xvalnew;
end % end of "if, length(Outbound)==0"
end % end of "while, Range==2"
D=D/2;
d=d/2;
end % end of "for, level"
% Correct and generalize SRA algorithm
D=2;
d=1;
delta_new=delta(Maxlevel);
while (delta_new/Sigma) >= Eps
delta_new=delta_new*0.5^(0.5*H);
Xval=Xval+delta_new.*randn(size(Xval));
end % end of "while, delta_new/Sigma"
Xval=exp(Xval);
Xval=Xval./1000;
fVol(iAper)=sum(Xval.*Awidth);
lXval=length(Xval);
fAper(1:lXval,iAper)=Xval;
end % end of "for, iAper"
fPoro=sum(fVol)/L^2*100;

