function [Cbox,SHub] =mrQ_boxScaleGlobLinear(ScaleMat)
%[Cbox,SHub] =mrQ_boxScaleGlobLinear(ScaleMat) 
%
% We will find the scale of each box so it will agree with the values in
% the overlap boxes (a scalar term). We assume that the boxPD  are free
% from gain and the only difference is a scalar that needs to be fit. All
% the different scalars or each two boxes are the input ScaleMat. 
%                boxII= boxJJ*Ratio
%                scaleFactor(ii,jj)=1./Ratio;
%                scaleFactor(jj,ii)=Ratio;
%                
% We will write a minimization equation so all the scalars together will be
% adjusted. The idea is that box(i)-box(j)*scalar equals 0. We will build
% such a set of equations for all the overlap boxes and solve them together
% in the end.
%
% One box will be a reference box. In the case of the reference, we decide
% that scalar is 1.
%
% A simulation of the problem we solve can be found at SimLinBoxJoin.m
%
%
% AM (C) Stanford University, VISTA

%% Find the biggest network of connected boxes

% what is the typical ratio between two boxes?
MScale=median(ScaleMat(ScaleMat>0));
% If a ratio is more than twice the usual, it might be wrong. It might be
% better to get rid of those boxes. Typically this is only a small fraction
% then the box >>1%. This operation is useful for solving the linear
% system without outliers.
ScaleMat(ScaleMat>MScale*2)=0;
ScaleMat(ScaleMat<MScale*0.5)=0;

%  [S, C]= graphconncomp(double(sparse(logical(ScaleMat))));
 [S, C]= conncomp(double(sparse(logical(ScaleMat))));

 for ii=1:S
 NudeN(ii)= length(find(C==ii));
 end
  [N,ind]=sort(NudeN);
  
% All the boxes that are connected in the biggest network:
boxT0Use=(C==ind(end)) ;

% Let's save the clean Scale matrix that contains only the boxes we would
% like to work with.
ScaleMatNew=zeros(size(ScaleMat));
ScaleMatNew(boxT0Use,boxT0Use)=ScaleMat(boxT0Use,boxT0Use);

% Build the matrix for the scale linear calculation
LinScaleMat=zeros(size(ScaleMatNew));

%

% Calculate how many connections each box got
NConect=sum(logical(ScaleMatNew),1);
 Hub= (NConect==max(NConect)) & boxT0Use;

% Find the node that is highly connected and select one
 
 SHub=find(Hub);SHub=SHub(1,1);



%% Build the linear equations to solve for the XXX with respect to one box.

% We have removed the coil sensitivity from each of the boxes before
% getting here. The only difference between the boxes is the unknown scale
% factor, which we created above and stored in ScaleMat.

% Therefore, the relationship between the data in the boxes is, say,
%
%   box1 = s12*box2, and perhaps box1 = s1J*boxJ
%
% We do not have an sij for every box relationship, but we do have it for
% many and an interlocking set of boxes.
%
% So, to solve for the sij given the boxes, we can set up the linear
% equation
%
%   N*box1 = \sum_j s1j*boxj
%
% where the sum is over all the boxes that overlap with box 1.

% So, we can set up a linear equation that has a row for every box and a
%
% We expect that after scaling the entries in the box
for   ii=find(boxT0Use)
   % if (ii ~=SHub)
       
 %           boxII= boxJJ*Ratio
 %               scaleFactor(ii,jj)=1./Ratio;
 %              scaleFactor(jj,ii)=Ratio;

     
 LinScaleMat(ii,:)=-ScaleMatNew(ii,:);
 LinScaleMat(ii,ii)=NConect(ii);
%    end
end
LinScaleMat(end+1,SHub)=1;



%% Solve for y = [SIJ] * [BoxMean]
%
% 0 = [K, -s1J .... (Matrix)] * [boxMean1 .... boxMeanN]
%

BoxLocation=[find(boxT0Use)   ];

y=zeros(size(LinScaleMat,1),1);

y=y(BoxLocation);
Mat=LinScaleMat(BoxLocation,BoxLocation);
% Let's add one more equation that will make the Hub box to have a scale
% coefficient of one.
% 
Mat(end+1,:)=0;
Mat(:,end+1)=0;
Mat(end,find(BoxLocation==SHub))=1;
y(end+1)=1;

% Solve it as multi-linear equation
C=pinv(Mat'*Mat)*Mat'*y;
% The C we need is one over the one we fit


Cbox=zeros(length(BoxLocation),1);

Cbox(BoxLocation)=C(1:end-1);
end

%-----------------------------------------------------------------------%
% Originally, in line 34, the function "graphconncomp" was called. This
% function is part of the Bioinformatics Toolbox. Since this is the only
% call of the Bioinformatics Toolbox in all of mrQ, we wanted to rewrite
% the function so that users wouldn't need an additional toolbox just for
% one function.
% Fortunately, such a rewrite already exists. We use "conncomp" from Alec
% Jacobson's gptoolbox (https://github.com/alecjacobson/gptoolbox).
% 
% The following code is Copyright Alec Jacobson 2015.
%
% Edited 02-Dec-2015


function [S,C] = conncomp(G)
  % CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
  % toobox. G is an n by n adjacency matrix, then this identifies the S
  % connected components C. This is also an order of magnitude faster.
  %
  % [S,C] = conncomp(G)
  %
  % Inputs:
  %   G  n by n adjacency matrix
  % Outputs:
  %   S  scalar number of connected components
  %   C  
  [p,q,r] = dmperm(G+speye(size(G)));
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;
end

