function [s2 Incoil] = mrQ_calculateCoilWeights(s1,opt,s,Incoil)
% 
% [s2 Incoil] = mrQ_calculateCoilWeights(s1,opt,s,Incoil)
% 
% Determine the channels with the best signal and place the data from those
% channels in 's2'. Those channels' indices are returned in 'Incoil' sorted
% in descending order of their mean (rate). *** VERIFY
% 
% INPUTS:   s1  - xformed/aligned multichannel data structure 
% 
%           opt - Determine the center of the boxs we will fit and create
%                 the grids and set grid options for the coil weighting
% 
% OUPUTS:
%           s2 - xformed/aligned data structure containing only data from
%                the top 8 channels (8 set by opt.numIn).  * Needs more
% 
%           Incoil - Coil indices. A MxN matrix where size(M) is the number of
%                    channels indexed and sorted in descending order by
%                    their mean for each grid element (N). *** Needs more
% 
% 
% WEB RESOURCES
%       http://white.stanford.edu/newlm/index.php/Quantitative_Imaging
% 
% 
% (C) Stanford University, VISTA Lab
% 
% 

%%

s2 = s; 
s2.imData = zeros(size(s.imData));

% Run over the grid you'd like to fit and check that there is data there
for i=1:numel(opt.X)
    
    [fb(1) fb(2) fb(3)] = ind2sub(size(opt.X),i);
    [Xx Yy Zz, skip]    = MrQPD_boxloc(opt,fb);
    
    
    if Xx(2) > opt.sz(1); Xx(2) = opt.sz(1);end;
    if Yy(2) > opt.sz(2); Yy(2) = opt.sz(2);end;
    if Zz(2) > opt.sz(3); Zz(2) = opt.sz(3);end;
    
    
    % *** Not sure about this code
    sz(1) = Xx(2) - Xx(1) + 1;
    sz(2) = Yy(2) - Yy(1) + 1;
    sz(3) = Zz(2) - Zz(1) + 1;
    sz(4) = size(s1,2);
    Dat   = zeros(sz);
    
    
    % For each of the channels populate Dat with imData
    for j=1:sz(4)      
        Dat(:,:,:,j) = s1(j).imData(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
    end;
    
    
    rate = squeeze(nanmean(nanmean(nanmean(Dat)))); % Rate is *** the non-zero mean for each channel in Dat
    SD   = squeeze(std(std(std(Dat))));             % The standard deviation of each channel *** Not used
    [~, In(:,i)] = sort(rate,'descend');            % 'In' is the index of sorted channel means in descending order
    
    
    % Assign the coil-weighted data to the 's2' structure.  ***
    % 'Incoil' will be the indices of the coils arranged in descending order
    % by their mean. This will only be set and used after the first volume
    % (ie. flip angle). 
    if ~exist('Incoil','var') || isempty(Incoil)   
        s2.imData(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2)) = sqrt(sum(Dat(:,:,:,In(1:opt.numIn,i)).^2,4));   
    else
        s2.imData(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2)) = sqrt(sum(Dat(:,:,:,Incoil(1:opt.numIn,i)).^2,4)); 
    end
    
end


if ~exist('Incoil','var') || isempty(Incoil)
    Incoil = In;
end


return




%% OLD CODE
% 
% Xx(1)=80;Xx(2)=100;
%         Yy(1)=135;Yy(2)=155;
%         Zz(1)=80;Zz(2)=100;
% 
% 
% D1=sum(Dat.^2,4);
% D2=sum(Dat(:,:,:,In(1:12)).^2,4);
% D3=sum(Dat(:,:,:,In(1:4)).^2,4);
% D3=zeros(size(D2));
% for i=1:sz(4);
% 
%     D3=D3+Dat(:,:,:,i)*SD(i)^2;
% end
% 
% 
% 
% receiver_weights1 =  receiver_weights(:,1).^(-2);
% receiver_weights2 =  receiver_weights(:,j).^(-2);
% coil_weights0=receiver_weights1./receiver_weights2;
% 
% endx=opt.X(end,end,end);
% endy=opt.Y(end,end,end);
% endz=opt.Z(end,end,end);
% for i=1:prod(size(opt{1}.X))
%     [fb(1) fb(2) fb(3)]=ind2sub(size(opt{1}.X),i);
%     Xx(1)=opt.X(fb(1),fb(2),fb(3))-opt.jump(1);
%     Xx(2)=opt.X(fb(1),fb(2),fb(3))+opt.jump(1);
%     Yy(1)=opt.Y(fb(1),fb(2),fb(3))-opt.jump(2);
%     Yy(2)=opt.Y(fb(1),fb(2),fb(3))+opt.jump(2);
%     Zz(1)=opt.Z(fb(1),fb(2),fb(3))-opt.jump(3);
%     Zz(2)=opt.Z(fb(1),fb(2),fb(3))+opt.jump(3);
%     if Zz(2)==endz, Zz(2)=opt.sz(3);end;
%     if Yy(2)==endy, Yy(2)=opt.sz(2);end;
%     if Xx(2)==endx, Xx(2)=opt.sz(1);end;
% 
%     box(:,:,:,:)=cat(4,s11(:).imData( Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2) ) );
%     noiseWeight=std(box,[],4);
% 
%     for recv=1:nrecv
%         sosrecdata = sosrecdata+abs(recdata(:,:,:,:,recv,:)) .^2 * coil_weights0(recv);
%     end
