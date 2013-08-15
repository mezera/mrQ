function  [ Segmask, C]= mrQ_localT1Seg(t1,T1cut,K)
%
%  [ Segmask C]= mrQ_localT1Seg(t1,T1cut)
% segmetation of a given T1 values (in secound) by there values to three
% regiond. the region are acourdind to the median with +-T1cut (0.5sec
% defult ). the T1cut will grow if most of the T1 value does not fall into
% three T1 semetation value.

if notDefined('T1cut')
        T1cut=0.1;
end
if notDefined('K')
        K=6;
end

      doneSeg=0;
      doneKmean=0
        Segmask=zeros(size(t1));
        
         [IDX,C] =kmeans(t1,K);
       %  cheack if the center are close. if so  we will use less clusters 
          while doneKmean==0
       for ii=1:K
        for jj=2:K
       if  abs(1-C(ii)/C(jj))<0.1
           Kless=Kless+1;
           
         if abs(1-C(1)/C(2))<0.1  &&  abs(1-C(1)/C(3))<0.1
             [IDX,C] =kmeans(t1,1);
         elseif abs(1-C(1)/C(2))<0.1
             [IDX,C] =kmeans(t1,2);
         elseif abs(1-C(1)/C(3))<0.1
             [IDX,C] =kmeans(t1,2);
         elseif abs(1-C(2)/C(3))<0.1
             [IDX,C] =kmeans(t1,2);
         end
         
         % take only the voxel that are close to the centers
        while doneSeg==0
        for ii=1:max(IDX)
        
        maskT= Segmask==0 & t1>(C(ii)*(1-T1cut) )  & t1<(C(ii)* (1+T1cut) );
        if  length(find(maskT))>100
        Segmask(maskT)=ii;
        end
        
        end
        
        % check if most of the values were segmented if not we will permit
        % greater distance form the center
        if  length(find(Segmask))>(length(Segmask)*0.6)
            doneSeg=1;
        else
            T1cut=T1cut*1.1;
                    Segmask=zeros(size(t1));
        end
        end
        