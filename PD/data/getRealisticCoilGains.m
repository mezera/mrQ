% fig42_60_50
% fig55_48_40
% fig60_60_50
% 
% [K1,G1x,G1xx,G1y,G1yy,G1xy]
% 
% [ x.^2      2.*x.*y   2.*x.*z      2.*x        y.^2       2.*y.*z          2.*y     z.^2            2.*z             ones(size(x)) ]'
 % 1             2             3             4               5              6                     7          8               9                  10
STR={ 'k'     'X '   ' X^2 '    'Y '    'Y^2 '    'XY'     'Z'    ' Z^2 '    'XZ '    'YZ ' };

  fits={'/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/test/simulations/PolyOrder/fig60_60_50/Params_Box_Loc60_60_50_Size_42_42_42.mat'...
      '/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/test/simulations/PolyOrder/fig55_48_40/Params_Box_Loc55_48_40_Size_42_42_42.mat'...
      '/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/test/simulations/PolyOrder/fig42_60_50/Params_Box_Loc42_60_50_Size_42_42_42.mat'}
  for i=1:3
      load (fits{i})
params=zeros(10,32);
params(1,:)= Polf{2}.X( 10,:);
params(2,:)= Polf{2}.X(4,:);params(2,:)=params(2,:)*2;
params(3,:)= Polf{2}.X(1,:);

params(4,:)= Polf{2}.X(7,:);params(4,:)=params(4,:)*2;
params(5,:)= Polf{2}.X(5,:);

params(6,:)= Polf{2}.X(2,:);params(6,:)=params(6,:)*2;


params(7,:)= Polf{2}.X(9,:);params(7,:)=params(7,:)*2;

params(8,:)= Polf{2}.X(8,:);

params(9,:)= Polf{2}.X(3,:);params(9,:)=params(9,:)*2;

params(10,:)= Polf{2}.X(6,:);params(10,:)=params(10,:)*2;



C{i}.params=params;
C{i}.STR=STR;
C{i}.Phantomsource=fits{i};
  end


  
  save('CoilGains','C')