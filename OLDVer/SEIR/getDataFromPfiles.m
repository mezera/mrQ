% getDataFromPfile.mat
%
% Generates the .mat from the Pfiles. 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University    

% !!!NB: no homodyne reconstruction is done so far!!!

clear all
close all

T1path = '../';

savename = [T1path 'data/' 'SE3T_brain_022310_Pfiles'];

Pfile1 = '../Pfiles/P022310_1'
[im1,hd1] = reconX(Pfile1); 
TI1 = hd1.opti/1000;

Pfile2 = '../Pfiles/P022310_2'
[im2,hd2] = reconX(Pfile2); 
TI2 = hd2.opti/1000;

Pfile3 = '../Pfiles/P022310_3'
[im3,hd3] = reconX(Pfile3); 
TI3 = hd3.opti/1000;

Pfile4 = '../Pfiles/P022310_4'
[im4,hd4] = reconX(Pfile4); 
TI4 = hd4.opti/1000;

[x,y,z] = size(im3);

data = zeros(x,y,z,4);

% The signs below may need to be inverted
data(:,:,:,1) = im1;
data(:,:,:,2) = im2; 
data(:,:,:,3) = im3; 
data(:,:,:,4) = im4;

% data = abs(data);

extra.tVec = [TI1,TI2,TI3,TI4].';
extra.T1Vec = 1:5000;
save(savename,'data','extra')
