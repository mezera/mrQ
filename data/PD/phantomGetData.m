function [M01 SZ meanVal ]= phantomGetData(boxsize,loc)

M0=readFileNifti( '/home/avivm/mrQ/PD/data/AligncombineCoilsM0.nii.gz');
pixdim=M0.pixdim(1:3);
M0=M0.data;


boxNeibhors = [-boxsize  -boxsize  -boxsize;  boxsize  boxsize    boxsize];


     
        if loc==1
            Cvoxel=[55 48 40];
        elseif loc==2
            Cvoxel=[55 48 30];
        elseif loc==3
            Cvoxel=[42 60 50];
        elseif loc==4
            Cvoxel=[60 60 50];
        elseif loc==5
            Cvoxel=[62 71 60];
        end
        %get te data
        clear M01 SZ meanVal XX YY ZZ
        [M01 SZ meanVal ]=GetLocalM0_1(M0,boxNeibhors,Cvoxel);
