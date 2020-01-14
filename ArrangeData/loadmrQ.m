function mrQ = loadmrQ(mrQpath,saveStruct)

% INPUT: path to mrQ
% OUTPUT: mrQ structure with all of the updated paths 

%% check the stucture exists:

if notDefined('mrQpath')
    error('No mrQ path given as input')
elseif ~exist(mrQpath,'file')    
    error([mrQpath, ' does not exist'])
end

if notDefined('saveStruct')
    saveStruct=1;
    disp('resaving..')
end

%% load mrQ

load(mrQpath);
s=mrQ;

%% correct input directory

inOld = '/Mezer-lab';
inNew = '/ems/elsc-labs/mezer-a/Mezer-Lab';

s.RawDir = regexprep(s.RawDir,inOld,inNew);
s.inputdata_spgr.rawDir = regexprep(s.inputdata_spgr.rawDir,inOld,inNew);
s.inputdata_seir.rawDir = regexprep(s.inputdata_seir.rawDir,inOld,inNew);
if isfield(s,'InitSPGR_NtSp')
    for ll = 1:length(s.InitSPGR_NtSp.SPGR_niiFile);
        s.InitSPGR_NtSp.SPGR_niiFile{ll} = regexprep(  s.InitSPGR_NtSp.SPGR_niiFile{ll},inOld,inNew);
    end
end
if ~isempty(s.refIm)
    s.refIm  = regexprep(s.refIm,inOld,inNew);
end
%% check the analysis head directory

BaseDir = fileparts(fileparts(mrQpath));
mrQOutDir = fileparts(s.outDir);
if ~strcmp(mrQOutDir, BaseDir)
    disp('the general mrQ directory was wrong')
    flds = fields(s);
    for jj=1:length(flds)
        if ischar(s.(flds{jj}))  % if the field is a path
            s.(flds{jj}) = regexprep(s.(flds{jj}),mrQOutDir,BaseDir);
        elseif isstruct(s.(flds{jj})) % if it is a structure
            flds2 = fields(s.(flds{jj}));
            for ff = 1:length(flds2)
                if ischar(s.(flds{jj}).(flds2{ff}))
                    s.(flds{jj}).(flds2{ff}) = regexprep(s.(flds{jj}).(flds2{ff}),mrQOutDir,BaseDir);
                elseif iscell(s.(flds{jj}).(flds2{ff}))
                    for cc = 1:length(s.(flds{jj}).(flds2{ff}))
                        s.(flds{jj}).(flds2{ff}){cc} = regexprep(s.(flds{jj}).(flds2{ff}){cc},mrQOutDir,BaseDir);
                    end
                end
            end
        elseif iscell(s.(flds{jj})) % if it is a cell
            for cc = 1:length(s.(flds{jj}))
                if ischar(s.(flds{jj}){cc})
                    s.(flds{jj}){cc} = regexprep(s.(flds{jj}){cc},mrQOutDir,BaseDir);
                elseif isstruct(s.(flds{jj}){cc})
                     flds2 = fields(s.(flds{jj}){cc});
                     for ff = 1:length(flds2)
                         if ischar(s.(flds{jj}){cc}.(flds2{ff}))
                             s.(flds{jj}){cc}.(flds2{ff}) = regexprep(s.(flds{jj}){cc}.(flds2{ff}),mrQOutDir,BaseDir);
                         end
                     end
                end
            end
        end
    end
end

%% correct the head SPGR directory
d = dir([BaseDir,'/*/*/*/PD.nii.gz']);
if ~isempty(d)
ActualSPGRDir = d.folder;
mrqSPGRDir    = s.InitSPGR.spgr_initDir;

if ~strcmp(ActualSPGRDir, mrqSPGRDir)
    disp('the SPGR directory was wrong')
    flds = fields(s);
    for jj=1:length(flds)
        if strcmp(flds{jj},'InitSPGR_NtSp')
            continue
        end
        if ischar(s.(flds{jj}))  % if the field is a path
            s.(flds{jj}) = regexprep(s.(flds{jj}),mrQOutDir,BaseDir);
        elseif isstruct(s.(flds{jj})) % if it is a structure
            flds2 = fields(s.(flds{jj}));
            for ff = 1:length(flds2)
                if ischar(s.(flds{jj}).(flds2{ff}))
                    s.(flds{jj}).(flds2{ff}) = regexprep(s.(flds{jj}).(flds2{ff}),mrQOutDir,BaseDir);
                elseif iscell(s.(flds{jj}).(flds2{ff}))
                    for cc = 1:length(s.(flds{jj}).(flds2{ff}))
                        s.(flds{jj}).(flds2{ff}){cc} = regexprep(s.(flds{jj}).(flds2{ff}){cc},mrQOutDir,BaseDir);
                    end
                end
            end
% % %         elseif iscell(s.(flds{jj})) % if it is a cell
% % %             for cc = 1:length(s.(flds{jj}))
% % %                 if ischar(s.(flds{jj}){cc})
% % %                     s.(flds{jj}){cc} = regexprep(s.(flds{jj}){cc},mrQOutDir,BaseDir);
% % %                 elseif isstruct(s.(flds{jj}){cc})
% % %                      flds2 = fields(s.(flds{jj}){cc});
% % %                      for ff = 1:length(flds2)
% % %                          if ischar(s.(flds{jj}){cc}.(flds2{ff}))
% % %                              s.(flds{jj}){cc}.(flds2{ff}) = regexprep(s.(flds{jj}){cc}.(flds2{ff}),mrQOutDir,BaseDir);
% % %                          end
% % %                      end
% % %                 end
% % %             end
        end
    end
end
end
%% save 

mrQ=s;

if saveStruct
    save(mrQ.name,'mrQ');
end

