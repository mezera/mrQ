
%%
 error('say what?!')
 
run('/home/shai.berman/Documents/MATLAB/startup.m'); % display comments

t=0; times=zeros(1,20);
for sub=501:510
    subject=num2str(sub);   subject=strcat(subject,'_s');
    for sc=1:2
        t=t+1;
        scan=num2str(sc);   fol=strcat(subject,scan);
        
        if strcmp(fol, '501_s1') || strcmp(fol, '507_s2') % one is missing one was already done
            continue
        end
        
        outDir=strcat('/home/shai.berman/Documents/Code/mrQ_test/',fol,'/output');
        dir = strcat('/home/shai.berman/Documents/Code/mrQ_test/',fol,'/input');
        tic
        error('say what?!')
        command=sprintf('qsub -cwd -j y -b y -N job "matlab -nodisplay -r ''mrQ_run_Ver2(%s,%s); exit'' >log"', dir,outDir);
        system(command)

        times(t)=toc;

    end
end
