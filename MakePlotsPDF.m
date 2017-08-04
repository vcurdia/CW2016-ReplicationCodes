% MakePlotsPDF
%
% Creates a Tex file, ready to load in Scientific Workplace with a
% collection of plots for each specification
%
% .........................................................................
% 
% Copyright 2008-2009 by Vasco Curdia and Michael Woodford
% Created: January 30, 2008
% Updated: July 30, 2010

%% ------------------------------------------------------------------------

%% Preamble
clear all
tic
setpath

bookmarkslevel = 1;
BenchmarkName = 'Pers_End';

Models = {'FF','NoFF','RepHH'};

nCP = 1;
% CmpPol(nCP).Policy = {'LQ','PiStab','Taylor'};
% CmpPol(nCP).Title = 'Optimal vs Inflation stab. vs Taylor';
% CmpPol(nCP).Label = 'LQTaylorPiStab';
CmpPol(nCP).Policy = {'LQ','PiStab','Taylor','FlexTarget','TaylorYn'};
CmpPol(nCP).Title = 'Optimal vs Inflation stab vs Taylor vs Targeting criterion vs TaylorYn';
CmpPol(nCP).Label = 'LQTaylorPiStabFlexTarget';
% nCP = nCP+1;
% CmpPol(nCP).Policy = {'LQ','Taylor','TaylorSP25','TaylorSP50','TaylorSP75','TaylorSP100'};
% CmpPol(nCP).Title = 'Taylor with spread';
% CmpPol(nCP).Label = 'TaylorSP';
% nCP = nCP+1;
% CmpPol(nCP).Policy = {'LQ','TaylorBm500','TaylorBm250','Taylor','TaylorBp250','TaylorBp500'};
% % CmpPol(nCP).Policy = {'LQ','TaylorBm050','TaylorBm025','Taylor','TaylorBp025','TaylorBp050'};
% CmpPol(nCP).Title = 'Taylor with credit';
% CmpPol(nCP).Label = 'TaylorB';
% nCP = nCP+1;
% CmpPol(nCP).Policy = {'LQ','FlexTarget','Taylor','TaylorRn','TaylorYn','TaylorRnYn'};
% CmpPol(nCP).Title = 'Taylor with nat vars';
% CmpPol(nCP).Label = 'TaylorNat';
% nCP = nCP+1;
% CmpPol(nCP).Policy = {'LQ','FlexTarget','Taylor','TaylorRnNoDFF','TaylorYnNoDFF','TaylorRnNoDFFYnNoDFF'};
% CmpPol(nCP).Title = 'Taylor with nat vars (NoDFF)';
% CmpPol(nCP).Label = 'TaylorNatNoDFF';
% nCP = nCP+1;
% CmpPol(nCP).Policy = {'LQ','FlexTarget','Taylor','TaylorRsBW','TaylorYsBW','TaylorRsBWYsBW'};
% CmpPol(nCP).Title = 'Taylor with nat vars (sBW)';
% CmpPol(nCP).Label = 'TaylorNatsBW';
% nCP = nCP+1;
% CmpPol(nCP).Policy = {'LQ','FlexTarget','Taylor','TaylorYn','TaylorYnNoDFF','TaylorYsBW'};
% CmpPol(nCP).Title = 'Taylor with output gap (different definitions)';
% CmpPol(nCP).Label = 'TaylorGapDD';
% nCP = nCP+1;
% CmpPol(nCP).Policy = {'LQ','FlexTarget','Taylor','TaylorRnYn','TaylorRnNoDFFYnNoDFF','TaylorRsBWYsBW'};
% CmpPol(nCP).Title = 'Taylor with Rn, Yn (different definitions)';
% CmpPol(nCP).Label = 'TaylorNatDD';
% nCP = nCP+1;
% CmpPol(nCP).Policy = ...
%     {'LQ','TaylorRnYn','TaylorRnYnSP25','TaylorRnYnSP50','TaylorRnYnSP75','TaylorRnYnSP100'};
% CmpPol(nCP).Title = 'Taylor with Rn, Yn and spread';
% CmpPol(nCP).Label = 'TaylorRnYnSP';
% nCP = nCP+1;
% CmpPol(nCP).Policy = ...
%     {'LQ','TaylorRnYnBm500','TaylorRnYnBm250','TaylorRnYn','TaylorRnYnBp250','TaylorRnYnBp500'};
% %     {'LQ','TaylorRnYnBm050','TaylorRnYnBm025','TaylorRnYn','TaylorRnYnBp025','TaylorRnYnBp050'};
% CmpPol(nCP).Title = 'Taylor with Rn, Yn and credit';
% CmpPol(nCP).Label = 'TaylorRnYnBm';
% nCP = nCP+1;
% CmpPol(nCP).Policy = ...
%     {'LQ','FlexTarget','FlexTargetNew50','FlexTargetNew100','FlexTargetNew150','FlexTargetNew200'};
% CmpPol(nCP).Title = 'Alternative Flexible Targeting Rule';
% CmpPol(nCP).Label = 'FlexTargetNew';
% nCP = nCP+1;
% CmpPol(nCP).Policy = ...
%     {'LQ','TaylorYn','TaylorYnSP25','TaylorYnSP50','TaylorYnSP75','TaylorYnSP100'};
% CmpPol(nCP).Title = 'Taylor with Yn and spread';
% CmpPol(nCP).Label = 'TaylorYnSP';
% nCP = nCP+1;
% CmpPol(nCP).Policy = ...
%     {'LQ','TaylorYnBW','TaylorYnBWSP25','TaylorYnBWSP50','TaylorYnBWSP75','TaylorYnBWSP100'};
% CmpPol(nCP).Title = 'Taylor with YnBW and spread';
% CmpPol(nCP).Label = 'TaylorYnBWSP';

Pers = 1;
End = 0:1;
NoRes = 0;
NoSpread = 0;
NoDist = 0;
GDebt = 0;
SmSigma = 0:1;
LowSigma = 0;

nM = length(Models);
nP = length(CmpPol(1).Policy);

%% set the specifications to consider
jE = 0;
for isPers=Pers
    for isEnd=End
        for isNoRes=NoRes
            for isNoSpread=NoSpread
                for isNoDist=NoDist
                    for isGDebt=GDebt
                        for isSmSigma=SmSigma
                            for isLowSigma=LowSigma
                                if isPers
                                    FileNameSuffix = '_Pers';
                                    TName = 'Plots:\\Persistence';
                                else
                                    FileNameSuffix = '_NoPers';
                                    TName = 'Plots:\\No Persistence';
                                end
                                if isEnd
                                    FileNameSuffix = sprintf('%s_End',FileNameSuffix);
                                    TName = [TName,', Endogenous Spread'];
                                else
                                    FileNameSuffix = sprintf('%s_Exo',FileNameSuffix);
                                    TName = [TName,', Exogenous Spread'];
                                end
                                if isNoRes
                                    FileNameSuffix = sprintf('%s_NoRes',FileNameSuffix);
                                    TName = [TName,', No Resources Used'];
                                end
                                if isNoSpread
                                    FileNameSuffix = sprintf('%s_NoSpread',FileNameSuffix);
                                    TName = [TName,', No Steady State Spread'];
                                end
                                if isNoDist
                                    FileNameSuffix = sprintf('%s_NoDist',FileNameSuffix);
                                    TName = [TName,', No Steady State Distortions'];
                                end
                                if isGDebt
                                    FileNameSuffix = sprintf('%s_GDebt',FileNameSuffix);
                                    TName = [TName,', Government Debt'];
                                end
                                if isSmSigma
                                    FileNameSuffix = sprintf('%s_SmSigma',FileNameSuffix);
                                    TName = [TName,', Same Intertemporal Elasticities'];
                                end
                                if isLowSigma
                                    FileNameSuffix = sprintf('%s_LowSigma',FileNameSuffix);
                                    TName = [TName,', High Intertemporal Elasticity'];
                                end
                                jE = jE+1;
                                ExerciseName{jE} = FileNameSuffix;
                                TitleName{jE} = TName;
                            end
                        end
                    end
                end
            end
        end
    end
end
nE = jE;

%% ------------------------------------------------------------------------

%% Loop over specifications
disp(' ')
for jE=1:nE
    TexFileName = sprintf('Plots%s',ExerciseName{jE});
    fprintf('%s\n',TexFileName)
%     load(['Output_FF',ExerciseName{jE}],'omega_ss')
    
%% Prepare sections
jS = 0;

%% Single plot for all models and policies
% for jM=1:1%nM %do it only for FF, we don't really care about the others...
%     load(sprintf('Output_%s%s',Models{jM},ExerciseName{jE}),'csi')
%     for jP=1:nP
%         jS = jS+1;
%         Section(jS).Title = sprintf('%s: %s',Models{jM},strrep(CmpPol(1).Policy{jP},'LQ','Optimal'));
%         Section(jS).Label = sprintf('%s_%s',Models{jM},CmpPol(1).Policy{jP});
%         csiP = csi;
%         if ismember(CmpPol(1).Policy{jP},{'LQ','PiStab','FlexTarget'})
%             [tf,idx] = ismember('xi_i',csi);
%             csiP(idx)=[];
%         end
%         Section(jS).csi = csiP;
%         Section(jS).Models = Models{jM};
%         Section(jS).Policy = CmpPol(1).Policy{jP};
%     end
% end

%% Set shock list
load(sprintf('Output_FF%s',ExerciseName{jE}),'csi')
csi = csi(1:end-2);
csi(end+1) = {'hCbar'};

%% Multiple models, single policy
for jP=1:nP
    jS = jS+1;
    Section(jS).Title = sprintf('%s: all models',CmpPol(1).Policy{jP});
    Section(jS).Label = sprintf('%s_AllMod',CmpPol(1).Policy{jP});
    csiP = csi;
    if ismember(CmpPol(1).Policy{jP},{'LQ','PiStab','FlexTarget'})
        [tf,idx] = ismember('xi_i',csi);
        csiP(idx)=[];
    end
    Section(jS).csi = csiP;
    Section(jS).Models = [Models{:}];
    Section(jS).Policy = CmpPol(1).Policy{jP};
end

%% FF, multiple policies
for jM=1%nM:-1:1
    for jCP=1:nCP
        if jCP>1&&~strcmp(Models(jM),'FF'), break, end
        if any(ismember({'TaylorSP50Bm50','TaylorBSPm50'},CmpPol(jCP).Policy))&&...
                ~isempty(strfind(ExerciseName{jE},'SmSigma'))
            continue
        end
        jS = jS+1;
        Section(jS).Title = sprintf('%s: %s',Models{jM},CmpPol(jCP).Title);
        Section(jS).Label = sprintf('%s_%s',Models{jM},CmpPol(jCP).Label);
        load(sprintf('Output_%s%s',Models{jM},ExerciseName{jE}),'csi')
        csiP = csi;
        if all(ismember(CmpPol(jCP).Policy,{'LQ','PiStab','FlexTarget'}))
            [tf,idx] = ismember('xi_i',csi);
            csiP(idx)=[];
        end
        Section(jS).csi = csiP;
        Section(jS).Models = Models{jM};
        Section(jS).Policy = [CmpPol(jCP).Policy{:}];
    end
end

%% set number of sections
nS = jS;

%% Begin tex file
fid=fopen(sprintf('%s.tex',TexFileName),'wt');
fprintf(fid,'\n\\documentclass[12pt]{article}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\usepackage{indentfirst}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage[dvipsnames,usenames]{color}\n');
fprintf(fid,'\\usepackage[pdftex,pdfstartview=FitH,bookmarksopen,');
    fprintf(fid,'bookmarksopenlevel=%.0f,colorlinks,linkcolor=MyDarkBlue,',bookmarkslevel);
    fprintf(fid,'citecolor=MyDarkBlue,urlcolor=MyDarkBlue,filecolor=MyDarkBlue,');
    fprintf(fid,'naturalnames]{hyperref}\n');
    
fprintf(fid,'\\usepackage[round,longnamesfirst]{natbib}\n');
fprintf(fid,'\\usepackage{fancyhdr}\n');

fprintf(fid,'\\hypersetup{pdftitle={%s},\n',TitleName{jE});
fprintf(fid,'            pdfauthor={Vasco Curdia}}\n');
fprintf(fid,'\\setlength{\\oddsidemargin}{0.0in}\n');
fprintf(fid,'\\setlength{\\evensidemargin}{0.0in}\n');
fprintf(fid,'\\setlength{\\topmargin}{0cm}\n');
fprintf(fid,'\\setlength{\\textheight}{8.5in}\n');
fprintf(fid,'\\setlength{\\textwidth}{6.5in}\n');
fprintf(fid,'\\setlength{\\hoffset}{0.0in}\n');
fprintf(fid,'\\setlength{\\headheight}{14.5pt}\n');
fprintf(fid,'\\pagestyle{fancy}\n');
fprintf(fid,'\\fancyhf{}\n');
fprintf(fid,'\\fancyhead[L]{\\leftmark}\n');
% fprintf(fid,'\\fancyhead[R]{\\textsc{%s}}\n',TexFileName);
fprintf(fid,'\\fancyfoot[C]{\\thepage}\n');
fprintf(fid,'\\renewcommand{\\headrulewidth}{0pt}\n');
fprintf(fid,'\\definecolor{MyDarkBlue}{rgb}{0,0.08,0.45}\n');
GraphicsPath = ['{',strrep(strrep(pwd,'\','/'),'D:','C:'),'/',strrep(ExerciseName{jE},'_',''),'/}'];
GraphicsPath = [GraphicsPath,strrep(GraphicsPath,'C:','D:')];
fprintf(fid,'\\graphicspath{%s}\n',GraphicsPath);

fprintf(fid,'\\begin{document}\n');

%% title page
fprintf(fid,'\\title{%s}\n',TitleName{jE});
fprintf(fid,'\\author{{Vasco C\\''{u}rdia}\\thanks{\\textit{E-mail address}: vasco.curdia@ny.frb.org} \\\\\n');
fprintf(fid,'{Federal Reserve Bank of New York} \\and\n');
fprintf(fid,'{Michael Woodford}\\thanks{\\textit{E-mail address}: michael.woodford@columbia.edu} \\\\\n');
fprintf(fid,'{Columbia University}}\n');
fprintf(fid,'\\maketitle\n');

fprintf(fid,'\\thispagestyle{empty}\n');
fprintf(fid,'\\newpage\n');

%% table of contents
fprintf(fid,'\\tableofcontents\n');
fprintf(fid,'\\newpage\n');

%% Insert plots in each section
for jS=1:nS
    fprintf(fid,'\\section{%s}\n',Section(jS).Title);
    csi = Section(jS).csi;
    ncsi = length(csi);
%     shock_size = ones(1,ncsi);
%     shock_size(ismember(csi,'hCsitil')) = (1+omega_ss)/omega_ss;
    for jcsi=1:ncsi
        csij = csi{jcsi};
        csijSection = strrep(csij,'_','');
        csijCaption = sprintf('%s_t',csij);
        csijCaption = strrep(csijCaption,'xi_i','\protect\xi^i');
        csijCaption = strrep(csijCaption,'hZ','Z');
        csijCaption = strrep(csijCaption,'hmu_w','\protect\mu^w');
        csijCaption = strrep(csijCaption,'htau','\protect\tau');
        csijCaption = strrep(csijCaption,'hG','G');
        csijCaption = strrep(csijCaption,'hb_g','\protect b^g');
        csijCaption = strrep(csijCaption,'hHbar','\protect\bar{H}');
        csijCaption = strrep(csijCaption,'hCbar_b','\protect\bar{C}^b');
        csijCaption = strrep(csijCaption,'hCbar_s','\protect\bar{C}^s');
        csijCaption = strrep(csijCaption,'hCbar','\protect\bar{C}');
        csijCaption = strrep(csijCaption,'hchitil','\protect\tilde{\chi}');
        csijCaption = strrep(csijCaption,'hXitil','\protect\tilde{\Xi}');
        csijCaption = strrep(csijCaption,'hbstar','\protect b^*');
        fprintf(fid,'\\subsection{Shock: %s}\n',csijSection);
        fprintf(fid,'\\begin{figure}[htbp] \\centering\n');
        fprintf(fid,'\\caption{Responses to ');
        if ismember(csij,{'htau'})
            fprintf(fid,'a shock to $%s$ of 1\\%%',csijCaption);
        elseif ismember(csij,{'hG'})
            fprintf(fid,'a shock to $%s$ equivalent to 1\\%% of steady state output',csijCaption);
        elseif ismember(csij,{'hb_g'})
            fprintf(fid,'a shock to $%s$ equivalent to 1\\%% of annual steady state output',csijCaption);
        elseif ismember(csij,{'hchitil','hXitil'})
            fprintf(fid,'a shock to $%s$ equivalent to 4\\%% increase in the spread',csijCaption);
        else
            fprintf(fid,'a 1\\%% shock to $%s$',csijCaption);
        end
        fprintf(fid,' (inflation and interest rate annualized)}');
        fprintf(fid,'\\label{%s_%s}\n',Section(jS).Label,csijSection);
        FigName = sprintf('IRF_%s_%s_%s.pdf',[Section(jS).Models,ExerciseName{jE}],Section(jS).Policy,csij);
        FigName = strrep(FigName,'Taylor','T');
        FigName = strrep(FigName,'FlexTarget','FT');
        FigName = strrep(FigName,'PiStab','Pi');
        fprintf(fid,'\\includegraphics[scale=1]{%s}\n',FigName);
        fprintf(fid,'\\end{figure}\n');
        fprintf(fid,'\\newpage \n');
    end
end

%% finish tex file
fprintf(fid,'\\end{document}\n');
fclose(fid);

%% Compile
% for j=1:3
%     system(['pdflatex ',TexFileName,'.tex']);
% end
% eval(['!',repmat([' pdflatex ',TexFileName,'.tex &'],1,3),' exit &'])
pdflatex(TexFileName)

%% end loop
end %jE

%% Compile all of them
% pdflatexall

%% count time
disp(' '),toc,disp(' ')

%% ------------------------------------------------------------------------
