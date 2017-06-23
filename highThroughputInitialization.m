%% highThroughputInitialization
% Unfortunately these flowsheet problems require tuned initialization
% procedures. This script offers a method to interate many combinations of 
% initialization parameters for each optimization problem.

%% Setup

addpath('/usr/gams/24.0.2/');
%addpath('/usr/gams/24.0.1/');
addpath(pwd);

% clear all; close all;

% Debuging options
disableParFor = false;
testDegenHunter = false;
fixedNumberOfWorkers = -1;
loadGAMSSolutions = false;
skipWritingInitData = false;

runsToSkip = [ ];
% runsToSkip = setdiff(1:288,[2,223,33,49,65,229,202,129,211,192,182]);

% Setup threads
if(~disableParFor)
  matlabpool(8);
spmd
    if(~isdir(sprintf('worker%d', labindex)))
        mkdir(sprintf('worker%d', labindex));
    end
    cd(sprintf('worker%d', labindex));
end

else
    cd('worker1');
end

% Specify target parameters
initParams = {'NumStageInit',...  % 1
              'PhiLBCEOS',...
              'UBASFact',...
              'LBASFact',...
              'Tmin',... % 5
              'Prune',...
              'alpha',...
              'HRATSimple',...
              'ProbSpecParam1',...
              'ProbSpecParam2',... % 10
              'HRATPlusBD',...
              'FirstSolverSimple',...
              'FirstSolverCEOS',...
              'FirstSolverTray',...
              'HeatIntegModel',... % 15
              'FirstDerivEpsilon',...
              'Inter1Epsilon',...
              'CascInterLB',...
              'PlusStages', ...
              'FixEffForInit', ... % 20
              'SelectCEOSModel', ...
              'PruneConfig', ...
              'ProbSpecParam3', ...
              'HRATCEOS',...
              'SimpleCscInitProb',... %25
              'LoadInitPoint',...
              };
          
PosHtIntegModel = 15;

% Specify values for each param
initVals = cell(size(initParams));
% initVals{1} = [10, 15];

% Shortening
% initVals{1} = [15, 20, 25];
%initVals{1} = [10, 20];
initVals{1} = [15, 25];


% initVals{2} = [ -6, -7, -8, -9 ];
initVals{2} = [ -6, -7, -8];

initVals{3} = [10 ];

initVals{4} = [0.0001 ];

initVals{5} = [65];

initVals{6} = [0];

% initVals{7} = [-0.05, -0.1, -0.15];

% Shortening
initVals{7} = [-0.1, -0.15];
% initVals{7} = -0.15;

% initVals{8} = [2, 4];

% Shortening
initVals{8} = [4, 6];
% initVals{8} = 4;

initVals{9} = [0 ];
% initVals{9} = [0, -0.05, -0.1];
% initVals{9} = [0.7, 0.8, 0.9];

% initVals{10} = [0.8];
initVals{10} = [0.55];
% initVals{10} = [0.55, 0.65];

initVals{11} = [0];

initVals{12} = [1];
initVals{13} = [1];
initVals{14} = [1];
initVals{15} = 1;
% initVals{16} = [ -4, -5, -6 ];
initVals{16} = [-5, -6, -7];
% initVals{17} = [ -4, -5, -6 ];
initVals{17} = [-7];
% initVals{18} = [ -4, -5, -6 ];
initVals{18} = [-7];

initVals{19} = [10];
%initVals{19} = [10, 15];
initVals{20} = [0];
initVals{21} = [3];

% Shortening
initVals{22} = [0, 4];
% initVals{22} = [1];

% initVals{23} = 0.9825:0.0025:0.995;
% initVals{23} = 0.85:0.01:0.98;
initVals{23} = [0.95];

% parfor testing
% initVals{10} = 0.75;
% initVals{14} = 1;
% initVals{19} = 15;
% initVals{20} = 0;

initVals{24} = [1.5];

initVals{25} = [0, 1];
initVals{26} = [1];

% Testing Mode... only runs 3 jobs
% initVals{1} = 20;
% initVals{2} = [ -8];
% initVals{7} = -0.1;
% initVals{8} = 4;
% initVals{16} = -6;
% initVals{22} = 4;
% initVals{25} = 0;

% Specify results gdx files

resultsGDX = {'Sec1_SimpleFlowsheetResults',...
                'Sec2_InitCEOSResults',...
                'Sec3_CEOSFlowsheetResults',...
                'Sec4_MESH_FlowsheetResults',...
                'Sec5_HeatCheckerResults',...
                'Sec6_MESH_HtExDecomp_Results'};

% Specify active problems - not fully implemented but important
actProb = [1 1 1 1 1 1];

% Fields to parse
flds = {'o2rec','o2pure','Qs','Qw','liqPen','vapPen','Z9'};

% Extra GDX files to save
% saveFiles = {'matdata_connect','matdata_heatchckr','matdata_heatchckr2','matdata_tray','matdata_connect_final','matdata_tray_final'};
saveFiles = { } ;

%% Prepare

% Assemble Parameters Matrix
nParam = length(initParams);

nComb = 1;
nOptn = zeros(1,length(initVals));
for i = 1:length(initVals)
    nOptn(i) = length(initVals{i});
    nComb = nOptn(i)*nComb;
end

if(fixedNumberOfWorkers > 0)
    nComb = fixedNumberOfWorkers;
end

disp(['Considering ',num2str(nComb),' initialization combinations']);

paramM = zeros(nComb,nParam);
combM = paramM;
x = ones(1,nParam);
combM(1,:) = x;

for i = 2:nComb
    r = 1;
    j = nParam;
    while(r == 1)
        x(j) = x(j) + r;
        if(x(j) > nOptn(j))
            x(j) = 1;
        else
            r = 0;
        end
        j = j - 1;
    end
    combM(i,:) = x;
end

% Assemble Results Matrix
nStat = 3*sum(actProb);
nFlds = sum(actProb)*length(flds);
nPerProb = 3 + length(flds);
resultsM = zeros(nComb,nStat + nFlds);

tm = zeros(nComb,1);
wrk = zeros(nComb,1);
degenInfoMat = zeros(nComb,2);
degenInfoNames = cell(nComb);

for i = 1:nComb
    degenInfoNames{i} = cell(1);
end

bestObj = Inf;

%% Execute

parfor i = 1:nComb

%    if(disableParFor)
%        cd(sprintf('../worker%i',i));
%    end
    
    s = regexp(pwd,'worker(\d+)','tokens');
    wrk(i) = str2double(char(s{1}));
    
    % Load gdx file with selected parameters
    param = zeros(1,length(nOptn));
    for j = 1:length(nOptn)
        param(j) = initVals{j}( combM(i,j) );
    end
    paramM(i,:) = param;
    
    tmpCell = cell(1,length(param));
    tmp = struct();
    
    disp(' ')
%    disp(['Starting solve ',num2str(i),'/',num2str(nComb),' with'])
%    disp(['Starting to solve ',num2str(i),'/',num2str(nComb),' using worker ',num2str(labindex)])
%    drawnow('update');
    
    for j = 1:length(param)
        tmp.name = initParams{j};
        tmp.val = param(j);
        tmp.type = 'parameter';
        tmpCell{j} = tmp;
%        disp([initParams{j}, ' = ', num2str(param(j))]);
    end
    
    tmp.name = 'RunNumber';
    tmp.val = i;
    tmp.type = 'parameter';

    if(~skipWritingInitData)
    wgdx('initData',...
        tmp,...
        tmpCell{1}, ...
        tmpCell{2}, ...
        tmpCell{3}, ...
        tmpCell{4}, ...
        tmpCell{5}, ...
        tmpCell{6}, ...
        tmpCell{7}, ...
        tmpCell{8}, ...
        tmpCell{9}, ...
        tmpCell{10}, ...
        tmpCell{11}, ...
        tmpCell{12}, ...
        tmpCell{13}, ...
        tmpCell{14}, ...
        tmpCell{15}, ...
        tmpCell{16}, ...
        tmpCell{17}, ...
        tmpCell{18}, ...
        tmpCell{19}, ...
        tmpCell{20}, ...
        tmpCell{21}, ...
        tmpCell{22}, ...
        tmpCell{23}, ...
        tmpCell{24}, ...
        tmpCell{25}, ...
        tmpCell{26});
    end
    
    fail = false;
    obj = 0;
    tSav = 0;
    Sec1Stat = [];
    Sec2Stat = [];
    Sec3Stat = [];
    Sec4Stat = [];
    Sec5Stat = [];
    Sec6Stat = [];
    
    nNonPivot = [];
    nDegenEqn = [];
    namesDegEqns = [];

    % Solve optimization problem or load results
    
    if(~loadGAMSSolutions)
        tic
        try
            Sec1Stat = [];
            Sec2Stat = [];
            Sec3Stat = [];
            Sec4Stat = [];
            Sec5Stat = [];
            Sec6Stat = [];
            if(sum(i == runsToSkip))
              fail = true;
            else
%              disp('Attempting to run GAMS');
%              drawnow('update');
              [Sec1Stat, Sec2Stat, Sec3Stat, Sec4Stat, Sec5Stat, Sec6Stat] = gams('../OptimizeFlowsheet.gms');
%	      disp('Finished running GAMS');
%              drawnow('update');
            end
        catch
            fail = true;
        end
        t = toc;
        if(fail)
            disp(sprintf('Solved %i/%i \t\t Time %.2f seconds \t\t Obj: GAMS Failed',i,nComb,t));
            drawnow('update');
        end
%        disp(['Elapsed time: ',num2str(t),' seconds']);
        tm(i) = t;
        tSav = t;
    else
        gdxFile = 'gamsstat.gdx';
        tmp2.name = 'Sec1Stat';
        tmp2.form = 'sparse';
        Sec1Stat = rgdx(gdxFile,tmp2);
        tmp2.name = 'Sec2Stat';
        Sec2Stat = rgdx(gdxFile,tmp2);
        tmp2.name = 'Sec3Stat';
        Sec3Stat = rgdx(gdxFile,tmp2);
        tmp2.name = 'Sec4Stat';
        Sec4Stat = rgdx(gdxFile,tmp2);
        tmp2.name = 'Sec5Stat';
        Sec5Stat = rgdx(gdxFile,tmp2);
        tmp2.name = 'Sec6Stat';
        Sec6Stat = rgdx(gdxFile,tmp2);
        tm(i) = 0;
    end
    
    % Run Degeneracy Hunter
    % To Do: Modify Degeneracy Hunter to write to a file... can't rely on
    % evalc
    if(~fail)
        if(testDegenHunter)
            degenFile = '../';
        else
            degenFile = '../results/temp/';
        end
        try
%          disp('Attempting to run Degeneracy Hunter');
%          drawnow('update');
          [nNonPivot, nDegenEqn, namesDegEqns] = degeneracyHunter3(true,false,[degenFile,'DegenResults_',num2str(i),'.txt']);
%	  disp('Finished running Degeneracy Hunter');
        catch
          nNonPivot = -1;
          nDegenEqn = -1;
          namesDegEqns = 'Degen Hunter Failed!';
          disp('Degen Hunter Failed');
          drawnow('update');
        end
%    [T, nNonPivot, nDegenEqn, namesDegEqns] = evalc('degeneracyHunter3(true,false)');
    
 %       disp(['nNonPivot = ',num2str(nNonPivot)]);
 %       disp(['nDegenEqns = ',num2str(nDegenEqn)]);
 %       drawnow('update');
        try
%          disp(['Size degenInfoMat ',num2str(size(degenInfoMat))]);
          degenInfoMat(i,:) = [nNonPivot, nDegenEqn];
        catch
          degenInfoMat(i,:) = [-1, -1];
%          disp('Issue with degenInfoMat');
        end
%          drawnow('update');
        try
          degenInfoNames{i} = strjoin(namesDegEqns,'\t');
        catch
%          disp('Issue with strjoin(namesDegEqns)');
          degenInfoNames{i} = ' ';
        end
%        disp('Stored Degeneracy Hunter results');   
    else
        degenInfoMat(i,:) = NaN;
        degenInfoNames{i} = ' ';
    end
        
    % Parse solution fields
    
    tmpResultsM = zeros(1,nStat + nFlds);
    
    if(fail)
        resultsM(i,:) = NaN;
    else
        objSav = 0;
        m = 1;
        tmp2 = struct();
        stats = [];
        tmpGDX = [];
        tmpGDX2 = [];
        tmp2.form = 'full';
        for j = 1:length(resultsGDX)
            switch j
                case 1
                    stats = Sec1Stat.val(:,2)';
                case 2
                    stats = Sec2Stat.val(:,2)';
                case 3
                    stats = Sec3Stat.val(:,2)';
                case 4
                    stats = Sec4Stat.val(:,2)';
                case 5
                    stats = Sec5Stat.val(:,2)';
                case 6
                    stats = Sec6Stat.val(:,2)';
                otherwise
                    disp(['Warning: Something is wrong. Attempting to access results for Section ',num2str(j)]);
            end
            tmpResultsM(m:m+2) = stats;
%            tmpResultsM(i,((j-1)*nPerProb + 1):((j-1)*nPerProb + 2)) = stats;
            m = m + 3;
            
            for k = 1:length(flds)
                tmp2.name = flds{k};
%                switch tmpCell{PosHtIntegModel}.val
%                    case 1
                        %                        tmpGDX = [resultsGDX{j},'_U_p.gdx'];
                        tmpGDX = [resultsGDX{j},'.gdx'];
%                    case 2
%                        tmpGDX = [resultsGDX{j},'_R_p.gdx'];
%                end
                R = rgdx(tmpGDX,tmp2);
                tmpResultsM(m) = R.val;
%                tmpResultsM(i,(j-1)*nPerProb + 2 + k) = R.val;
                m = m + 1;
                if(k == length(flds) && j == length(resultsGDX))
%                    if(R.val < bestObj)
%                        bestObj = R.val;
%                    end
%                    display(['Tray-by-Tray Objective: ',num2str(R.val)]);
%                    display(['Best TbT Objective Thus Far: ',num2str(bestObj)]);

%                    disp(sprintf('Solved %i/%i \t\t Time %.2f seconds \t\t Obj: %e',i,nComb,t,R.val));
%                    drawnow('update');
                    objSav = R.val;
                end
                
            end
            
            % Copy GDX file for archieving
            switch tmpCell{PosHtIntegModel}.val
                case 1
                    tmpGDX2 = [resultsGDX{j},'_U_',num2str(i),'.gdx'];
                case 2
                    tmpGDX2 = [resultsGDX{j},'_R_',num2str(i),'.gdx'];
            end
            if(~loadGAMSSolutions)
              movefile(['./',tmpGDX],['../results/temp/',tmpGDX2]);
            end
            
        % resultsM... has previously before the "end"
        resultsM(i,:) = tmpResultsM;

        end
        
        for k = 1:length(saveFiles)
            tmpGDX = [saveFiles{k},'.gdx'];
            tmpGDX2 = [saveFiles{k},'_',num2str(i),'.gdx'];
            if(~loadGAMSSolutions)
              movefile(['./',tmpGDX],['../results/temp/',tmpGDX2]);
            end
        end
        
        estProg = 100*length(dir('../results/temp/Sec*'))/((length(resultsGDX))*nComb);
        disp(sprintf('Solved %i/%i \t\t Time %.2f seconds \t\t Obj: %e \t\t Est Progress %.1f %%',i,nComb,tSav,objSav,estProg));
        drawnow('update');
        
    end
    
    
    
end

if(~disableParFor)
  matlabpool close
end

%% Assemble Header File

hdr = cell(1,size(resultsM,2)+3);

m = 1;
for j = 1:length(resultsGDX)
    hdr(m:m+2) = {'Model','Solve','Timer'};
    m = m + 3;
    hdr(m:m + length(flds) - 1) = flds;
    m = m + length(flds);
end

hdr{end-2} = 'Time';
hdr{end-1} = 'NonPivot';
hdr{end} = 'NumberDegenEqns';

%% Create text log file
filename = 'hTI.csv';
fid = fopen(filename,'w');

fprintf(fid,'%s,','Run','Wrkr');

for i = 1:length(initParams)
    fprintf(fid,'%s,',initParams{i});
end

for i = 1:length(hdr)
    if(i ~= length(hdr))
        fprintf(fid,'%s,',hdr{i});
    else
        fprintf(fid,'%s\n',hdr{i});
    end
end
fclose(fid);

if(length(paramM) > nComb)
    paramM((nComb+1):end,:) = [];
    resultsM((nComb+1):end,:) = [];
    tm((nComb+1):end) = [];
end

dlmwrite(filename, [[1:nComb]', wrk, paramM, resultsM, tm, degenInfoMat], '-append','delimiter',',');

fid = fopen('DegenInfoSummary.txt','w');
fprintf(fid,'Run\tNumNonPivot\tNumDegenEqns\n');

for i = 1:nComb
    fprintf(fid,'%i\t%i\t%i',i,degenInfoMat(i,1),degenInfoMat(i,2));
    n = length(degenInfoNames{i});
    if(n > 1)
        fprintf(fid,'\t%s',degenInfoNames{i});
    end
    fprintf(fid,'\n');
end
fclose(fid);
