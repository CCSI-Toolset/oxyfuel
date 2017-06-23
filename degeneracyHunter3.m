function [nNonPivot, nDegenSets, degenEqnsNames] = degeneracyHunter3(varargin)
%% ***** Degeneracy Hunter *****
% This script helps identify rank deficiencies (degeneracies) in the 
% Jacobian extracted from GAMS. This is particularly useful for model
% refinement of NLPs.
%
% Created by Alex Dowling at Carnegie Mellon University
% awdowlin@andrew.cmu.edu
% Last Updated: June 10th, 2014
%
% Number of inputs
% No inputs: act like a script
%
% Three inputs: Classic Degeneracy Hunter, MILP mode
% Input 1: Consider Weakly Active Constraints?
% Input 2: Consider Variable Bounds?
% Input 3: Filename to write output to ([ ] means display to console)

%% Toggle on/off core modules

% Run Degeneracy Hunter Module
module.degenHunt.tog = true;

% Run Second Order Condition Checker Module
module.SOCcheck.tog = false;

% Run Other Modules
module.other.tog = false;

% Run Setup Module
if(module.degenHunt.tog || module.SOCcheck.tog || module.other.tog)
  module.setup.tog = true;    
end

%% General Settings

% Should the code use sparse or dense linear algebra routines?
module.sparse = true;

% Tolerance for classifying multipliers
module.multTol = 1E-10;

% Tolerance for rank command
module.rankTol = 1E-10;

% Verbose output
module.verbose = false;

% Path to rank helper function
myrank = @(A) rankHelper(A,module.sparse,module.rankTol);

%% Degeneracy Hunter Module Specific Settings

% Should weakly active constraints be considered when checking for
% degeneracies?
module.degenHunt.weakAct = false;

% Should variables bounds be considered when checking for degeneracies?
module.degenHunt.varBounds = true;

% Use the heuristic search to identify degenerate equations
module.degenHunt.heurstic = false;

% Solve MILPs to identify degenerate equations
module.degenHunt.optimal = true;

% Tolerance for displaying equations in MILP degeneracy hunter
module.degenHunt.tol = 1E-6;

% Specify equations for Hueristic mode - only for debugging/special cases
% module.suspectEquations = [ ];
% module.suspectEqnGroup = [];


%% Second Order Condition Checker Settings

% Should weakly active constraints be considered when checking second
% order conditions?
module.SOCcheck.weakAct = true;

%% Other Module Specific Settings

% % Module 4 - This module identifies the largest and smallest elements of
% % the Jacobian.

module.other.four = true;
% % Order of magnitude for listing the largest (pos 1) and smallest (pos 2)
% % elements of the Jacobian
module.other.fourTol = [1, 6];


% % Module 5 - This module conducts analysis using the SVD of the Jacobian
% % to identify equations and variables contributing to ill-conditioning.

module.other.five = false;
% 
% module.fiveVerbose = true;
% 
% module.fiveCheckSmallSV = true;
% 
% % Module 6 - Print all non-zero elements for a specified variable or
% % equation

module.other.six = false;
module.other.sixSearch = {'F(SV2)'};
module.other.sixType = {'v'};

tocM = @() toc;
disp = @(s) fprintf('%s\n',s);

%% Parse Function Input
if(nargin == 3)
    module.degenHunt.tog = true;
    module.SOCcheck.tog = false;
    module.other.tog = false;
    
    module.degenHunt.heurstic = false;
    module.degenHunt.optimal = true;
    
    module.degenHunt.weakAct = varargin{1};
    module.degenHunt.varBounds = varargin{2};
    
    if(~isempty(varargin{3}))
        % Open file
        fid = fopen(varargin{3},'w');
        
        % Overload disp
        disp = @(s) fprintf(fid,'%s\n',s);
        
        % Redefine tocM
        tocM = @() disp(['Elapsed time is ',num2str(toc),' seconds.']);
    end
end

nNonPivot = -1;
nDegenSets = -1;
degenEqnsNames = [];

%% Setup Module
%%%%%%%%%%%%%%%

if(module.setup.tog)
    
    % % % % % % % % % % % % % % % % %
    % Classify equations and bounds %
    % % % % % % % % % % % % % % % % %
    
    % Extract equation names from GAMS
    s.name = 'i';
    eName = rgdx('hessian',s);
    
    nEq = length(eName.val);
    
    s.name = 'j';
    vName = rgdx('hessian',s);
    nVar = length(vName.val);
    
    % Extract equation multipliers and bounds from GAMS
    s.name = 'e';
    s.form = 'full';
    s.field = 'm';
    eqnMult = rgdx('hessian',s);
    
    s.field = 'lo';
    eqnLow = rgdx('hessian',s);
    
    s.field = 'up';
    eqnUp = rgdx('hessian',s);
    
    s.field = 'l';
    eqnLvl = rgdx('hessian',s);
    
    % Locate equality constraints
    eqlCnst = find(eqnUp.val(1:nEq) - eqnLow.val(1:nEq) < module.multTol);
    
    % Locate strongly active equality constraints
    eqlStrg = intersect( eqlCnst, find(abs(eqnMult.val(1:nEq)) > module.multTol) );
    
    % Locate weakly active equality constraints
    eqlWeak = setdiff(eqlCnst, eqlStrg);
    
    % Locate inequality constraints
    ineqlCnst = setdiff(1:nEq, eqlCnst);
    
    % Locate active inequality constraints
    actv = intersect(ineqlCnst, ...
        find(eqnLvl.val(1:nEq) - eqnLow.val(1:nEq) < module.multTol | eqnUp.val(1:nEq) - eqnLvl.val(1:nEq) < module.multTol ));
    
    % Locate strongly active inequality constraints
    actvS = intersect(actv, find(abs(eqnMult.val(1:nEq)) > module.multTol));
    
    % Locate weakly active inequality constraints
    actvW = setdiff(actv, actvS);
    
    % Determine sign on inequality constraints. Want h(x) <= 0
    ineqlSign = 1.0*(isfinite(eqnUp.val(1:nEq))) - 1.0*(isfinite(eqnLow.val(1:nEq)));
    
    % Extract variable multiplers and bounds from GAMS
    s.name = 'x';
    s.form = 'full';
    s.field = 'm';
    varMult = rgdx('hessian',s);
    
    s.field = 'lo';
    varLow = rgdx('hessian',s);
    
    s.field = 'up';
    varUp = rgdx('hessian',s);
    
    s.field = 'l';
    varLvl = rgdx('hessian',s);
    
    % Locate active variable bounds
    actvB = find(varLvl.val(nEq+1:end) - varLow.val(nEq+1:end) < module.multTol | varUp.val(nEq+1:end) - varLvl.val(nEq+1:end) < module.multTol );
    
    % Locate strongly active variable bounds
    actvBS = intersect(actvB, find(abs(varMult.val(nEq+1:end)) > module.multTol));
    
    % Locate weakly active variable bounds
    actvBW = setdiff(actvB, actvBS);
    
    % Determine which active variable bounds are lower
    actvBLow = intersect(actvB, find(varLvl.val(nEq+1:end) - varLow.val(nEq+1:end) < module.multTol));
    
    % Determine which active variables bounds are upper
    actvBUp = intersect(actvB, find(varUp.val(nEq+1:end) - varLvl.val(nEq+1:end) < module.multTol));
    
    % Note: With fixed variables there are two bounds and upper and lower.
    % Only one may be strongly active. Both may be weakly active. For the
    % purposes of null space calculations these will be treated as strongly
    % active upper bounds. This is equivalent to treating them like
    % equality constraints
    actvBBoth = intersect(actvBLow,actvBUp);
    
    actvBS = union(actvBS, actvBBoth);
    actvBUp = union(actvBUp, actvBBoth);
    actvBLow = setdiff(actvBLow, actvBBoth);
    
    % % Determine which active variable bounds are upper
    % actvBUp = setdiff(actvB, actvBLow);
    
    % % % % % % % % % % % % %
    % Assemble Jacobian(s)  %
    % % % % % % % % % % % % %
    
    clear s;
    
    % Extract the Jacobian from GAMS (GDX file)
    s.name = 'A';
    Jac = rgdx('hessian',s);
    
    % Rows: Equations, Columns: Variables
    i = Jac.val(:,1);
    j = Jac.val(:,2) - nEq;
    
    k = i > 0 & j > 0;
    
    J = sparse(i(k), j(k), Jac.val(k,3), nEq, nVar);
        
    % Assemble Jacobian matricies for bounds
    Jbnd_actv = zeros(length(actvB),nVar);
    
    Jbnd_Sactv = zeros(length(actvBS),nVar);
    
    Jbnd_Wactv = zeros(length(actvBW),nVar);
    
    if(~isempty(actvB) && (module.SOCcheck.tog || (module.degenHunt.tog && module.degenHunt.varBounds)))
        % Note: For the SOC Checker Module needs two bound Jacobians: one with
        % all active bounds. The other with only strongly active bounds.
        %
        % The DH Module bound matrix (active vs strongly active) depends on the settings
        
        % Note: This code won't work. Need to preserve variable order for DH module.
        %     countA = 0;
        %     countAS = 0;
        %     if(~isempty(actvBLow))
        %         for i = 1:length(actvBlow)
        %             Jbnd_actv(countA + i,actvBLow(i)) = -1;
        %             if(~isempty(intersect(actvBLow(i), actvBS)));
        %                 Jbnd_Sactv(countAS + i, actvBLow(i)) = -1;
        %             end
        %         end
        %         countA = countA + length(actvBLow);
        %         countAS = countAS + length(intersect(actvBLow, actvBS));
        %     end
        %     if(~isempty(actvBUp))
        %         for i = 1:length(actvBlow)
        %             Jbnd_actv(countA + i,actvBLow(i)) = 1; % Check this. Should this be negative one?
        %             if(~isempty(intersect(actvBLow(i), actvBS)));
        %                 Jbnd_Sactv(countAS + i, actvBLow(i)) = 1;
        %             end
        %         end
        %         countA = countA + length(actvBLow);
        %         countAS = countAS + length(intersect(actvBLow, actvBS));
        %     end
        
        count = 1;
        countW = 1;
        chckS = false;
        for i = 1:length(actvB)
            
            if(isempty(intersect(actvBS,actvB(i)))) % Bound is weakly active
                chckS = false;
            else % Bound is strongly active
                chckS = true;
            end
            
            if(~isempty(intersect(actvB(i), actvBLow))) % Lower bound
                Jbnd_actv(i,actvB(i)) = -1;
%                 if(isempty(intersect(actvB(i), actvBBoth))) % Fixed variable. Both upper and lower bounds
%                     if(varMult.val(actvB(i)) > module.multTol)
%                         % Lower bound is strongly active, upper bound is
%                         % weakly active
%                         chckS = true;
%                     end
%                 end
                if(chckS)
                    Jbnd_Sactv(count,actvB(i)) = -1;
                    count = count + 1;
                else
                    Jbnd_Wactv(countW,actvB(i)) = -1;
                    countW = countW + 1;
                end
            end
            if(~isempty(intersect(actvB(i), actvBUp))) % Upper bound
                Jbnd_actv(i,actvB(i)) = 1;
                
%                 if(isempty(intersect(actvB(i), actvBBoth))) % Fixed variable. Both upper and lower bounds
%                     if(varMult.val(actvB(i)) > module.multTol)
%                         % Lower bound is strongly active, upper bound is
%                         % weakly active
%                         chckS = false;
%                     end
%                 end
                
                if(chckS)
                    Jbnd_Sactv(count,actvB(i)) = 1;
                    count = count + 1;
                else
                    Jbnd_Wactv(countW,actvB(i)) = 1;
                    countW = countW + 1;
                end
            end
        end
    end
    
    % % % % % % % % % % % % %
    % Parse dictionary file %
    % % % % % % % % % % % % %
    % dict.txt is created by GAMS. It containts variable and equation names.
    
    % Read in file to string
    text = fileread('dict.txt');
    
    % Use regular expressions to parse string
    str = regexp(text,'  (x|e)(\d+)  (\S+)','tokens');
    
    % Manipulate to maintain similar structure as textscan results
    dic = cell(1,3);
    for i = 1:length(str)
        for j = 1:3
            dic{j}(i) = str{i}(j)';
        end
    end
    
    np = length(dic{2});
    
    % Create cell of equation and variable names
    eqn = cell(np,1);
    var = cell(np,1);
    eqnCount = 0;
    varCount = 0;
    
    for i = 1:np
        if(strcmpi(dic{1}(i), 'e'))
            eqn(eqnCount+1) = dic{3}(i);
            eqnCount = eqnCount + 1;
        elseif(strcmpi(dic{1}(i), 'x'))
            var(varCount+1) = dic{3}(i);
            varCount = varCount + 1;
        end
    end
    
    % Check eqnCount matches Hessian/Jacobian file(s)
    if(eqnCount ~= nEq)
        disp('Warning. Mismatch between number of equations is data');
        disp('and dictionary files');
    end
    if(varCount ~= nVar)
        disp('Warning. Mismatch between number of variables is data');
        disp('and dictionary files');
    end
    
    eqn = eqn(1:eqnCount);
    var = var(1:varCount);
    
    disp('***************************************');
    disp('********** Basic Information **********');
    disp('***************************************');
    disp(' ');
    disp('********** Equations **********');
    disp(['Total number: ', num2str(nEq)]);
    disp(' ');
    disp(['Number of equality constraints: ',num2str(length(eqlCnst))]);
    disp(['Number of STRONGLY active equality constraints: ',num2str(length(eqlStrg))]);
    disp(['Number of WEAKLY active equality constraints: ',num2str(length(eqlWeak))]);
    disp(' ');
    disp(['Number of inequality constraints: ',num2str(length(ineqlCnst))]);
    disp(['Number of ACTIVE inequality constraints: ',num2str(length(actv))]);
    disp(['Number of STRONGLY ACTIVE inequality constraints: ',num2str(length(actvS))]);
    disp(['Number of WEAKLY ACTIVE inequality constraints: ',num2str(length(actvW))]);
    disp(['Number of INACTIVE inequality constraints: ',num2str(length(ineqlCnst) - length(actv))]);
    disp(' ');
    disp('********** Variables **********');
    disp(['Total number: ',num2str(varCount)]);
    disp(['Number of ACTIVE variable bounds: ',num2str(length(actvB))]);
    disp(['Number of STRONGLY ACTIVE variable bounds: ',num2str(length(actvBS))]);
    disp(['Number of WEAKLY ACTIVE variable bounds: ',num2str(length(actvBW))]);
    disp(' ');
    disp('********** Degrees of Freedom **********');
    disp('At analyzed point, considering...');
    disp(['   all ACTIVE inequalities and bounds: ',num2str(varCount - length(eqlCnst) - length(actv) - length(actvB))]);
    disp(['   only STRONGLY ACTIVE inequalities and bounds: ',num2str(varCount - length(eqlCnst) - length(actvS) - length(actvBS))]);
    disp(' ');
end


%% Classic Degeneracy Hunter Module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(module.degenHunt.tog)
    
    disp('***********************************************');
    disp('********** Classic Degeneracy Hunter **********');
    disp('***********************************************');
    
    %% Part A: Rank Check
    % (1) Assemble the "active Jacobian": Adh
    % (1) Compute the rank of Adh
    % (3) Identify linearly dependent equations
    
    % Assemble the "active Jacobian" (Adh) and set of indicies (Adh_i) to map
    % back to equations and bounds
    
    if(module.degenHunt.weakAct)
        if(module.degenHunt.varBounds)
            Adh = [J(union(eqlCnst, actv),:); Jbnd_actv];
            Adh_i = union(union(eqlCnst,actv),actvB + nEq);
            disp('Analyzed Jacobian contains all active inequality constraints');
            disp('  and all active variable bounds');
        else
            Adh = J(union(eqlCnst, actv),:);
            Adh_i = union(eqlCnst, actv);
            disp('Analyzed Jacobian contains all active inequality constraints');
            disp('  and NO variable bounds');
        end
    else
        if(module.degenHunt.varBounds)
            Adh = [J(union(eqlCnst, actvS),:); Jbnd_Sactv];
            Adh_i = union(union(eqlCnst, actvS), actvBS + nEq);
            disp('Analyzed Jacobian contains only STRONGLY active inequality constraints');
            disp('  and only STRONGLY active variable bounds');
        else
            Adh = J(union(eqlCnst, actvS),:);
            Adh_i = union(eqlCnst, actvS);
            disp('Analyzed Jacobian contains only STRONGLY active inequality constraints');
            disp('  and NO variable bounds');
        end
    end
    
    disp(' ');
    
    if(sum(Adh_i > nEq) > 0)
        Adh_names = [eqn(Adh_i(Adh_i <= nEq)) ; var(Adh_i(Adh_i > nEq) - nEq)];
    else
        Adh_names = eqn(Adh_i);
    end
    
    nRowAdh = size(Adh,1);
    
    
    % Computer the rank of Adh
    r = myrank(Adh);
    disp(['Rank of analyzed Jacobian: ',num2str(r)]);
    
    nDE = size(Adh,1) - r;
    disp(['Estimated number of dependent equations (using rank): ',num2str(nDE)]);
    if(module.sparse && max(size(Adh)) < 1000)
        disp(['Estimated number of dependent equations (using dense rank calc.): ',num2str(size(Adh,1) - rank(full(Adh)))]);    
    end
    disp(' ');
    
    
    if(nDE > 0 || module.degenHunt.optimal)
        disp('********** Suspect Equations (and Bounds) **********');
        
        % Calculate the reduced row echelon form of the "small" Jacobian
        % This also gives the pivot points, which are the linearly independent
        % columns. Note: frref ignores the specified tolerance when operating
        % on a sparse matrix.
        tic
        [R, licol] = frref(Adh',1E-6,'s');
        tocM;
        
        % Determine the linearly dependent columns in the "small" Jacobian
        dcol = setdiff(1:size(Adh,1), licol);
        
        nNonPivot = length(dcol);
        
        if(isempty(dcol))
            disp('No degenerate equations detected during factorization.');
        else
            
            % Display the linearly dependent equation number and names
            st1 = [];
            st2 = [];
            disp(sprintf('Index\tActive\tType\tName\n'));
            for i = 1:length(dcol)
                if(sum(i == eqlCnst))
                    st2 = 'Eqlty.';
                    if(sum(i == eqlStrg))
                        str1 = 'Strong';
                    else
                        str1 = 'Weak';
                    end
                elseif sum(i == ineqlCnst)
                    st2 = 'Ineql.';
                else
                    str2 = 'V Bnd.';
                end
                if(sum(i == actvS) || sum(i == actvBS + nEq) || sum(i == eqlStrg))
                    st1 = 'Strong';
                else
                    st1 = 'Weak';
                end
                disp(sprintf('%i\t%s\t%s\t%s\n',Adh_i(dcol(i)),st1,st2,Adh_names{dcol(i)}));
            end
        end
    else
        dcol = [];
    end
    disp(' ');
    
    %% Part B: Heuristic Degeneracy Hunter
    % Parts 1 - 3
    % Overview - This module expands the set of suspect equations in search of
    % degeneracy. This is done iteratively.
    %
    % (1) Identify variables in suspect equations with non-zero entries in the
    %     Jacobian.
    % (2) Identify other equations with a non-zero entry in the Jacobian for
    %     the variables found in step 1. Add these equations to the suspect
    %     equation set if they are part of the active Jacobian (not inactive
    %     equality constraints.)
    % (3) If the expanded set of suspect equations is rank deficiency, stop.
    %     Other go to step 1.
    %
    % Part (4)
    % This module attempts to identify the smallest subset of the initial
    % suspect equations required to maintain rank decifiency.
    
    module.degenHunt.suspectEquations = [];
    
    module.degenHunt.suspectEqnGroup = [ ];
    
    if(module.degenHunt.heurstic)
        
        % Subparts 1 - 3
        % Process supplied suspect equations or get them from previous results
        if(isempty(module.degenHunt.suspectEquations))
            if(~isempty(dcol))
       %         problemEqns = Adh_i(dcol);
                problemEqns = dcol;
            else
                problemEqns = [];
            end
        else
            problemEqns = module.degenHunt.suspectEquations;
        end
        
        % How is problemEqns indexed?
        % Option 1: These are the indicies of Adh
        % Option 2: These are the indicies relating to all equations
        % Let's try option 1...
        
        if(isempty(problemEqns))
            if(module.verbose)
                disp('Warning: No linearly dependent equations. Skipping expansion step');
            end
        else
            
            flag = true;
            iter = 1;
            if(nDE ~= length(problemEqns) && isempty(module.degenHunt.suspectEquations))
                disp(['Due to numerical reasons the "rank" command predicted ',num2str(nDE),' linearly dependent equations']);
                disp(['whereas sparse QR factorization identified ',num2str(length(problemEqns)),]);
                disp(['For speed this module uses the "rank" command and searches for a rank deficiency of ',num2str(nDE)]);
                goalRankDef = nDE;
            else
                goalRankDef = length(problemEqns);
            end
            
            disp('Beginning to Expand Set of Suspect Equations (and Bounds).');
            
            while(flag)
                
                disp(['Iteration ',num2str(iter),' ...']);
                
                % Determine variables with non-zero entries for problem equations
                problemVar = zeros(varCount,1);
                for i = 1:length(problemEqns)
                    j = (Adh(problemEqns(i),:) ~= 0);
                    problemVar = logical(problemVar + j');
                    
                    if(module.verbose)
                        disp('Problem equation/bound...');
                        disp(eqn(problemEqns(i)));
                        disp('Variables with nonzero entries...');
                        disp(var(j));
                    end
                end
                disp(['Number of variables identified: ',num2str(sum(problemVar))]);
                
                % Need to expand this section to include logic for variable
                % bounds
                
                % Determine equations with non-zero entries for the identified variables
                pV = find(problemVar);
                suspectEqn = zeros(nRowAdh,1);
                
                for i = 1:length(pV)
                    k = (Adh(:,pV(i)) ~= 0);
                        suspectEqn = logical(suspectEqn + k);
                    
                    if(module.verbose)
                        disp('Identified variable...');
                        disp(var(pV(i)));
                        disp('Equations/bounds with nonzero entries for this variable');
                        disp(Adh_names(k));
                    end
                end
                
                nSE = sum(suspectEqn);
                disp(['Number of suspect equations (and bounds): ', num2str(nSE)]);
                
                % Check "reduced" Jacobian of only suspect equatons for rank deficiency
                redJac = Adh(suspectEqn,:);
                r = myrank(redJac');
                
                disp(['Rank: ',num2str(r)]);
                problemEqns = find(suspectEqn);
                iter = iter + 1;
                disp(' ');
                
                % Check if rank deficiency identified
                if(r + goalRankDef <= nSE)
                    flag = false;
                end
            end
            
            if(module.verbose)
                disp('Final Suspect Equations (and Bounds)');
                disp(eqn(problemEqns));
            end
        end
        
        % Subpart 4
        % Process supplied suspect equations or get them from Module 2 results
        if(~isempty(module.degenHunt.suspectEqnGroup))
            problemEqns = module.suspectEqnGroup;
        end
        
        if(isempty(problemEqns))
            if(module.verbose)
            disp('Warning: No linearly dependent equations. Skipping shrinking step');
            end
        else
            
            incldEqns = zeros(length(eqn),1);
            incldEqns(problemEqns) = 1;
            incldEqns = logical(incldEqns);
            
            nPE = sum(incldEqns);
            r = myrank(Adh(problemEqns,:)');
            flag = true;
            iter = 1;
            
            disp('Beginning to Shrink Set of Degenerate Equations (and Bounds)...');
            disp(['Number of equations (and bounds): ',num2str(nPE)]);
            disp(['Rank: ',num2str(r)]);
            disp(' ');
            
            while(flag)
                flag = false;
                nPE = sum(incldEqns);
                for i = 1:length(problemEqns)
                    incldEqns(problemEqns(i)) = false;
                    tr = myrank(Adh(incldEqns,:));
                    if(tr == r - 1)
                        % Remove equation
                        r = tr;
                        flag = true;
                    else
                        % Keep equation... add it back in
                        incldEqns(problemEqns(i)) = true;
                    end
                end
                nnPE = sum(incldEqns);
                disp(['Iteration ',num2str(iter),' ...']);
                disp(['New number of equation: ',num2str(nnPE)]);
                disp(['Number of equations dropped this iteration: ',num2str(nPE - nnPE)]);
                disp(['Final rank: ',num2str(r)]);
                disp(' ');
                
                problemEqns = find(incldEqns);
                iter = iter + 1;
            end
            
            disp('Smallest set of rank deficient equations...');
            disp(Adh_names(problemEqns));
            
            % Find multipliers for each equation
            
            miniJac = Adh(problemEqns,:);
            [mR, mLiCol] = frref(miniJac',1E-6,'s');
            
            % Determine the linearly dependent columns in the "mini" Jacobian
            mDCol = setdiff(1:size(miniJac,1), mLiCol);
            
            % Identify multipliers for each equation
            
            disp('***** Smallest Sets of Linearly Dependent Equations *****');
            for i = 1:length(mDCol)
                j = find(abs(mR(:,mDCol(i))) > 1E-6);
                disp(Adh_names(problemEqns(mDCol(i))));
                count = 1;
                for k = 1:length(j)
                    %               disp(eqn(problemEqns(j)));
                    
                    ii = problemEqns(j(k));
                    
                    if(sum(ii == eqlCnst))
                        st2 = 'Eqlty.';
                        if(sum(ii == eqlStrg))
                            str1 = 'Strong';
                        else
                            str1 = 'Weak';
                        end
                    elseif sum(ii == ineqlCnst)
                        st2 = 'Ineql.';
                    else
                        str2 = 'V Bnd.';
                    end
                    if(sum(ii == actvS) || sum(ii == actvBS + nEq) || sum(ii == eqlStrg))
                        st1 = 'Strong';
                    else
                        st1 = 'Weak';
                    end
                    
                    fprintf('%f\t%s\t%s\t%s\n',full(mR(j(k), mDCol(i))),st1,st2, eqn{problemEqns(j(k))});
                    count = count + 1;
                end
                disp(['Set Size: ',num2str(count),' Eqns (and Bounds)']);
                disp('***');
            end
        end
    end
    
    %% Part C: MILP Degeneracy Hunter
    
    if(module.degenHunt.optimal)
        
        % Setup function output
        nDegenSets = 0;
        degenEqnsNames = cell(1,nNonPivot);
        
        % Write active equations to GDX file
        actE.name = 'actE';
        actE.type = 'Set';
        actE.uels = cell(1,length(Adh_i));
        
        for i = 1:length(Adh_i)
            actE.uels{1,i} = ['e',num2str(Adh_i(i))];
        end
        
        for j = 1:length(dcol)
            SE.name = 'SE';
            SE.type = 'Set';
            SE.uels = {['e',num2str(Adh_i(dcol(j)))]};
            wgdx('degenData',actE,SE);
            
            disp(['Consider degeneracy with ',Adh_names{dcol(j)}]);
            tic
            gamso.form = 'full';
            [x, y, Z, MIPstat] = gams('../GenerateMinDegenerateSet_MIP2.gms');
            tocM;
            switch MIPstat.val
                case 1
                    disp(['Optimal: This equation is part of a degenerate set with ',num2str(Z.val),' equations (total)']);
%                    n = find(y.val(Adh_i) > 0);
                    n = find(abs(x.val(Adh_i)) > module.degenHunt.tol);
                    
                    for k = 1:length(n)
                        disp(sprintf('%f \t\t %s',x.val(Adh_i(n(k))),Adh_names{n(k)}));
                    end
                    if(length(n) < round(Z.val))
                        disp(['Only ',num2str(length(n)),' equation displayed. The remaining ', num2str(Z.val - length(n)),' equations`']);
                        disp(['singular vector components are smaller than ',num2str(module.degenHunt.tol),' in magnitude.']);
                    end
                    
                    nDegenSets = nDegenSets + 1;
                    degenEqnsNames{nDegenSets} = Adh_names{dcol(j)};
                otherwise
                    disp('Infeasible: This equation is not part of a degenerate set');
            end
            
            disp(' ');
                
        
        end
        
        if(nDegenSets < nNonPivot)
            degenEqnsNames = degenEqnsNames(1:nDegenSets);    
        end
        
    end
end

%% Second Order Condition Checker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(module.SOCcheck.tog)
    
    disp('*****************************************************');
    disp('********** Second Order Conditions Checker **********');
    disp('*****************************************************');
    
    
    % Calculate null space
    
    if(module.SOCcheck.weakAct)
        A2 = [J(union(eqlCnst, actv),:); Jbnd_actv];
    else
        A2 = [J(union(eqlCnst, actvS),:); Jbnd_Sactv];
    end
    
    Z = nulls(A2);
%    Z = null(full(A2));
    
    % Extract objective coefficient (min or max) from GAMS
    s.name = 'objcoef';
    objcoef = rgdx('hessian',s);
    
    % Extract the Hessian from GAMS
    s.name = 'W';
%    s.form = 'sparse';
    Hes = rgdx('hessian',s);

    i = Hes.val(:,1) - nEq;
    j = Hes.val(:,2) - nEq;
    
    k = i > 0 & j > 0;
    
    W = sparse(i(k), j(k), Hes.val(k,3), nVar, nVar);
       
%    W = sparse(Hes.val(:,1),Hes.val(:,2),Hes.val(:,3));
    
%    W = W(nEq+1:end,nEq+1:end);
    
    % This is needed b/c only an upper triangle is extracted from GAMS
    W = objcoef.val*(W + W' - diag(diag(W)));
    
    if(isempty(Z))
        disp('The null space is empty, and thus');
        disp('second order conditions are vacuously true.');
    else
        
        Wr = Z'*W*Z;
        
        [V, D] = eigs(Wr);
        
%        Z = nulls(A2);
        
        d = diag(D);
        i = find(d < 0);
        
        disp('Smallest Eigenvalue of Reduced Hessian');
        disp(min(real(d)));
        
        disp('Largest Eigenvalue of Reduced Hessian');
        disp(max(real(d)));
        
        if(~isempty(i))
            if(isempty(actvW) && isempty(actvBW))
                disp('This point violates second order conditions.');
                disp('There are no weakly active constraints or bounds');
                disp('and the reduced Hessian has a negative eigenvalue.');
            else
                
                disp('This point may violate second order conditions.');
                disp(' ');
                disp('Beginning to search for an allowable direction with');
                disp('negative curvature...');
                disp(' ');
                
                if(isempty(actvW) && ~isempty(actvBW))
                    A2w = Jbnd_Wactv;
                elseif(~isempty(actvW) && isempty(actvBW))
                    A2w = diag(ineqlSign(actvW))*J(actvW,:);
                else
                    A2w = [diag(ineqlSign(actvW))*J(actvW,:); Jbnd_Wactv];
                end
            end
            
            f = zeros(length(Wr),1);
            b = zeros(size(A2w,1),1);
            x0 = V(:,real(d) == min(real(d)));
            bnd = ones(size(x0));
            
            optns = optimset('algorithm','active-set');
            
            [x, fval] = quadprog(full(Wr),f,full(A2w*Z),b,[],[],-bnd,bnd,x0,optns);
            
            xFull = Z*x;
            
            if(fval < 0)
                disp('Second Order Conditions (SOC) are violated. Located');
                disp('a feasible direction with negative curvature:');
                disp(' ');
                for k = 1:nVar
                    if(xFull(k) ~= 0)
                    disp(sprintf('%s \t%f',var{k},full(xFull(k))));
                    end
                end
            else
                disp('Inconclusive. Failed to find feasible direction with');
                disp('negative curvature. Need to solve non-convex QP to');
                disp('global optimality to prove second order conditions are');
                disp('satified.');
            end
        else
            disp('There are no negative eigenvalues of the reduced Hessian');
            disp('at this point, thus second order conditions are satisfied.');
            
% The LP version does not work!!!!            
%             % Time to solve LPs
%             np = find(real(d) <= 0); % Locate non-positive eigenvalues
%             f = zeros(length(np)+1,1);
%             f(end) = 1;
%             Alp = [A2w(np,:)*Z*V(np,:), -1*ones(length(np),1)];
%             blp = zeros(1,size(A2w,1));
%             
%             beq = 1;
%             
%             % Verify this tomorrow
%             flag = true;
%             i = 1;
%             optn = optimset('Display','off');
%             while(flag)
%                 
%                 Aeq = zeros(1,length(np)+1);
%                 Aeq(ceil(i/2)) = 1;
%                 beq = 2*mod(i,2) - 1;
%                 
%                 v_lo = [-Inf(length(np),1); -1];
%                 v_up = [ ];
%                 
%                 vWght = linprog(f,Alp,blp,Aeq,beq,v_lo,v_up,f,optn);
%                 
%                 deltaX = Z*V(np,:)*vWght(1:end-1);
%                 
%                 if(vWght(end) < 0)
%                     disp('Allowable negative curvature direction found.')
%                     disp('Second order conditions violated.');
%                     disp('Direction stored in "deltaX"');
%                     flag = false;
%                 else
%                     i = i + 1;
%                     if(i > 2*length(np))
%                         flag = false;
%                         disp('Failed to find an allowable direction with negative curvature')
%                         disp('Second order conditions are satisfied.')
%                     end
%                 end
%                 
%             end

        end
        
    end
    
end


%% Other Module
%%%%%%%%%%%%%%%

if(module.other.tog)

    if(module.other.four)
        disp('******************************');
        disp('********** Module 4 **********');
        disp('******************************');
        disp(' ');
        
        [nzr, nzc] = find(J ~= 0);
        nz = find(J ~= 0);
        
        mm = [full(max(max(abs(J(nz))))), full(min(min(abs(J(nz)))))];
        logJ = log10(abs(J(nz)));
        txt1 = {'large','small'};
        txt2 = {'Largest','Smallest'};
        
        %    disp('Creating histogram of Jacobian elements');
        %    hist(logJ);
        %    xlabel('Log10 of abs(Jacobian)');
        %    ylabel('Frequency');
        %    title('Degeneracy Hunter');
        
        for i = 1:length(mm)
            disp(['Analyzing ',txt1{i},' values of the Jacobian']);
            disp([txt2{i},' value: ',num2str(mm(i),'%10.3e')]);
            
            disp(['Elements of the Jacobian within ',num2str(module.other.fourTol(i)),' order(s) of magnitude:']);
            switch i
                % Rows - equations
                % Cols - variables
                case 1
                    ii = find(logJ > log10(mm(i)) - module.other.fourTol(i));
                case 2
                    ii = find(logJ < log10(abs(mm(i))) + module.other.fourTol(i));
            end
            logJs = zeros(length(ii),1);
            for j = 1:length(ii)
                logJs(j) = logJ(ii(j));
            end
            [tmp, isrt] = sort(logJs);
            if(i == 1)
                isrt = flipud(isrt);
            end
            for k = 1:length(ii)
                j = isrt(k);
                fprintf('%s \t \t %s \t \t %e\n',eqn{nzr(ii(j))},var{nzc(ii(j))},full(J(nz(ii(j)))));
            end
            disp(' ');
        end
    end
if(module.five)
    disp('******************************');
    disp('********** Module 5 **********');
    disp('******************************');
    disp(' ');
        
    disp('Computing Singular Value Decomposition (SVD) of Jacobian...');
    drawnow('update');
    
    if(module.sparse)
        [lrgU,lrgS,lrgV] = svds(jacSmall,10);
        options.tol = 1E-12;
        options.maxit = 5000;
        options.disp = 0;
        [smlU, smlS, smlV] = svds(jacSmall,10,1E-5,options);
        sTmp = diag(smlV);
        i = sTmp > 1E-10;
        U = [lrgU, smlU(:,i)];
        V = [lrgV, smlV(:,i)];
        s = [diag(lrgS); sTmp(i)];
        for j = 1:length(i)
            if(~i(j))
                disp(['Rejecting singular value ',num2str(sTmp(j),'%10.3e'),' as numerical noise']);
            end
        end
    else
        tic
        [U,S,V] = svd(full(jacSmall));
        tocM;
        s = diag(S);
    end
    
    disp(' ');
    disp(['Condition number: ',num2str(max(s)/min(s),'%10.3e')]);
    disp('Creating histogram of singular values...');
    
    logS = log10(s);
    
    hist(logS);
    xlabel('Log10 of Singular Values')
    ylabel('Frequency')
    title('Degeneracy Hunter');
    
    flag = true;
    
    %i = length(s);
    
    sm = logS < logS(end) + 0.1;
    lrg = logS > logS(1) - 0.1;
    
    ii = find(sm | lrg);
    
    for m = 1:length(ii);
        
        i = ii(m);
        
        sml = logical(sm(i));
        disp(' ');
        disp(' ***** ');
        disp(['Analyzing SV ',num2str(i),'   value: ',num2str(s(i),'%10.3e')]);

        rsv = U(:,i);
        srsv = rsv.^2;
        ssrsv = sum(srsv);
        
%         if(sml)
%             tmp2 = find(abs(srsv) > 0);
%             [tmp, isrt] = sort(srsv(tmp2),'ascend');
%             isrt = tmp2(isrt);
%         else
        [tmp, isrt] = sort(srsv,'descend');
%         end
        
        disp('Equations most prominent in the right-vector for the SV:');
        if(sml)
            if(module.fiveVerbose)
                fprintf('Cm. Ttl. Proj. \t KKT Mult. \t Run. Sml SV \t Equation \n');  
            else
                fprintf('Cm. Ttl. Proj. \t Sm. Var. \t Sm Elmt. \t Equation \n');
            end
        else
            if(module.fiveVerbose)
                fprintf('Cm. Ttl. Proj. \t KKT Mult. \t Run. Lrg SV \t Equation\n');  
            else
                fprintf('Cm. Ttl. Proj. \t Lrg. Var. \t Lrg Elmt. \t Equation\n');
            end
        end        
        cumTotProj = 0;
        cumSrsv = 0;
        j = 0;
        flag = true;
        while flag
            j = j + 1;
            cumSrsv = srsv(isrt(j)) + cumSrsv;
            cumTotProj = sqrt(cumSrsv/ssrsv);
            
            if(module.fiveVerbose)
                if(module.sparse)
                    if(sml)
                        [Ut, St, Vt] = svds(jacSmall(isrt(1:j),:),1,1E-6);
                        mSt = St;
                    else
                        [Ut, St, Vt] = svds(jacSmall(isrt(1:j),:),1);
                        mSt = St;
                    end
                else
                    [Ut, St, Vt] = svd(full(jacSmall(isrt(1:j),:)));
                    if(sml)
                        mSt = min(diag(St(1:j,1:j)));
                    else
                        mSt = max(diag(St(1:j,1:j)));
                    end
                end
                fprintf('%10.3e \t %10.3e \t %10.3e \t %s\n',cumTotProj,marg(jacSmallEqn(isrt(j))),mSt, eqn{jacSmallEqn(isrt(j))});
                if(j == length(isrt) || mSt < 10*s(i))
                flag = false;
                end
            else
                % Lookup smallest non-zero element for this equation (row) in
                % the Jacobian
                kk = find(jacSmall(isrt(j),:) ~= 0);
                if(sml)
                    k = find(min(abs(jacSmall(isrt(j),kk))) == abs(jacSmall(isrt(j),kk)),1);
                else
                    k = find(max(abs(jacSmall(isrt(j),kk))) == abs(jacSmall(isrt(j),kk)),1);
                endV 1   value: 1.351e+03
                fprintf('%10.3e \t \t %s \t %10.3e \t \t %s\n',cumTotProj,var{kk(k)},full(jacSmall(isrt(j),kk(k))),eqn{jacSmallEqn(isrt(j))});
                end
                if(j == length(isrt) || cumTotProj > 1 - 1E-6 || j > 25)
                    flag = false;
                end
            end
        end
%         [Ut, St, Vt] = svd(full(jacSmall(isrt(1:j),:)));
%         disp(['Smallest SV with subset of printed equations: ',num2str(min(diag(St(1:j,1:j))))]);
%         disp(' ');
    end
    
    if(module.fiveCheckSmallSV)
        disp('Checking smallest singular value...');
        opt = optimset('GradObj','on');
        test = @(x) svOpt(x, jacSmall');
%        x = fminunc(test, 0.001*ones(size(jacSmall,2),1));
        x = fminunc(test, U(:,end));
        if(module.sparse)
            disp('Sparse mode active.');
        else
            disp('Sparse mode inactive.');
        end
        disp(['Smallest SV from svd/svds: ',num2str(min(s))]);
        disp(['Smallest SV from opt. alg: ',num2str(sqrt(test(x)))]);
        disp(['Diagnostic - norm(x): ',num2str(norm(x))]);
        
        rsv = x;
        srsv = rsv.^2;
        ssrsv = sum(srsv);
        
        [tmp, isrt] = sort(srsv,'descend');
        
        cumTotProj = 0;
        cumSrsv = 0;
        j = 0;
        flag = true;
        fprintf('Equation \t Cum. Total. Proj. \t KKT Mult. \n');
        while flag
            j = j + 1;
            cumSrsv = srsv(isrt(j)) + cumSrsv;
            cumTotProj = sqrt(cumSrsv/ssrsv);
            
                fprintf('%s \t %10.3e \t \t %10.3e \n',eqn{jacSmallEqn(isrt(j))},cumTotProj,marg(jacSmallEqn(isrt(j))));
            if(cumTotProj > 0.9999 || j > 15)
                flag = false;
            end
        end
    end
    
%    disp(['SV ',num2str(i),' is the next smallest with a value of ',num2str(s(i),'%10.3e')]);
end

%% Module 6
if(module.six)
    disp('******************************');
    disp('********** Module 6 **********');
    disp('******************************');
    disp(' ');
    
    for i = 1:length(module.sixSearch)
        type = -1;
        switch module.sixType{i}
            case 'v'
                ii = find(strncmpi(module.sixSearch{i},var,length(module.sixSearch{i})));
                type = 1;
                nms = var(ii);
                typeT = 'variable';
            case 'e'
                ii = find(strncmpi(module.sixSearch{i},eqn,length(module.sixSearch{i})));
                type = 2;
                nms = eqn(ii);
                typeT = 'equation';
            otherwise
                disp('Hint: Specify either "v" or "e" as type');
                disp(['Search term "',module.sixSearch{i},'" ignored']);
        end
        
        if(type > 0)
            disp('**********');
            disp(['For search parameters "',module.sixSearch{i},'" ...']);
            for j = 1:length(ii)
            	disp(['Identified ',typeT,' ',nms{j},' with the following non-zero elements']);
                
                % Recall: rows are equations, columns are variables
                
                switch type
                    case 1 % Search for variables
                        jj = find(jac(:,ii(j)) ~= 0);
                        for k = 1:length(jj)
                            fprintf('%s\t\t%f\n',eqn{jj(k)}, full(jac(jj(k),ii(j))));
                        end
                    case 2 % Search for equations
                        jj = find(jac(ii(j),:) ~= 0);
                        for k = 1:length(jj)
                            fprintf('%s\t\t%f\n',var{jj(k)}, full(jac(ii(j),jj(k))));
                        end
                end
            end
        end
    end
end
    
end

if(nargin == 3)    
    if(~isempty(varargin{3}))
        fclose(fid);
    end
end

end

function [A,jb] = frref(A,tol,type)
%FRREF   Fast reduced row echelon form.
%   R = FRREF(A) produces the reduced row echelon form of A.
%   [R,jb] = FRREF(A,TOL) uses the given tolerance in the rank tests.
%   [R,jb] = FRREF(A,TOL,TYPE) forces frref calculation using the algorithm
%   for full (TYPE='f') or sparse (TYPE='s') matrices.
%   
% 
%   Description: 
%   For full matrices, the algorithm is based on the vectorization of MATLAB's
%   RREF function. A typical speed-up range is about 2-4 times of 
%   the MATLAB's RREF function. However, the actual speed-up depends on the 
%   size of A. The speed-up is quite considerable if the number of columns in
%   A is considerably larger than the number of its rows or when A is not dense.
%
%   For sparse matrices, the algorithm ignores the tol value and uses sparse
%   QR to compute the rref form, improving the speed by a few orders of 
%   magnitude.
%
%   Authors: Armin Ataei-Esfahani (2008)
%            Ashish Myles (2012)
%
%   Revisions:
%   25-Sep-2008   Created Function
%   21-Nov-2012   Added faster algorithm for sparse matrices
%
%
% From http://www.mathworks.com/matlabcentral/fileexchange/21583-fast-reduced-row-echelon-form/content/frref.m
%
% Copyright (c) 2008, Armin Ataei, Ashish Myles
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

[m,n] = size(A);

switch nargin
  case 1,
    % Compute the default tolerance if none was provided.
    tol = max(m,n)*eps(class(A))*norm(A,'inf');
    if issparse(A)
      type = 's';
    else
      type = 'f';
    end
  case 2,
  if issparse(A)
    type = 's';
  else
    type = 'f';
  end
  case 3,
    if ~ischar(type)
      error('Unknown matrix TYPE! Use ''f'' for full and ''s'' for sparse matrices.')
    end
    type = lower(type);
    if ~strcmp(type,'f') && ~strcmp(type,'s')
      error('Unknown matrix TYPE! Use ''f'' for full and ''s'' for sparse matrices.')
    end
end


do_full = ~issparse(A) || strcmp(type,'f');

if do_full
    % Loop over the entire matrix.
    i = 1;
    j = 1;
    jb = [];
    % t1 = clock;
    while (i <= m) && (j <= n)
       % Find value and index of largest element in the remainder of column j.
       [p,k] = max(abs(A(i:m,j))); k = k+i-1;
       if (p <= tol)
          % The column is negligible, zero it out.
          A(i:m,j) = 0; %(faster for sparse) %zeros(m-i+1,1);
          j = j + 1;
       else
          % Remember column index
          jb = [jb j];
          % Swap i-th and k-th rows.
          A([i k],j:n) = A([k i],j:n);
          % Divide the pivot row by the pivot element.
          Ai = A(i,j:n)/A(i,j);    
          % Subtract multiples of the pivot row from all the other rows.
          A(:,j:n) = A(:,j:n) - A(:,j)*Ai;
          A(i,j:n) = Ai;
          i = i + 1;
          j = j + 1;
       end
    end
else
    % Non-pivoted Q-less QR decomposition computed by Matlab actually
    % produces the right structure (similar to rref) to identify independent
    % columns.
    R = qr(A);

    % i_dep = pivot columns = dependent variables
    %       = left-most non-zero column (if any) in each row
    % indep_rows (binary vector) = non-zero rows of R
    [indep_rows, i_dep] = max(R ~= 0, [], 2);
    indep_rows = full(indep_rows); % probably more efficient
    i_dep = i_dep(indep_rows);
    i_indep = setdiff(1:n, i_dep);

    % solve R(indep_rows, i_dep) x = R(indep_rows, i_indep)
    %   to eliminate all the i_dep columns
    %   (i.e. we want F(indep_rows, i_dep) = Identity)
    F = sparse([],[],[], m, n);
    F(indep_rows, i_indep) = R(indep_rows, i_dep) \ R(indep_rows, i_indep);
    F(indep_rows, i_dep) = speye(length(i_dep));

    % result
    A = F;
    jb = i_dep;
end

end
function H = nulls(A, thresh)
% NULLS   Null space of a sparse matrix.
%    Z = NULLS(A) is a basis for the null space of A.  That is, A*Z has
%    negligible elements and size(Z,2) is the nullity of A.  The algorithm
%    used [1] is optimized for sparse matrices.
%
%    Z = NULLS(A,thresh) controls pivoting between maximum preservation of
%    sparsity (thresh=0) and maximum accuracy (thresh=1). The default is
%    thresh=0.1.
%
%    [1] M. Khorramizadeh and N. Mahdavi-Amiri, "An efficient algorithm for
%    sparse null space basis problem using ABS methods," Numerical
%    Algorithms, vol. 62, no. 3, pp. 469–485, Jun. 2012.
%
%    Example:
%    A = sprand(1000, 1500, 0.001);
%    Z = nulls(A);
%    size(Z,2)    % at least 500, typically more, depending on A
%    norm(A*Z,1)  % should give a small value

%   Copyright 2013 Martin Holters

% From http://www.mathworks.com/matlabcentral/fileexchange/42922-null-space-for-sparse-matrix/content/nulls.m
%
% Copyright (c) 2013, Martin Holters
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
if nargin==1
    thresh = 0.1;
end

[m,n] = size(A);

% get ordering of the rows of A by ascending count of non-zero elements
r = sum(spones(A),2);
[~, t] = sort(r, 'ascend');

% initialize H1
H = speye(n);

% A is transposed once beforehand as column accesses are cheaper than row accesses
At = A.';
for i=1:m
    % note that H is transposed compared to [1], and hence, so is s
    s = At(:,t(i)).' * H;
    if nnz(s) == 0
        continue
    end
    % only consider non-zero entries in s
    jnz = find(s);
    % filter based on the pivoting threshold
    j = jnz(abs(s(jnz)) >= thresh*max(abs(s)));
    % from the permissible columns in H, find the one with the lowest number of non-zero elements
    [~, jj] = min(sum(spones(H(:,j))));
    j = j(jj);
    % multiplication with G from [1] can be represented by the following matrix manipulation
    H = [H(:,1:j-1), H(:,j+1:end)] - H(:,j)*[s(1:j-1), s(j+1:end)]/s(j);
end

end

function [ r ] = rankHelper( A, sparse, tol )

if(sparse)
    r = sprank(A);
else
    if(isempty(tol))
        r = rank(full(A));
    else
        r = rank(full(A),tol);
    end
end

end

function joinedStr = strjoin(c, aDelim)
%STRJOIN  Join cell array of strings into single string
%   S = STRJOIN(C) constructs the string S by linking each string within
%   cell array of strings C together with a space.
%
%   S = STRJOIN(C, DELIMITER) constructs S by linking each element of C
%   with the elements of DELIMITER. DELIMITER can be either a string or a
%   cell array of strings having one fewer element than C.
%
%   If DELIMITER is a string, then STRJOIN forms S by inserting DELIMITER
%   between each element of C. DELIMITER can include any of these escape
%   sequences:
%       \\   Backslash
%       \0   Null
%       \a   Alarm
%       \b   Backspace
%       \f   Form feed
%       \n   New line
%       \r   Carriage return
%       \t   Horizontal tab
%       \v   Vertical tab
%
%   If DELIMITER is a cell array of strings, then STRJOIN forms S by
%   interleaving the elements of DELIMITER and C. In this case, all
%   characters in DELIMITER are inserted as literal text, and escape
%   characters are not supported.
%
%   Examples:
%
%       c = {'one', 'two', 'three'};
%
%       % Join with space.
%       strjoin(c)
%       % 'one two three'
%
%       % Join as a comma separated list.
%       strjoin(c, ', ')
%       % 'one, two, three'
%
%       % Join with a cell array of strings DELIMITER.
%       strjoin(c, {' + ', ' = '})
%       % 'one + two = three'
%
%   See also STRCAT, STRSPLIT.

%   Copyright 2012 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/11/15 13:56:49 $

narginchk(1, 2);

% Check input arguments.
if ~isCellString(c)
    error(message('MATLAB:strjoin:InvalidCellType'));
end
if nargin < 2
    aDelim = ' ';
end

% Allocate a cell to join into - the first row will be C and the second, D.
numStrs = numel(c);
joinedCell = cell(2, numStrs);
joinedCell(1, :) = reshape(c, 1, numStrs);
if isString(aDelim)
    if numStrs < 1
        joinedStr = '';
        return;
    end
    escapedDelim = strescape(aDelim);
    joinedCell(2, 1:numStrs-1) = {escapedDelim};
elseif isCellString(aDelim)
    numDelims = numel(aDelim);
    if numDelims ~= numStrs - 1
        error(message('MATLAB:strjoin:WrongNumberOfDelimiterElements'));
    end
    joinedCell(2, 1:numDelims) = aDelim(:);
else
    error(message('MATLAB:strjoin:InvalidDelimiterType'));
end

% Join.
joinedStr  = [joinedCell{:}];

end