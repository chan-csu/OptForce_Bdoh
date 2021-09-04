function [v,v_fit,stat,sol] = fitExpFlux(CbModel,rxns_exp,v_exp,weight,minNorm)
%[v,v_fit,stat,sol] = FitExpFlux(CbModel,rxns_exp,v_exp,weight,minNorm)
%Given a set of experimentally measured flux v_exp of reactions rxns_exp,
%fitExpFlux find the a sum-of-squared-error solution flux distribution v by
%solving the following quadratic program (QP) using SolveCobraQP in the COBRA toolbox:
%   min sum(w_j * (v_j - v_exp,j)^2)
%   subject to Sv = 0
%   lb <= v <= ub
%
%Input: 
%CbModel: must contain the field 'S', 'lb' and 'ub'
%rxns_exp: the set of reactions with fluxes experimentally measured.
%   Can either be:
%   a vector (n x 1) of numbers indicating the indices of the reactions 
%   in the column of CbModel.S; or
%   a cell (n x 1) containing the name of reactions in CbModel.rxns 
%   (must contain the field 'rxns' in this case), e.g. {'Ex_glc(e)','Ex_ac(e)'}
%v_exp: experimentally measured flux values (n x 1 vector)
%
%Optional input:
%weight: the weight for the minimization problem. 
%   default all equal to 1.
%   e.g. For the best linear unbiased estimator (BLUE), one can use the
%   reciprocal of the variance of measurement for each flux as weight
%   (1/sd^2)
% minNorm       (COBRA parameter, if not provided, get from COBRA
%               parameters, if not exist, default to be 1)
%               (modified from optimizeCbModel.m) 
%               {(0), 'one', > 0 , n x 1 vector}, where [m,n]=size(S);
%                0      normal LP
%                'one'  Minimise the Taxicab Norm using LP.
%                                 min |v|
%                                   s.t. S*v = b
%                                        c'v = f
%                                        lb <= v <= ub
%                -----
%                The remaining options work only with a valid QP solver:
%                -----
%                > 0    (default situation)
%                       Minimises the Euclidean Norm of internal fluxes.
%                       Typically 1e-6 works well.
%                                 min ||v||
%                                   s.t. S*v = b
%                                        c'v = f
%                                        lb <= v <= ub
%               n x 1   Forms the diagonal of positive definiate
%                       matrix F in the quadratic program
%                               min 0.5*v'*F*v
%                               st. S*v = b
%                                   c'*v = f
%                                   lb <= v <= ub
%
%Output:
%v: the whole fitted flux distribution 
%v_fit: the fitted fluxes for the reacitons with provided flux
%state is a structure reporting the following statistics:
%   res: residue for each fluxes
%   sse: the weighted sum of squared error (the error minimized by this function)
%sol: original solution returned by SolveCobraQP
%
%
%Siu Hung Joshua Chan       26/11/2014


%% Consistency checking for function inputs
%Check the weight vector
if ~exist('weight','var')
    weight = ones(numel(v_exp),1);
elseif numel(weight) ~= numel(v_exp)
    warning('Weight vector not in the same size as v_exp, not used')
    weight = ones(numel(v_exp),1);
elseif any(weight<0)
    warning('Weight vector contains negative entries, not used')
    weight = ones(numel(v_exp),1);
end
weight = weight(:);
% Check minNorm 
if ~exist('minNorm','var') 
    minNorm = 1;
end
%Convert string of rxns_exp into cell array if needed
v_exp = v_exp(:);
if ischar(rxns_exp)
    rxns_exp = cellstr(rxns_exp);
end
rxns_exp = rxns_exp(:);
if numel(rxns_exp) ~= numel(v_exp)
    error('rxns_exp and v_exp not in the same size')
end
%Check the 'rxns_exp', 'v_exp' and 'CbModel' compatibility
n = size(CbModel.S,2);
if numel(CbModel.lb) ~= n || numel(CbModel.ub)~= n
    error('Size of lb or ub not the same as the number of columns of S in CbModel.')
end
if isfield(CbModel,'rxns')
    if numel(CbModel.rxns) ~= n
        error('Size of rxns in CbModel not correct.')
    end
end
m = size(CbModel.S,1);
if isfield(CbModel,'b')
    if numel(CbModel.b) ~= m
        error('Size of b and the number of rows of S in CbModel not the same.')
    end
end

%% Pre-processing
%   Convert reaction names in rxns_exp into indices rxns_ind
if iscell(rxns_exp)
    rxns_ind = zeros(size(rxns_exp));
    if ~isfield(CbModel,'rxns')
            error('rxns is not a field in CbModel')
    else
        for j = 1:numel(rxns_exp);
            f = find(strcmp(rxns_exp{j},CbModel.rxns),1);
            if ~isempty(f)
                rxns_ind(j) = f;
            end
        end
    end
    if ~any(rxns_ind)
        error('rxns_exp does not have any match in CbModel.rxns')
    elseif ~all(rxns_ind)
        warning(['The following reaction fluxes are ignored due to no '...
            'any match in CbModel.rxn: ' strjoin(rxns_exp(rxns_ind==0)')])
        v_exp = v_exp(rxns_ind~=0);
        rxns_ind = rxns_ind(rxns_ind~=0);
        weight = weight(rxns_ind~=0);
    end
else
    rxns_ind = rxns_exp;
    OutOfRange = rxns_ind > size(CbModel.S,2);
    if any(OutOfRange)
        warning(['The following reaction fluxes are ignored because they '...
            'are out of the range of CbModel.S: ' num2str(rxns_ind(OutOfRange)')])
    end
    v_exp = v_exp(~OutOfRange);
    rxns_ind = rxns_ind(~OutOfRange);
    weight = weight(~OutOfRange);
end


%% Problem formulation for minimizing the error
QPproblem = struct();
QPproblem.A = sparse(CbModel.S); %assign constraint matrix
if isfield(CbModel,'b')
    b = CbModel.b;
else
    b = zeros(size(CbModel.S,1),1);
end
QPproblem.b = b;
%assign the quadratic objective matrix
QPproblem.F = sparse(n,n);
for j = 1:numel(rxns_ind)
    QPproblem.F(rxns_ind(j),rxns_ind(j)) = 2*weight(j);
end
%assign the linear objective vector
QPproblem.c = zeros(n,1);
QPproblem.c(rxns_ind) = -2*weight.*v_exp;
%assign bounds
QPproblem.lb = CbModel.lb;
QPproblem.ub = CbModel.ub;
%assign sense
QPproblem.osense = 1;
QPproblem.csense = char('E'*ones(1,m));
%solve QP 
sol = solveCobraQP(QPproblem);
%flux vector
v = sol.full;

%% Minimizing the norm if needed (modified from optimizeCbModel.m)
%Cannot implement the objective funtion value as constraint because
%quadratic constraint is not allowed in COBRA.
%Instead, fix the fitted fluxes for reactions with provided values, and
%minimize.
if strcmp(minNorm, 'one')
    rev = CbModel.lb < 0 & CbModel.ub > 0; %true reversible reactions
    nRE = sum(rev); %number of reversible reactions
    revRE = zeros(n,1); 
    revRE(rev) = n+(1:nRE)'; %indices for reverse direction of reversible reactions
    LPproblem = struct();
    LPproblem.A = [CbModel.S -CbModel.S(:,rev)];
    LPproblem.c = ones(n,1);
    LPproblem.c(CbModel.ub < 0) = -1; %for reactions always in reverse directon
    LPproblem.c = [LPproblem.c; ones(nRE,1)];
    LPproblem.lb = CbModel.lb;
    LPproblem.lb(rev) = 0; %update lower bound
    LPproblem.lb = [LPproblem.lb; zeros(nRE,1)];
    LPproblem.ub = [CbModel.ub; -CbModel.lb(rev)]; %upper bound
    rxns_indL = false(n,1);
    rxns_indL(rxns_ind) = true;
    rxnsIR_indL = rxns_indL & ~rev;
    LPproblem.lb(rxnsIR_indL) = v(rxnsIR_indL) - 1e-12; %satisfied the fitted reaction fluxes
    LPproblem.ub(rxnsIR_indL) = v(rxnsIR_indL) + 1e-12;
    rxnsRE_indL = rxns_indL & rev;
    if any(rxnsRE_indL)
        rxnsREpos_indL = rxnsRE_indL & v >= 0;
        if any(rxnsREpos_indL)
            LPproblem.lb(rxnsREpos_indL) = v(rxnsREpos_indL) - 1e-12;
            LPproblem.ub(rxnsREpos_indL) = v(rxnsREpos_indL) + 1e-12;
            LPproblem.ub(revRE(rxnsREpos_indL)) = 0;
        end
        rxnsREneg_indL = rxnsRE_indL & v < 0;
        if any(rxnsREneg_indL)
            LPproblem.ub(rxnsREneg_indL) = 0;
            LPproblem.lb(revRE(rxnsREneg_indL)) = abs(v(rxnsREneg_indL)) - 1e-12;
            LPproblem.ub(revRE(rxnsREneg_indL)) = abs(v(rxnsREneg_indL)) + 1e-12;
        end
    end
    LPproblem.b = b;
    LPproblem.osense = 1;
    LPproblem.csense = char('E'*ones(1,m));
    sol = solveCobraLP(LPproblem);
    v0 = sol.full;
    v = v0(1:n);
    v(rev) = v(rev) - v0(n+1:end);
elseif length(minNorm)> 1 || minNorm > 0
    % quadratic minimization of the norm.
    % set previous optimum as constraint.
    if length(minNorm)==1
        minNorm=ones(n,1)*minNorm;
    end
    QPproblem.F = spdiags(minNorm,0,n,n);
    QPproblem.c = zeros(n,1);
    QPproblem.lb(rxns_ind) = v(rxns_ind) - 1e-12;
    QPproblem.ub(rxns_ind) = v(rxns_ind) + 1e-12;
    sol = solveCobraQP(QPproblem);
    %flux vector
    v = sol.full;
end
%% Function output
v_fit = v(rxns_ind);
%statistics:
res = v_exp - v(rxns_ind);
sseW = res' * (weight .* res);
stat=struct('res',res,'sse',sseW);

end