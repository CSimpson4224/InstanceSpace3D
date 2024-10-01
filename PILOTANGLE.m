function ViewAngle = PILOTANGLE(X, Y, opts)
% -------------------------------------------------------------------------
% PILOT.m
% -------------------------------------------------------------------------
%
% By: Mario Andres Munoz Acosta
%     School of Mathematics and Statistics
%     The University of Melbourne
%     Australia
%     2020
%
% -------------------------------------------------------------------------

errorfcn = @(alpha,Xbar,n,n2) nanmean(nanmean((Xbar(:,(n+1):end)-(reshape(alpha((2*n)+1:end),n2,2)*... % C
                                                    reshape(alpha(1:2*n),2,n)...        % A
                                                   *Xbar(:,1:n)')').^2,1),2)...
                                                   +0.2*abs(dot(alpha(1+2*(0:(n-1))),alpha((1:n)*2)))...
                                                   + (0-1);                                               
                                               
                                               
n = size(X, 2); % Number of features
n2 = size(Y, 2);
Xbar = [X Y];
m = size(Xbar, 2);
Hd = pdist(X)';
if exist('gcp','file')==2
    mypool = gcp('nocreate');
    if ~isempty(mypool)
        nworkers = mypool.NumWorkers;
    else
        nworkers = 0;
    end
else
    nworkers = 0;
end



if isfield(opts,'X0') && isnumeric(opts.X0) && ...
        size(opts.X0,1)==2*m+2*n && size(opts.X0,2)>=1
    disp('  -> PILOT is using a user defined starting points for BFGS.');
    X0 = opts.X0;
    opts.ntries = size(opts.X0,2);
else
    disp('  -> PILOT is using a random starting points for BFGS.');
    state = rng;
    rng('default');
    X0 = 2*rand(2*n2+2*n, opts.ntries)-1;
    rng(state);
end
alpha = zeros(2*n2+2*n, opts.ntries);
eoptim = zeros(1, opts.ntries);
perf = zeros(1, opts.ntries);
disp('-------------------------------------------------------------------------');
disp('  -> PILOT is solving numerically the projection problem.');
disp('  -> This may take a while. Trials will not be run sequentially.');
disp('-------------------------------------------------------------------------');
parfor (i=1:opts.ntries,nworkers)
    [alpha(:,i),eoptim(i)] = fminunc(errorfcn, X0(:,i), ...
                                     optimoptions('fminunc','Algorithm','quasi-newton',...
                                                            'Display','off',...
                                                            'UseParallel',false,...
                                                            'MaxIterations',30000,...
                                                            'FunctionTolerance',1e-20),...
                                     Xbar, n, n2);
    aux = alpha(:,i);
    A = reshape(aux(1:2*n),2,n);
    smA = sum(A.^2');
    A(1,:) = A(1,:)/sqrt(smA(1))
    A(2,:) = A(2,:)/sqrt(smA(2))
    Z = X*A';
    perf(i) = corr(Hd,pdist(Z)');
    disp(['    -> PILOT has completed trial ' num2str(i)]);
end
out.X0 = X0;
out.alpha = alpha;
out.eoptim = eoptim;
out.perf = perf;
[~,idx] = max(out.perf);

out.A = reshape(out.alpha(1:2*n,idx),2,n);
smA = sum(out.A.^2');
out.A(1,:) = out.A(1,:)/sqrt(smA(1));
out.A(2,:) = out.A(2,:)/sqrt(smA(2));
out.Z = X*out.A';
ViewAngle = cross(out.A(1,:),out.A(2,:));


disp('-------------------------------------------------------------------------');
disp('  -> PILOT has completed. The angle vector is:');
disp(ViewAngle);

end