function [res] = myExpsForItersBttldBsctn(fname)

ind = ssget;

minmn = min(ind.nrows', ind.ncols');

%% ORG larges for run time comparison

probs = minmn == ind.sprank' & ind.nrows' >=100000 & ind.nrows' < 50000000 & ind.nnz' < 250000000 & ind.nzero'==0;


probs = minmn == ind.sprank' & ind.nrows' == ind.ncols' & ind.nrows' >=100000 & ind.nrows' < 50000000 & ind.nnz' < 250000000 & ind.nzero'==0;

%% TEST mediums for correctness.
%probs = minmn == ind.sprank' & ind.nrows' >=10000 & ind.nrows'< 100000 & ind.nnz' < 250000000 & ind.nzero'==0;

%% TEST smalls for correctness.
%probs = minmn == ind.sprank' & ind.nrows' >=1000 & ind.nrows'< 10000 & ind.nnz' < 250000000 & ind.nzero'==0;

%% TEST tiny for correctness.
%probs = minmn == ind.sprank' & ind.nrows' >=100 & ind.nrows'< 1000 & ind.nnz' < 250000000 & ind.nzero'==0;


%% TEST minuscule for correctness.
%probs = minmn == ind.sprank' & ind.nrows' >=10 & ind.nrows'< 100 & ind.nnz' < 250000000 & ind.nzero'==0;
%probs = ind.nrows' ~= ind.ncols' &  minmn<5000 & ind.nzero'==0;
probs = find(probs);

if exist(fname, 'file')
    rawData = load(fname);
    solved = rawData(:, 1);
    probs = setdiff(probs, solved);
end

res = zeros(length(probs), 7);
for ii = 1:length(probs)
    fprintf('\t\t%d: solving %d of %d --- %s\n', probs(ii), ii, length(probs), deal(ind.Name{probs(ii)}));
    prob = ssget(probs(ii));
    A = abs(prob.A);
%     rdegs = spones(A) * ones(size(A, 2), 1);
%     cdegs = spones(A') * ones(size(A, 1), 1);
%     A = A(rdegs > 0, cdegs > 0);
%     
%     [m, n] = size(A);
%     if(m < n)
%         A = A';
%         fprintf('this is rectanglar n < m\n');
%     end
%     
    [pp, tt, itBttlA, itBsctA] = pureBisectionBasedOnMC64J3_runner(A);
    [r, c, As, numIters, err] = buScaleSK( A, 20, 1.0e-10);
    [pp, tt, itBttlAs, itBsctAs] = pureBisectionBasedOnMC64J3_runner(As);
    
    [r, c, As, numIters, err] = buScaleSK(spones(A), 20, 1.0e-20);
    [pp, tt, itBttlPattAs, itBsctPattAs] = pureBisectionBasedOnMC64J3_runner(As);
    
    res(ii, 1:7) = [probs(ii), itBttlA, itBsctA, itBttlAs, itBsctAs, itBttlPattAs, itBsctPattAs];
    fop = fopen(fname, 'a');
    fprintf(fop, '%d ', probs(ii));
    fprintf(fop, '%d %d %d %d %d %d ', itBttlA, itBsctA, itBttlAs, itBsctAs, itBttlPattAs, itBsctPattAs);
    fprintf(fop,'%%%s\n', deal(prob.name));
    fclose(fop);
    
end