## Load aggregated data
load Bone-Data.mat

UlnaData = {'Variables', 'Males', 'Females', 'M_Lo', 'M_Hi', 'F_Lo', 'F_Hi'};
M_bs = []; M_bs_out = {};
F_bs = []; F_bs_out = {};
Pval_CI95 = [];
## Go through all variables
for i = 6:66
  ## Males both sides
  Msample = Ulna(Ulna(:,4) == 1, [1, 2, 3, i]);
  s = statistics (Msample(:,end));
  IQR = (s(4) - s(2)) * 1.5;
  M_outlo = s(2) - IQR;
  M_outhi = s(4) + IQR;
  idxlo = Msample(:,end) < M_outlo;
  idxhi = Msample(:,end) > M_outhi;
  M_bs_out = [M_bs_out, {[Msample(idxlo,[1:3]); Msample(idxhi,[1:3])]}];
  M_bs = [M_bs; s([6, 7])', sum([idxlo; idxhi])];
  ## Keep samples without outliers
  M = Msample(! idxlo & ! idxhi, :);

  ## Females both sides
  Fsample = Ulna(Ulna(:,4) == 2, [1, 2, 3, i]);
  s = statistics (Fsample(:,end));
  IQR = (s(4) - s(2)) * 1.5;
  F_outlo = s(2) - IQR;
  F_outhi = s(4) + IQR;
  idxlo = Fsample(:,end) < F_outlo;
  idxhi = Fsample(:,end) > F_outhi;
  F_bs_out = [F_bs_out, {[Fsample(idxlo,[1:3]); Fsample(idxhi,[1:3])]}];
  F_bs = [F_bs; s([6, 7])', sum([idxlo; idxhi])];
  ## Keep samples without outliers
  F = Fsample(! idxlo & ! idxhi, :);

  ## T-test excluding outliers
  [~, pval, ci, stats] = ttest2 (M(:,end), F(:,end));
  Pval_CI95 = [Pval_CI95; pval, ci', stats.tstat];

  ## Aggregate Ulna data sample without outliers only from variables
  ## that exhibit sexual dimorphism
  if (pval < 0.05)
    UlnaData = [UlnaData; Header(i), {M}, {F}, ...
                          {M_outlo}, {M_outhi}, {F_outlo}, {F_outhi}];
  endif

endfor
varNames = {'Male Mean & SD', 'Female Mean & SD', 'Diagnostic', ...
            'p-value', 'CI', 't-statistic', 'Sample Size'};

UlnaOutlierIDs = [{'Variables', 'Males', 'Females'}; ...
                  Header(6:66)', M_bs_out(:), F_bs_out(:)];

#UlnaStats = table (M_bs(:,[1:2]), F_bs(:,[1,2]), Pval_CI95(:,1) < 0.05, ...
#                   Pval_CI95(:,1), Pval_CI95(:,[2:3]), Pval_CI95(:,4), ...
#                   size (Ulna, 1) - sum ([M_bs(:,3), F_bs(:,3)], 2), ...
#                   'VariableNames', varNames, 'RowNames', Header(6:66))

## Create a matrix for multivariate analysis by aggregating variables per sample
## Find unique samples
Midx = [];
Fidx = [];
for i = 2:size (UlnaData, 1)
  Midx = [Midx; UlnaData{i,2}(:,[1:3])];
  Fidx = [Fidx; UlnaData{i,3}(:,[1:3])];
endfor
Midx = unique (Midx, 'rows');
Fidx = unique (Fidx, 'rows');

## Append each variable to respective samples, remaining samples get NaN
for i = 1:size(Midx, 1)
  for v = 2:size (UlnaData, 1)
    tmp = UlnaData{v,2};
    idx = Midx(i,1) == tmp(:,1) & Midx(i,2) == tmp(:,2) & Midx(i,3) == tmp(:,3);
    if (sum (idx) > 0) # found
      Midx(i,v+2) = tmp(idx,4);
    else
      Midx(i,v+2) = NaN;
    endif
  endfor
endfor
for i = 1:size(Fidx, 1)
  for v = 2:size (UlnaData, 1)
    tmp = UlnaData{v,3};
    idx = Fidx(i,1) == tmp(:,1) & Fidx(i,2) == tmp(:,2) & Fidx(i,3) == tmp(:,3);
    if (sum (idx) > 0) # found
      Fidx(i,v+2) = tmp(idx,4);
    else
      Fidx(i,v+2) = NaN;
    endif
  endfor
endfor

## Merge measurment variables of males and females into a single data matrix X
## and a corresponding group vector Y, and compute Z-scores for training/testing
Xu = [Midx(:,4:end); Fidx(:,4:end)];
Yu = [ones(size(Midx(:,1))); 2*ones(size(Fidx(:,1)))];
U_mu = mean (Xu, 'omitnan');
Ustd = std (Xu, 'omitnan');
Zu = (Xu - U_mu) ./ Ustd;

## Save variables and clean up workspace
clear -x UlnaData UlnaOutlierIDs UlnaStats Xu Yu Zu U_mu Ustd
save ('-binary', 'Ulna-Data.mat', 'UlnaData', 'UlnaOutlierIDs', ...
      'Xu', 'Yu', 'Zu', 'U_mu', 'Ustd');
