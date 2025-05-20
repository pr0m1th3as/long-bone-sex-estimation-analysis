## Load data sets
load Femur-Data.mat
load Humerus-Data.mat
load Tibia-Data.mat
load Ulna-Data.mat

## Fit classification models for each variable of each bone and construct a
## table with univariate classification rankings based on prediction outcomes
Fvars = {};
order = [];
for i = 1:size (Xf, 2)
  tmp = Zf(:,i);
  mdl = fitcdiscr (tmp, Yf, "Gamma", 0.3, "Prior", "uniform");
  ## Remove samples with missing values
  TF = ! logical (sum (isnan (tmp), 2));
  ns = sum (TF);
  out = predict (mdl, tmp(TF,:));
  out = sum (out == Yf(TF)) / ns;
  Fvars = [Fvars; FemurData(i+1,1), {out}, {ns}, {F_mu(i)}, {Fstd(i)}, FemurData(i+1,4:7)];
  order = [order; out];
endfor
[~, Fvidx] = sort (order, "descend");
Fvars = Fvars(Fvidx,:);

Hvars = {};
order = [];
for i = 1:size (Xh, 2)
  tmp = Zh(:,i);
  mdl = fitcdiscr (tmp, Yh, "Gamma", 0.3, "Prior", "uniform");
  ## Remove samples with missing values
  TF = ! logical (sum (isnan (tmp), 2));
  ns = sum (TF);
  out = predict (mdl, tmp(TF,:));
  out = sum (out == Yh(TF)) / ns;
  Hvars = [Hvars; HumerusData(i+1,1), {out}, {ns}, {H_mu(i)}, {Hstd(i)}, HumerusData(i+1,4:7)];
  order = [order; out];
endfor
[~, Hvidx] = sort (order, "descend");
Hvars = Hvars(Hvidx,:);

Tvars = {};
order = [];
for i = 1:size (Xt, 2)
  tmp = Zt(:,i);
  mdl = fitcdiscr (tmp, Yt, "Gamma", 0.3, "Prior", "uniform");
  ## Remove samples with missing values
  TF = ! logical (sum (isnan (tmp), 2));
  ns = sum (TF);
  out = predict (mdl, tmp(TF,:));
  out = sum (out == Yt(TF)) / ns;
  Tvars = [Tvars; TibiaData(i+1,1), {out}, {ns}, {T_mu(i)}, {Tstd(i)}, TibiaData(i+1,4:7)];
  order = [order; out];
endfor
[~, Tvidx] = sort (order, "descend");
Tvars = Tvars(Tvidx,:);

Uvars = {};
order = [];
for i = 1:size (Xu, 2)
  tmp = Zu(:,i);
  mdl = fitcdiscr (tmp, Yu, "Gamma", 0.1, "Prior", "uniform");
  ## Remove samples with missing values
  TF = ! logical (sum (isnan (tmp), 2));
  ns = sum (TF);
  out = predict (mdl, tmp(TF,:));
  out = sum (out == Yu(TF)) / ns;
  Uvars = [Uvars; UlnaData(i+1,1), {out}, {ns}, {U_mu(i)}, {Ustd(i)}, UlnaData(i+1,4:7)];
  order = [order; out];
endfor
[~, Uvidx] = sort (order, "descend");
Uvars = Uvars(Uvidx,:);

## Pad cell arrays with extra rows so that they can be concatenated
maxrows = max (cellfun (@(x) size (x, 1), {Fvars, Hvars, Tvars, Uvars}));
padrows = cellfun (@(x) maxrows - size (x, 1), {Fvars, Hvars, Tvars, Uvars});
Fvars = [Fvars; repmat({"", NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}, padrows(1), 1)];
Hvars = [Hvars; repmat({"", NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}, padrows(2), 1)];
Tvars = [Tvars; repmat({"", NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}, padrows(3), 1)];
Uvars = [Uvars; repmat({"", NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN}, padrows(4), 1)];
Table = [Fvars, Hvars, Tvars, Uvars];

varNames = {'Measurement', 'Accuracy', 'Samples', 'Mean', 'StD', 'M_lo', 'M_hi', 'F_lo', 'F_hi'};
Femur = cell2table (Fvars, 'VariableNames', varNames);
Humerus = cell2table (Hvars, 'VariableNames', varNames);
Tibia = cell2table (Tvars, 'VariableNames', varNames);
Ulna = cell2table (Uvars, 'VariableNames', varNames);
T = table (Femur, Humerus, Tibia, Ulna);

## Save to CSV file
cell2csv ("univar_stats.csv", [varNames, varNames, varNames, varNames; Table]);

## Keep 31 best performing variables (accuracy: >70%)
Vidx = [Fvidx(1:31), Hvidx(1:31), Tvidx(1:31), Uvidx(1:31)];

## Save a few new variables
save ('-binary', 'univar-index.mat', 'Fvars', 'Hvars', 'Tvars', 'Uvars', 'Vidx');

## Find measurement indices corresponding to DATA returned by longbone_Geometry
Header = {'Max Distance', ...
        'Area 20%', 'Perimeter 20%', 'ArPerIndex 20%', 'Ix 20%', 'Iy 20%', ...
        'Ixy 20%', 'Ix/Iy 20%', 'Imin 20%', 'Imax 20%', 'Imax/Imin 20%', ...
        'theta 20%', 'Dihedral angle 20_35', ......
        'Area 35%', 'Perimeter 35%', 'ArPerIndex 35%', 'Ix 35%', 'Iy 35%', ...
        'Ixy 35%', 'Ix/Iy 35%', 'Imin 35%', 'Imax 35%', 'Imax/Imin 35%', ...
        'theta 35%', 'Dihedral angle 35_50', ...
        'Area 50%', 'Perimeter 50%', 'ArPerIndex 50%', 'Ix 50%', 'Iy 50%', ...
        'Ixy 50%', 'Ix/Iy 50%', 'Imin 50%', 'Imax 50%', 'Imax/Imin 50%', ...
        'theta 50%', 'Dihedral angle 50_65', ...
        'Area 65%', 'Perimeter 65%', 'ArPerIndex 65%', 'Ix 65%', 'Iy 65%', ...
        'Ixy 65%', 'Ix/Iy 65%', 'Imin 65%', 'Imax 65%', 'Imax/Imin 65%', ...
        'theta 65%', 'Dihedral angle 65_80', ...
        'Area 80%', 'Perimeter 80%', 'ArPerIndex 80%', 'Ix 80%', 'Iy 80%', ...
        'Ixy 80%', 'Ix/Iy 80%', 'Imin 80%', 'Imax 80%', 'Imax/Imin 80%', ...
        'theta 80%', 'Diaphyseal Bending'};
f = @(x) find (strcmpi (Header, x));
Fvars(:,1) = cellfun (f, Fvars(:,1), 'UniformOutput', false);
Hvars(:,1) = cellfun (f, Hvars(:,1), 'UniformOutput', false);
Tvars(:,1) = cellfun (f, Tvars(:,1), 'UniformOutput', false);
Uvars(:,1) = cellfun (f, Uvars(:,1), 'UniformOutput', false);

## Remove accuracy and sample size columns and keep 31 best performing variables
## and convert to numeric matrices
Fvars = cell2mat (Fvars(1:31,[1,4:9]));
Hvars = cell2mat (Hvars(1:31,[1,4:9]));
Tvars = cell2mat (Tvars(1:31,[1,4:9]));
Uvars = cell2mat (Uvars(1:31,[1,4:9]));

## Save to descriptives.mat file
filename = fullfile (pwd, 'private', 'descriptives.mat');
save ('-binary', filename, 'Fvars', 'Hvars', 'Tvars', 'Uvars');
