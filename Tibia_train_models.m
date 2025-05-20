## Load required variables
load Tibia-Data.mat
load ("-binary", "univar-index.mat", "Vidx");

## Keep best performing measurements and remove samples missing values
Zt = Zt(:,Vidx(:,3));
TF = ! logical (sum (isnan (Zt), 2));
Zt = Zt(TF,:);
Yt = Yt(TF);

## Print sample information
fprintf ("Training Tibia models: %d samples (%d males, %d females)\n\n", ...
         length (Yt), sum (Yt == 1), sum (Yt == 2));

## Save classification models into a cell array
models = {};

## Get scores from multivariate classification using all measurements
mdl = fitcdiscr (Zt, Yt, "Gamma", 0.1, "Prior", "uniform");
mdl = compact (mdl);
[out, Scores] = predict (mdl, Zt);
acc = sum (out == Yt) / sum (TF);
models = [models; {mdl}];
fprintf ("Multivariable Linear Discriminant Classification yields %0.2f%% accuracy.\n", acc * 100);

mdl = fitcknn (Zt, Yt, "NumNeighbors", 15);
[out, score] = predict (mdl, Zt);
acc = sum (out == Yt) / sum (TF);
models = [models; {mdl}];
fprintf ("Multivariable K-Nearest Neighbor Classification yields %0.2f%% accuracy.\n", acc * 100);

## Aggregate scores
Scores = [Scores, score];

## Get scores for each measurement
LDA = {};
KNN = {};
for i = 1:size (Vidx, 1)
  Z = Zt(:,i);
  mdl = fitcdiscr (Z, Yt, "Gamma", 0.1, "Prior", "uniform");
  mdl = compact (mdl);
  [out, score] = predict (mdl, Z);
  models = [models; {mdl}];
  LDA = [LDA, {[out, score, Yt]}];
  Scores = [Scores, score];

  mdl = fitcknn (Z, Yt, "NumNeighbors", 15);
  [out, score] = predict (mdl, Z);
  models = [models; {mdl}];
  KNN = [KNN, {[out, score, Yt]}];
  Scores = [Scores, score];
endfor

## Train a FCNN model on aggregated scores
mdl = fitcnet (Scores, Yt, 'LayerSizes', [80,40,20]);
mdl = compact (mdl);
[out, s] = predict (mdl, Scores);
models = [models; {mdl}];
fprintf ("\nFCNN Classification on Scores yields %0.2f%% accuracy.\n\n", ...
         (sum (out == Yt) / sum (TF)) * 100);

## Group scores according to their magnitude between males and females
for i = 1:size (Scores, 1)
  m_1(i,:) = any (Scores(i,[1:2:end]) == 1 & Scores(i,[2:2:end]) != 1);
  f_1(i,:) = any (Scores(i,[2:2:end]) == 1 & Scores(i,[1:2:end]) != 1);

  m09id{i} = find (Scores(i,[1:2:end]) >= 0.9 & Scores(i,[1:2:end]) < 1);
  m09(i,:) = [any(Scores(i,[1:2:end]) >= 0.9 & Scores(i,[1:2:end]) < 1), numel(m09id{i})];

  f09id{i} = find (Scores(i,[2:2:end]) >= 0.9 & Scores(i,[2:2:end]) < 1);
  f09(i,:) = [any(Scores(i,[2:2:end]) >= 0.9 & Scores(i,[2:2:end]) < 1), numel(f09id{i})];

  m08id{i} = find (Scores(i,[1:2:end]) >= 0.8 & Scores(i,[1:2:end]) < 0.9);
  m08(i,:) = [any(Scores(i,[1:2:end]) >= 0.8 & Scores(i,[1:2:end]) < 0.9), numel(m08id{i})];

  f08id{i} = find (Scores(i,[2:2:end]) >= 0.8 & Scores(i,[2:2:end]) < 0.9);
  f08(i,:) = [any(Scores(i,[2:2:end]) >= 0.8 & Scores(i,[2:2:end]) < 0.9), numel(f08id{i})];

  m07id{i} = find (Scores(i,[1:2:end]) >= 0.7 & Scores(i,[1:2:end]) < 0.8);
  m07(i,:) = [any(Scores(i,[1:2:end]) >= 0.7 & Scores(i,[1:2:end]) < 0.8), numel(m07id{i})];

  f07id{i} = find (Scores(i,[2:2:end]) >= 0.7 & Scores(i,[2:2:end]) < 0.8);
  f07(i,:) = [any(Scores(i,[2:2:end]) >= 0.7 & Scores(i,[2:2:end]) < 0.8), numel(f07id{i})];

  m06id{i} = find (Scores(i,[1:2:end]) >= 0.6 & Scores(i,[1:2:end]) < 0.7);
  m06(i,:) = [any(Scores(i,[1:2:end]) >= 0.6 & Scores(i,[1:2:end]) < 0.7), numel(m06id{i})];

  f06id{i} = find (Scores(i,[2:2:end]) >= 0.6 & Scores(i,[2:2:end]) < 0.7);
  f06(i,:) = [any(Scores(i,[2:2:end]) >= 0.6 & Scores(i,[2:2:end]) < 0.7), numel(f06id{i})];

  m05id{i} = find (Scores(i,[1:2:end]) >= 0.5 & Scores(i,[1:2:end]) < 0.6);
  m05(i,:) = [any(Scores(i,[1:2:end]) >= 0.5 & Scores(i,[1:2:end]) < 0.6), numel(m05id{i})];

  f05id{i} = find (Scores(i,[2:2:end]) >= 0.5 & Scores(i,[2:2:end]) < 0.6);
  f05(i,:) = [any(Scores(i,[2:2:end]) >= 0.5 & Scores(i,[2:2:end]) < 0.6), numel(f05id{i})];
endfor

## Aggregate grouped scores
M = [Yt, m_1, f_1, m09, f09, m08, f08, m07, f07, m06, f06, m05, f05];

## Find samples that are scored with 100% certainty
M_1 = [find(m_1); find(f_1)];

sp1 = numel (M_1) / size (M, 1) * 100;
fprintf ("Sample Percentage with 100%% certainty: %0.2f%%\n\n", sp1);

## Remove these samples before assessing next score range
M1 = M;
M1(M_1,:) = [];
M1(:,2:3) = [];
s(M_1,:) = [];

## Find samples that score at least one >=0.9 for only one sex and FCNN is in accordance
M09 = find (M1(:,2) == 1 & M1(:,4) == 0 & s(:,1) > s(:,2));
F09 = find (M1(:,2) == 0 & M1(:,4) == 1 & s(:,2) > s(:,1));

sp9 = (numel (M09) + numel (F09)) / size (M, 1) * 100;
fprintf ("Sample Percentage with 90%% certainty: %0.2f%%\n", sp9);

sp9mt = sum (M1(M09,1) == 1) / (numel (M09) + numel (F09)) * 100;
fprintf (" of which %0.1f%% are correctly classified as males and %0.1f%% misclassified,\n", ...
         sp9mt, sum (M1(M09,1) == 2) / (numel (M09) + numel (F09)) * 100);

sp9ft = sum (M1(F09,1) == 2) / (numel (M09) + numel (F09)) * 100;
fprintf (" and %0.1f%% are correctly classified as females and %0.1f%% misclassified.\n\n", ...
         sp9ft, sum (M1(F09,1) == 1) / (numel (M09) + numel (F09)) * 100);

## Remove these samples before assessing next score range
M2 = M1;
M2([M09;F09],:) = [];
M2(:,[2:5]) = [];
s([M09;F09],:) = [];

## Remaining samples that score in FCNN s(:,1) > s(:,2) are males, those which
## score s(:,1) + 0.25 < s(:,2) are females, and the rest are categorized
## according to the number of score instances >=0.8,
## that is when s(:,2) - s(:,1) < 0.25 & s(:,2) - s(:,1) >= 0
s = round (s * 100) / 100;
idx = s(:,2) - s(:,1) < 0.25 & s(:,2) - s(:,1) >= 0;
M08S = find (s(:,1) > s(:,2) | (idx & (M2(:,3) > M2(:,5))));
F08S = find (s(:,2) >= s(:,1) + 0.25 | (idx & (M2(:,5) > M2(:,3))));

sp8 = (numel (M08S) + numel (F08S)) / size (M, 1) * 100;
fprintf ("Sample Percentage with 80%% certainty with FCNN scores: %0.2f%%\n", sp8);

sp8mt = sum (M2(M08S,1) == 1) / (numel (M08S) + numel (F08S)) * 100;
fprintf (" of which %0.1f%% are correctly classified as males and %0.1f%% misclassified,\n", ...
         sp8mt, sum (M2(M08S,1) == 2) / (numel (M08S) + numel (F08S)) * 100);

sp8ft = sum (M2(F08S,1) == 2) / (numel (M08S) + numel (F08S)) * 100;
fprintf (" and %0.1f%% are correctly classified as females and %0.1f%% misclassified.\n\n", ...
         sp8ft, sum (M2(F08S,1) == 1) / (numel (M08S) + numel (F08S)) * 100);

## Calculate and display full model accuracy
fms = sp1 + sp9 * (sp9mt + sp9ft) / 100 + sp8 * (sp8mt + sp8ft) / 100;
fprintf ("Total accuracy for Tibia model: %0.2f%%\n\n", fms);

## Save models into private folder
for i = 1:numel (models)
  name = fullfile (pwd, "private", sprintf ("T%02d.mat", i));
  savemodel (models{i}, name);
endfor
