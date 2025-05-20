## Calculate descriptive statistics for sample sizes, age, and sex per bone,
## side, and collection.  Aggregate data per bone into numeric arrays by
## substituting filenames with collection indices.

Collections = {'Athens', 'Crete', 'Granada'};
Bones = {'Femur', 'Humerus', 'Tibia', 'Ulna'};
Side = {'Left', 'Right'};
Sex = {'Males', 'Females'};
Header = {'Sample ID', 'Side', 'Collections', 'Sex', 'Age', 'Max Distance', ...
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

## load data
data = load ('CSG-Data.mat').T;

## useful function handles
nametoid = @(x) str2num (x([4:6]));
sidetoid = @(x) double (x(end-4));

LR = table ('Size', [4, 2], 'VariableTypes', {'double', 'double'}, ...
            'VariableNames', Side);

SEX = table ('Size', [5, 2], 'VariableTypes', {'double', 'double'}, ...
            'VariableNames', Sex);

AGE = table ('Size', [3, 3], 'VariableTypes', {'double', 'double', 'double'}, ...
             'VariableNames', {'Min', 'Median', 'Max'});

T1 = table (LR, LR, LR, LR, 'VariableNames', Bones, 'RowNames', [Collections, {'Total'}]);
T2 = table (SEX, SEX, SEX, 'VariableNames', Collections, 'RowNames', [Bones, {'Individuals'}]);
T3 = table (AGE, AGE, AGE, 'VariableNames', Collections, 'RowNames', [Sex, {'All'}]);

SSEX = {};
CSEX = {};
CAGE = {};
Femur = [];
Humerus = [];
Tibia = [];
Ulna = [];
ageIdx = 1;

for ic = 1:size (data, 2)
  header = data{1, ic};
  hnames = strsplit (header);
  subset = data{2, ic};
  subset = cell2mat (subset(2:end,2:end));
  snames = data{2, ic}(2:end,1);

  nid = cellfun (nametoid, snames);
  sid = cellfun (sidetoid, snames);
  sid(sid == 76) = 1; sid(sid == 82) = 2;
  cid = ones (size (subset, 1), 1) * find (strcmpi (hnames{1}, Collections));
  if (strcmpi (hnames{2}, 'Femur'))
    Femur = [Femur; nid, sid, cid, subset];
  elseif (strcmpi (hnames{2}, 'Humerus'))
    Humerus = [Humerus; nid, sid, cid, subset];
  elseif (strcmpi (hnames{2}, 'Tibia'))
    Tibia = [Tibia; nid, sid, cid, subset];
  else
    Ulna = [Ulna; nid, sid, cid, subset];
  endif

  idx = strcmpi (hnames{1}, Collections);
  T1.(hnames{2}).(hnames{3})(idx) = size (subset, 1);

  SSEX = [SSEX; data{2, ic}(2:end,1:3)];
  if (mod (ic, 2) == 0)
    names = SSEX(:,1);
    [~, idx] = unique (names);
    sex = cell2mat (SSEX(:,2));
    sex = sex(idx);
    idx = strcmpi (hnames{2}, Bones);
    T2.(hnames{1}).Males(idx) = sum (sex == 1);
    T2.(hnames{1}).Females(idx) = sum (sex == 2);
    CSEX = [CSEX; SSEX];

    if (mod (ic, 8) == 0)
      idnums = cell2mat (cellfun (nametoid, CSEX(:,1), "UniformOutput", false));
      [~, idx] = unique (idnums);
      sex = cell2mat (CSEX(:,2));
      sex = sex(idx);
      T2.(hnames{1}).Males(5) = sum (sex == 1);
      T2.(hnames{1}).Females(5) = sum (sex == 2);

      ## Handle age here
      age = cell2mat (CSEX(:,3));
      age = age(idx);
      T3.(hnames{1}).Min(1) = nanmin (age(sex==1));
      T3.(hnames{1}).Min(2) = nanmin (age(sex==2));
      T3.(hnames{1}).Min(3) = nanmin (age);
      T3.(hnames{1}).Median(1) = round (median (age(sex==1), 'omitnan'));
      T3.(hnames{1}).Median(2) = round (median (age(sex==2), 'omitnan'));
      T3.(hnames{1}).Median(3) = round (median (age, 'omitnan'));
      T3.(hnames{1}).Max(1) = nanmax (age(sex==1));
      T3.(hnames{1}).Max(2) = nanmax (age(sex==2));
      T3.(hnames{1}).Max(3) = nanmax (age);
      CAGE(ageIdx) = {age(sex==1)};
      ageIdx += 1;
      CAGE(ageIdx) = {age(sex==2)};
      ageIdx += 1;
      CSEX = {};
    endif
    SSEX = {};
  endif
endfor

## Add sums in T1
for ib = Bones
  for sd = Side
    T1.(ib{:}).(sd{:})(4) = sum (T1.(ib{:}).(sd{:}));
  endfor
endfor

disp (T1);
disp (T2);
disp (T3);

save ('-binary', 'Bone-Data.mat', 'Header', 'Femur', 'Humerus', 'Tibia', 'Ulna');

## Print figure with age distributions per sex and collection
Titles = {"Athens males", "Crete males", "Granada males", ...
          "Athens females", "Crete females", "Granada females"};
hcolor = {"cyan", "cyan", "cyan", "magenta", "magenta", "magenta"};
ageIdx = [1, 3, 5, 2, 4, 6];
hf = figure ("position", [400 0 1200 800]);
for i = 1:6
  subplot (2, 3, i);
  data = CAGE{ageIdx(i)};
  data(isnan (data)) = [];
  hist (data, [18:12:102] + 6, hcolor{i});
  set (gca, "xtick", [18 30 42 54 66 78 90 102], "fontsize", 14);
  xlim ([15, 105]);
  ylim ([0, 40]);
  xlabel ("age-at-death", "fontsize", 16, "fontangle", "italic");
  ylabel ("individuals", "fontsize", 16, "fontangle", "italic");
  title (Titles{i}, "fontsize", 18);
  hl1 = line (min (data) * [1, 1], [0, 40], "linestyle", "-.", "linewidth", 2, ...
                                   "color", "red", "displayname", "min");
  hl2 = line (median (data) * [1, 1], [0, 40], "linestyle", "-.", "linewidth", 2, ...
                                   "color", "green", "displayname", "median");
  hl3 = line (max (data) * [1, 1], [0, 40], "linestyle", "-.", "linewidth", 2, ...
                                   "color", "blue", "displayname", "max");
  legend ([hl1, hl2, hl3], "location", "northwest");
endfor
print (hf, "age.png");
