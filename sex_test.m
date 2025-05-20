## Process all 3D models and aggregate results per bone, side, and collection
Bones = {'Femur', 'Humerus', 'Tibia', 'Ulna'};

T = cell (2, 4);
for ib = 1:numel (Bones)
  ## Load data for specific bone
  DATA = load ('Bone-Data.mat').(Bones{ib});

  ## Process all samples of specific bone
  data = {};
  for is = 1:size (DATA, 1)
    ## Generate name for sample
    filename = sprintf ("%d_%d_%d", DATA(is,3), DATA(is,1), DATA(is,2));

    ## Analyze data of sample
    [out, grp, rep, vars, outliers] = SexModel (filename, Bones{ib}, DATA(is, [6:66]));

    ## Find sex for sample
    sex = DATA(is, 4);

    ## Aggregate results
    data = [data; filename, {sex, out}, {grp}, {rep}, {vars}, {outliers}];
  endfor

  ## Append to cell array
  T(1, ib) = Bones{ib};
  T(2, ib) = {data};

endfor

## Save cell array with all collected data
save ('-binary', 'Sex-Results.mat', 'T')
