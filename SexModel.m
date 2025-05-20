## Implementation of the final multi-level voting-based hierarchical prediction
## model, which is ported into the longbone_Sex function of the csg-toolkit

function [varargout] = SexModel (filename, BONE, DATA)

  ## Load measurements
  measurements = longbone_Measurements ();

  ## Select variables for sex estimation
  if (strcmpi (BONE, 'Femur'))
    [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives ('Femur');
    X = DATA(idx);
    Z = (X - mu) ./ sigma;
    measurements = measurements(idx);
    prob = {[1, 0.98, 0.94, 0.95], 0.62, 0.75, [1, 1, 0.8, 0.67], NaN};
  elseif (strcmpi (BONE, 'Humerus'))
    [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives ('Humerus');
    X = DATA(idx);
    Z = (X - mu) ./ sigma;
    measurements = measurements(idx);
    prob = {[0.97, 0.96, 0.95, 0.97], 0.67, 0.56, [0.96, 0.63, 0.63, 0.63], NaN};
  elseif (strcmpi (BONE, 'Tibia'))
    [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives ('Tibia');
    X = DATA(idx);
    Z = (X - mu) ./ sigma;
    measurements = measurements(idx);
    prob = {[0.99, 0.95, 0.92, 0.94], 0.49, 0.65, [1, 0.67, 0.67, 0.67], NaN};
  elseif (strcmpi (BONE, 'Ulna'))
    [idx, mu, sigma, Mlo, Mhi, Flo, Fhi] = load_descriptives ('Ulna');
    X = DATA(idx);
    Z = (X - mu) ./ sigma;
    measurements = measurements(idx);
    prob = {[0.99, 0.97, 0.96, 0.92], 0.56, 0.68, [0.96, 0.5 0.5, 0.5], NaN};
  else
    sex = categorical (0, [1, 2], {'Male', 'Female'}); # <undefined>
    vars = {};
    outliers = {};
    if (nargout == 1 || nargout == 0)
      varNames = {'Name', 'Bone', 'Sex', 'Prob', 'Typical', 'Used', 'Outliers'};
      T = table ({filename}, {BONE}, sex, NaN, false, {vars}, {outliers}, ...
                 'VariableNames', varNames);
      varargout{1} = T;
    else
      varargout{1} = double (sex);
      if (nargout > 1)
        varargout{2} = 5;
      endif
      if (nargout > 2)
        varargout{3} = false;
      endif
      if (nargout > 3)
        varargout{4} = vars;
      endif
      if (nargout > 4)
        varargout{5} = outliers;
      endif
    endif
    return;
  endif

  ## Load classifiers and compute prediction scores
  models = load_classifiers (BONE);
  Scores = [];
  [~, score] = predict (models{1}, Z);
  Scores = [Scores, score];
  [~, score] = predict (models{2}, Z);
  Scores = [Scores, score];
  idx = 3;
  for i = 1:31
    [~, score] = predict (models{idx}, Z(i));
    Scores = [Scores, score];
    idx += 1;
    [~, score] = predict (models{idx}, Z(i));
    Scores = [Scores, score];
    idx += 1;
  endfor
  [~, s] = predict (models{65}, Scores);

  ## For samples without any outliers in measurements
  Mrep = all (X > Mlo & X < Mhi);
  Frep = all (X > Flo & X < Fhi);
  if (Mrep || Frep)
    rep = true;
    ## Group scores according to their magnitude
    odd_cols = Scores([1:2:end]);
    evencols = Scores([2:2:end]);
    m10s = any (odd_cols == 1 & evencols != 1);
    f10s = any (evencols == 1 & odd_cols != 1);
    m09s = any (odd_cols >= 0.9 & odd_cols < 1);
    f09s = any (evencols >= 0.9 & evencols < 1);
    m08n = numel (find (odd_cols >= 0.8 & odd_cols < 0.9));
    f08n = numel (find (evencols >= 0.8 & evencols < 0.9));
    odd_LDAs = odd_cols(3:2:end);
    evenLDAs = evencols(3:2:end);
    m09 = odd_cols >= 0.9;
    f09 = evencols >= 0.9;
    m08 = odd_cols >= 0.8;
    f08 = evencols >= 0.8;

    ## Apply hierarchical decision model
    if (m10s && s(1) > s(2))
      sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
      grp = 1;
      idx = find (odd_cols == 1);
      vars = resolve_idx (idx, measurements);
      idx = find_prob_idx (vars);
    elseif (f10s && s(1) < s(2))
      sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
      grp = 1;
      idx = find (evencols == 1);
      vars = resolve_idx (idx, measurements);
      idx = find_prob_idx (vars);
    elseif (strcmpi (BONE, 'Femur') && sum (m09) > sum (f09))
      sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
      grp = 2;
      idx = find (odd_cols >= 0.9 & odd_cols < 1);
      vars = resolve_idx (idx, measurements);
      idx = 1;
    elseif (strcmpi (BONE, 'Femur') && sum (m09) < sum (f09))
      sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
      grp = 2;
      idx = find (evencols >= 0.9 & evencols < 1);
      vars = resolve_idx (idx, measurements);
      idx = 1;
    elseif (! strcmpi (BONE, 'Femur') && (m09s && ! f09s && s(1) > s(2)))
      sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
      grp = 2;
      idx = find (odd_cols >= 0.9 & odd_cols < 1);
      vars = resolve_idx (idx, measurements);
      idx = 1;
    elseif (! strcmpi (BONE, 'Femur') && (! m09s && f09s && s(1) < s(2)))
      sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
      grp = 2;
      idx = find (evencols >= 0.9 & evencols < 1);
      vars = resolve_idx (idx, measurements);
      idx = 1;
    else
      s = round (s * 100) / 100;
      gap = s(2) - s(1) < 0.25 & s(2) - s(1) >= 0;
      if (s(1) > s(2) || (gap && (m08n > f08n)))
        sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
        grp = 3;
        vars = {};
        if (gap && (m08n > f08n))
          idx = find (odd_cols >= 0.8 & odd_cols < 0.9);
          vars = [vars, resolve_idx(idx, measurements)];
        endif
        if (s(1) > s(2))
          vars = [vars, {'FCNN Scores'}];
        endif
      elseif (s(2) >= s(1) + 0.25 || (gap && (f08n > m08n)))
        sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
        grp = 3;
        vars = {};
        if (gap && (f08n > m08n))
          idx = find (evencols >= 0.8 & evencols < 0.9);
          vars = [vars, resolve_idx(idx, measurements)];
        endif
        if (s(2) >= s(1) + 0.25)
          vars = [vars, {'FCNN Scores'}];
        endif
      else
        sex = categorical (0, [1, 2], {'Male', 'Female'}); # <undefined>
        grp = 5;
        vars = {};
      endif
      idx = 1;
    endif
    outliers = {};

  else
    rep = false;
    ## For samples with outliers work only with univariate LDA models
    odd_cols = Scores([1:2:end]);
    evencols = Scores([2:2:end]);
    odd_cols = odd_cols(3:2:end);
    evencols = evencols(3:2:end);
    m09 = odd_cols > 0.9;
    f09 = evencols > 0.9;
    if (sum (m09) > sum (f09))
      sex = categorical (1, [1, 2], {'Male', 'Female'}); # Male
      grp = 4;
      vars = measurements(find (m09));
      idx = find_prob_idx (vars);
      outliers = measurements(m09  & (X < Mlo | X > Mhi));
    elseif (sum (m09) < sum (f09))
      sex = categorical (2, [1, 2], {'Male', 'Female'}); # Female
      grp = 4;
      vars = measurements(find (f09));
      idx = find_prob_idx (vars);
      outliers = measurements(f09 & (X < Flo | X > Fhi));
    else
      sex = categorical (0, [1, 2], {'Male', 'Female'}); # <undefined>
      grp = 5;
      vars = {};
      idx = 1;
      outliers = {};
    endif
  endif

  ## Prepare output
  if (nargout == 1 || nargout == 0)
    prob = prob{grp}(idx);
    varNames = {'Name', 'Bone', 'Sex', 'Prob', 'Typical', 'Used', 'Outliers'};
    T = table ({filename}, {BONE}, sex, prob, rep, {vars}, {outliers}, ...
               'VariableNames', varNames);
    varargout{1} = T;
    return;
  endif
  varargout{1} = double (sex);
  if (nargout > 1)
    varargout{2} = grp;
  endif
  if (nargout > 2)
    varargout{3} = rep;
  endif
  if (nargout > 3)
    varargout{4} = vars;
  endif
  if (nargout > 4)
    varargout{5} = outliers;
  endif

endfunction

function used = resolve_idx (idx, measurements)
  used = {};
  if (idx(1) == 1)
    used = [used, {'LDA All'}];
    idx(1) = [];
  endif
  if (! isempty (idx))
    if (idx(1) == 2)
      used = [used, {'KNN All'}];
      idx(1) = [];
    endif
  endif
  if (! isempty (idx))
    if (idx(1) == 3)
      idx(1) = 4;
    endif
    used = [used, measurements(unique (floor ((idx - 2) / 2)))];
  endif
endfunction

function idx = find_prob_idx (vars)
  nvars = length (vars);
  if (nvars <= 8)
    idx = 4;
  elseif (nvars <= 16)
    idx = 3;
  elseif (nvars <= 24)
    idx = 2;
  elseif (nvars <= 32)
    idx = 1;
  endif
endfunction
