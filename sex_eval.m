load Sex-Results.mat
Bones = {'Femur', 'Humerus', 'Tibia', 'Ulna'};

for ib = 1:numel (Bones)

  A = [cell2mat(T{2,ib}(:,[2:5])), cellfun(@length, T{2,ib}(:,6))];

  ## No outliers
  a1p = sum (A(:,1) == A(:,2) & A(:,3) == 1 & A(:,4) == 1);
  a1n = sum (A(:,1) != A(:,2) & A(:,3) == 1 & A(:,4) == 1);
  a2p = sum (A(:,1) == A(:,2) & A(:,3) == 2 & A(:,4) == 1);
  a2n = sum (A(:,1) != A(:,2) & A(:,3) == 2 & A(:,4) == 1);
  a3p = sum (A(:,1) == A(:,2) & A(:,3) == 3 & A(:,4) == 1);
  a3n = sum (A(:,1) != A(:,2) & A(:,3) == 3 & A(:,4) == 1);
  ## With outliers
  a4p = sum (A(:,1) == A(:,2) & A(:,3) == 4 & A(:,4) == 0);
  a4n = sum (A(:,1) != A(:,2) & A(:,3) == 4 & A(:,4) == 0);
  a6o = sum (A(:,3) == 5 & A(:,4) == 0);
  a6n = sum (A(:,3) == 5 & A(:,4) == 1);

  L0s = a6o + a6n;
  fprintf ("Total %s samples: %d\n", Bones{ib}, size (A, 1));

  L1s = a1p + a1n;
  L1p = L1s / size (A,1) * 100;
  L1a = a1p / L1s;
  fprintf ("\nLevel 1 contains %d samples: %0.1f%% @ %0.3f\n", L1s, L1p, L1a);
  for iv = 4:-1:1
    vs = (iv - 1) * 8 + 1; ve = iv * 8;
    tr = sum (A(:,3) == 1 & A(:,1) == A(:,2) & A(:,5) >= vs & A(:,5) <= ve);
    tf = sum (A(:,3) == 1 & A(:,1) != A(:,2) & A(:,5) >= vs & A(:,5) <= ve);
    fprintf ("   Range of decisive vars: %d - %d\n", vs, ve);
    fprintf ("   True: %d    False: %d    Accuracy: %0.2f\n", ...
             tr, tf, tr / (tr + tf));
  endfor

  L2s = a2p + a2n;
  L2p = L2s / size (A,1) * 100;
  L2a = a2p / L2s;
  fprintf ("\nLevel 2 contains %d samples: %0.1f%% @ %0.3f\n", L2s, L2p, L2a);
  for iv = 4:-1:1
    vs = (iv - 1) * 8 + 1; ve = iv * 8;
    tr = sum (A(:,3) == 2 & A(:,1) == A(:,2) & A(:,5) >= vs & A(:,5) <= ve);
    tf = sum (A(:,3) == 2 & A(:,1) != A(:,2) & A(:,5) >= vs & A(:,5) <= ve);
    fprintf ("   Range of decisive vars: %d - %d\n", vs, ve);
    fprintf ("   True: %d    False: %d    Accuracy: %0.2f\n", ...
             tr, tf, tr / (tr + tf));
  endfor

  L3s = a3p + a3n;
  L3p = L3s / size (A,1) * 100;
  L3a = a3p / L3s;
  fprintf ("\nLevel 3 contains %d samples: %0.1f%% @ %0.3f\n", L3s, L3p, L3a);
  for iv = 4:-1:1
    vs = (iv - 1) * 8 + 1; ve = iv * 8;
    tr = sum (A(:,3) == 3 & A(:,1) == A(:,2) & A(:,5) >= vs & A(:,5) <= ve);
    tf = sum (A(:,3) == 3 & A(:,1) != A(:,2) & A(:,5) >= vs & A(:,5) <= ve);
    fprintf ("   Range of decisive vars: %d - %d\n", vs, ve);
    fprintf ("   True: %d    False: %d    Accuracy: %0.2f\n", ...
             tr, tf, tr / (tr + tf));
  endfor

  L4s = a4p + a4n;
  L4p = L4s / size (A,1) * 100;
  L4a = a4p / L4s;
  fprintf ("\nLevel 4 contains %d samples: %0.1f%% @ %0.3f\n", L4s, L4p, L4a);
  for iv = 4:-1:1
    vs = (iv - 1) * 8 + 1; ve = iv * 8;
    tr = sum (A(:,3) == 4 & A(:,1) == A(:,2) & A(:,5) >= vs & A(:,5) <= ve);
    tf = sum (A(:,3) == 4 & A(:,1) != A(:,2) & A(:,5) >= vs & A(:,5) <= ve);
    fprintf ("   Range of decisive vars: %d - %d\n", vs, ve);
    fprintf ("   True: %d    False: %d    Accuracy: %0.2f\n", ...
             tr, tf, tr / (tr + tf));
  endfor

  L0p = L0s / size (A,1) * 100;
  fprintf ("\nNumber of undetermined samples: %d (%0.1f%%)\n", L0s, L0p);
  fprintf ("Total classification accuracy: %0.1f\n\n", ...
           L1p * L1a + L2p * L2a + L3p * L3a + L4p * L4a);

endfor
