## Load required packages
pkg load csg-toolkit
pkg load csg-dataset

## Aggregate CSG data for each bone
longbone_SampleStats
clear

## Process each bone and get cleaned data
Femur_stats; clear
Humerus_stats; clear
Tibia_stats; clear
Ulna_stats; clear

## Compute univariate performance
univar_stats;
clear

## Train/Test models for every bone
Femur_train_models
Femur_test_models
clear
Humerus_train_models
Humerus_test_models
clear
Tibia_train_models
Tibia_test_models
clear
Ulna_train_models
Ulna_test_models
clear

## Test and evaluate final model
sex_test
sex_eval
