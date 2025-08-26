python train.py \
  --pickle_descriptor \
    /Users/dihedral_features/feature_descriptions.pkl \
  --pickle_features \
    /Users/dihedral_features/best_frames_features.pkl \
  --tau 1

 python predict.py \
  --pickle_descriptor \
    /Users/dihedral_features/feature_descriptions.pkl \
  --pickle_features \
    /Users/dihedral_features/best_frames_features.pkl \
  --tau 1
