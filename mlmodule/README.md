## Installation
We recommend to use **[pixi](https://pixi.sh/dev/)** as package manager. Installation instruction can be found on their website. 
As soon as you have pixi, you can simply run the following command in this folder (where the pyproject.toml file is located):
```shell
pixi install
```

After that you can run 
```shell
pixi shell
```
and you have "activated" your local environment.By changing to the example folder you can then execute the training example by

```shell
python train.py \
  --pickle_descriptor \
    /Users/dihedral_features/feature_descriptions.pkl \
  --pickle_features \
    /Users/dihedral_features/best_frames_features.pkl \
  --tau 1
```
and have a look at the resulting files.
For further investigating the result, you can have a look at the predict.py file, which loads the model and makes a prediction for the provided files.
```shell
python predict.py \
  --pickle_descriptor \
    /Users/dihedral_features/feature_descriptions.pkl \
  --pickle_features \
    /Users/dihedral_features/best_frames_features.pkl \
  --tau 1
```
In the current state this should provide you with the prediction of the CV space, where each trajectory is colored differently.
