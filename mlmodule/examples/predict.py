import argparse
from timelaggedcv import TimeLaggedCV
import numpy as np
import matplotlib.pyplot as plt

def main(pickle_descriptor, pickle_features, tau, teacher_flag):
    options = TimeLaggedCV.DEFAULT
    options['create_dataset_options']['pickle_descriptor'] = pickle_descriptor
    options['create_dataset_options']['pickle_features'] = pickle_features
    options['dataset']['tau'] = tau
    options['general']['teacher'] = teacher_flag

    tlcv = TimeLaggedCV(options)
    model = tlcv.get_model()

    if not teacher_flag:
        model.load('model_final.pt')
    else:
        model.load('model.pt')

    cv_space = model.estimate_cv_from_dataset(tlcv.dataset_weight, batchsize=100000)
    list_cv = np.split(cv_space, tlcv.dataset_weight.total_length)

    if not teacher_flag:
        np.savez('cv_pred.npz', *list_cv)
        print('Saved teacher prediction')
    else:
        np.savez('cv_pred_student.npz', *list_cv)
        print('Saved student prediction')

    skip = 1
    for start, end in zip(tlcv.dataset_weight.total_length[:-1], tlcv.dataset_weight.total_length[1:]):
        plt.plot(*cv_space[start:end:skip].T, '.')
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run TimeLaggedCV with specified pickle descriptors, features, and tau.'
    )
    parser.add_argument(
        '--pickle_descriptor',
        type=str,
        nargs='+',
        required=True,
        help='Space-separated list of paths for pickle_descriptor'
    )
    parser.add_argument(
        '--pickle_features',
        type=str,
        nargs='+',
        required=True,
        help='Space-separated list of paths for pickle_features'
    )
    parser.add_argument(
        '--tau',
        type=int,
        required=True,
        help='Integer value for tau'
    )
    parser.add_argument(
        '--teacher_flag',
        action='store_true',
        help='Set this flag if teacher mode is desired (default is False)'
    )

    args = parser.parse_args()

    main(args.pickle_descriptor, args.pickle_features, args.tau, args.teacher_flag)
