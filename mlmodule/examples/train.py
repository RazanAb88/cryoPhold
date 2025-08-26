import argparse
from timelaggedcv import TimeLaggedCV

def main(pickle_descriptor, pickle_features, tau):
    options = TimeLaggedCV.DEFAULT
    options['create_dataset_options']['pickle_descriptor'] = pickle_descriptor
    options['create_dataset_options']['pickle_features'] = pickle_features
    options['dataset']['tau'] = tau
    model = TimeLaggedCV(options)
    model.run()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='TimeLaggedCV command line tool with configurable pickle files and tau parameter'
    )
    parser.add_argument(
        '--pickle_descriptor',
        nargs='+',  # Accept one or more file paths
        type=str,
        default=[
            './data-fgfr2/AC-loop/dfsincos.pkl',
            './data-fgfr2/AC-helix/dfsincos.pkl'
        ],
        help='One or more file paths for the pickle descriptor (e.g., file1.pkl [file2.pkl ...])'
    )
    parser.add_argument(
        '--pickle_features',
        nargs='+',  # Accept one or more file paths
        type=str,
        default=[
            './data-fgfr2/AC-loop/featuressincos.pkl',
            './data-fgfr2/AC-helix/featuressincos.pkl'
        ],
        help='One or more file paths for the pickle features (e.g., file1.pkl [file2.pkl ...])'
    )
    parser.add_argument(
        '--tau',
        type=int,
        default=10,
        help='Time lag parameter (integer)'
    )
    args = parser.parse_args()

    # The argparse with nargs='+' ensures that even a single file is wrapped in a list.
    main(args.pickle_descriptor, args.pickle_features, args.tau)
