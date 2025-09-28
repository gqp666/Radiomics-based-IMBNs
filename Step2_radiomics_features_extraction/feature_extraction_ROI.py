import six
import os
import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor
from scipy.io import savemat
import time
from multiprocessing import Pool, cpu_count  # Import for parallel processing
import logging
# set level for all classes
logger = logging.getLogger("radiomics")
logger.setLevel(logging.ERROR)
# ... or set level for specific class
logger = logging.getLogger("radiomics.glcm")
logger.setLevel(logging.ERROR)

def initialize_extractor():
    """Initializes and returns a PyRadiomicsFeatureExtractor."""
    settings = {}
    settings['sigma'] = [ 3, 5]
    settings['binCount'] = 80
    settings['normalize'] = False
    settings['normalizeScale'] = 1

    extractor = featureextractor.RadiomicsFeatureExtractor(**settings)

    # Disable all and enable specific image types for efficiency
    extractor.disableAllImageTypes()
    extractor.enableImageTypeByName('Original')
    extractor.enableImageTypeByName('Wavelet')
    extractor.enableImageTypeByName('LoG')
    extractor.enableAllFeatures()

    return extractor


def extract_features_for_label(args):
    """
    Extracts radiomics features for a given image, mask, and label.
    This function is designed to be called by a multiprocessing pool.
    """
    image_path, mask_path, label_reg = args

    # Initialize extractor inside the child process to avoid pickling issues
    # and ensure each process has its own extractor instance.
    extractor = initialize_extractor()

    print(f"  Processing ROI: {label_reg} for {os.path.basename(os.path.dirname(image_path))}")

    try:
        result = extractor.execute(image_path, mask_path, label=label_reg)
    except Exception as e:
        print(f"Error extracting features for label {label_reg} in {image_path}: {e}")
        return label_reg, None, None  # Return None for data and names if extraction fails

    feature_cur = []
    feature_name = []
    for key, value in six.iteritems(result):
        feature_name.append(key)
        feature_cur.append(value)

    # Ensure feature_cur values are floats and get the relevant subset
    processed_features = [float(f) for f in feature_cur[36:]]
    # Get the corresponding feature names
    processed_names = feature_name[36:]

    return label_reg, processed_features, processed_names


# 解析输入的区间，单个区域后跟英文逗号，区间用:连接，如3:5,左右区间都包含
def parse_reg_values(reg_values):
    # 解析 reg_values，支持单个数、多个数、一个区间、多个区间的组合
    result = []
    for item in reg_values.split(','):
        item = item.strip()  # Remove whitespace
        if ':' in item:
            start, end = map(int, item.split(':'))
            result.extend(range(start, end + 1))
        else:
            result.append(int(item))
    return result


if __name__ == "__main__":

    reg_values_str = "1:246"

    # Determine the number of CPU cores to use for multiprocessing
    num_processes = 7  # Use all but one core to avoid freezing your system
    print(f"Using {num_processes} processes for parallel feature extraction.")
    parsed_reg_values = parse_reg_values(reg_values_str)

    base_data_paths_list = [
        r"E:\Data\ASD\ASD",
        r"E:\Data\ASD\TD",
    ]
    output_base_folders_list = [

        r"D:\project\pycharm\brain\ROI\project_1\ASD\ASD",
        r"D:\project\pycharm\brain\ROI\project_1\ASD\TD",
    ]
    prefix_list = [
        'mwp1',
        'mwp1',
    ]

    mask_path_list = [

        r'D:\project\pycharm\brain\atlas\BN_Atlas_246\BN_Atlas_246_1.5mm.nii',
        r'D:\project\pycharm\brain\atlas\BN_Atlas_246\BN_Atlas_246_1.5mm.nii',
    ]

    for base_data_path, output_base_folder, prefix, mask_path in zip(base_data_paths_list, output_base_folders_list, prefix_list, mask_path_list ):
        print(f"  Processing {base_data_path} for {prefix}")
        if not os.path.exists(output_base_folder):
            os.makedirs(output_base_folder)
        # Traverse directories and process files
        for first_level_dir in os.listdir(base_data_path):
            first_level_dir_path = os.path.join(base_data_path, first_level_dir)
            mri_dir_path = os.path.join(first_level_dir_path, 'mri')

            for filename in os.listdir(mri_dir_path):
                if filename.startswith(prefix) and filename.endswith('.nii'):
                    file_path = os.path.join(mri_dir_path, filename)


                    start_time = time.time()
                    print(f"Processing folder: {first_level_dir_path}")

                    # Define the output folder for the current subject, preserving the relative path

                    current_subject_output_folder = os.path.join(output_base_folder, first_level_dir)
                    os.makedirs(current_subject_output_folder, exist_ok=True)

                    # Prepare arguments for parallel processing
                    tasks = [(file_path, mask_path, label_reg) for label_reg in parsed_reg_values]

                    reg_id_list = []
                    ft_list_data = []
                    name_features = None  # To store feature names (they should be consistent across labels)

                    # Use a multiprocessing Pool to extract features in parallel for different labels
                    with Pool(processes=num_processes) as pool:
                        results = pool.map(extract_features_for_label, tasks)

                    for label_reg, save_curdata, current_name_features in results:
                        if save_curdata is not None:
                            reg_id_list.append(label_reg)
                            ft_list_data.append(save_curdata)
                            if name_features is None:  # Only set feature names once
                                name_features = current_name_features

                    if ft_list_data:  # Only proceed if features were successfully extracted
                        # Create DataFrame, set 'LabelReg' as index
                        # Ensure name_features is not None before using it
                        if name_features is None:
                            print(f"Warning: No feature names found for {first_level_dir}. Skipping DataFrame creation.")
                            continue

                        name_df = pd.DataFrame(data=np.column_stack([reg_id_list, ft_list_data]),
                                               columns=['LabelReg'] + list(name_features))
                        name_df.set_index('LabelReg', inplace=True)

                        # Transpose DataFrame
                        name_df_transposed = name_df.T

                        # Create filename: 'Radiomics' + last sub-directory name + '.mat'

                        mat_filename = f'Radiomics_{first_level_dir}.mat'
                        mat_output_path = os.path.join(current_subject_output_folder, mat_filename)

                        savemat(mat_output_path, {'data': name_df_transposed.to_numpy()})
                        print(f"Saved .mat file to: {mat_output_path}")
                    else:
                        print(f"No features extracted for {first_level_dir}. Skipping .mat file creation.")

                    end_time = time.time()
                    time_taken = end_time - start_time
                    print(f"Time taken to process {first_level_dir}: {time_taken:.2f} seconds\n")



