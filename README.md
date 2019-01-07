# LabelForest
Semi-supervised learning framework for personalized activity recognition

## Contents
Framework implementation: "LabelForest.R"

Test function: "testCase(filename, numPerLabel)" where "filename" is the name of data file and "numPerLabel" is the number of ground truth data for each class label.

Test output: precision, recall and F1 score averaging over all the classes, as well as the overall labeling accuracy.

Sample data is from a publicly available dataset of 6 activities \[1\], and a laboratory dataset gathered by our research group for 12 activities \[2\].

## Citation 
If you will use this algorithm or the code, please cite the following paper:
@inproceedings{labelforest, 
title={LabelForest: Non-Parametric Semi-Supervised Learning for Activity Recognition}, 
author={Ma, Yuchao and Ghasemzadeh, Hassan}, 
booktitle={The Thirty-Third AAAI Conference on Artificial Intelligence (AAAI'19)},
year={2019}, 
} 

## Reference
\[1\] Jorge-L. Reyes-Ortiz, Luca Oneto, Albert Samà, Xavier Parra, and Davide Anguita. 2016. Transition-Aware Human Activity Recognition Using Smartphones. Neurocomput. 171, C (Jan. 2016), 754–767. https://doi.org/10.1016/j.neucom.2015.07.085
\[2\] R. Fallahzadeh, M. Pedram, and H. Ghasemzadeh. 2016. SmartSock: A wearable platform for context-aware assessment of ankle edema. In 2016 38th Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC). 6302–6306. https://doi.org/10.1109/EMBC.2016.7592169
