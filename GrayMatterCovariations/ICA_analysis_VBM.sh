#!/bin/bash

#author:T Mei,last edited:18-09-2023

datapath=/.../VBM_ICA

#generate a 4d VBM file
fslmerge -t ${datapath}/4d_VBM.nii.gz `cat ${datapath}/VBM_sub.txt` #the VBM data path of each participant
#run Melodic ICA on this 4d file
melodic -i ${datapath}/4d_VBM.nii.gz -o ${datapath}/results_ICA_100 -d 100 --Oall --report --nobet 




