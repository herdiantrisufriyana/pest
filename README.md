# Blood biomarkers representing maternal-fetal interface tissues used to predict early-and late-onset preeclampsia but not COVID-19 infection

Herdiantri Sufriyana,a,b Hotimah Masdan Salim,c Akbar Reza Muhammad,d Yu-Wei 
Wu,a,e Emily Chia-Yu Sua,e,f, *

a Graduate Institute of Biomedical Informatics, College of Medical Science and 
Technology, Taipei Medical University, 250 Wu-Xing Street, Taipei 11031, Taiwan.

b Department of Medical Physiology, Faculty of Medicine, Universitas Nahdlatul 
Ulama Surabaya, 57 Raya Jemursari Road, Surabaya 60237, Indonesia.

c Department of Molecular Biology, Faculty of Medicine, Universitas Nahdlatul 
Ulama Surabaya, 57 Raya Jemursari Road, Surabaya 60237, Indonesia.

d Faculty of Medicine, Universitas Nahdlatul Ulama Surabaya, 57 Raya Jemursari 
Road, Surabaya 60237, Indonesia.

e Clinical Big Data Research Center, Taipei Medical University Hospital, 250 
Wu-Xing Street, Taipei 11031, Taiwan.

f Research Center for Artificial Intelligence in Medicine, Taipei Medical 
University, 250 Wu-Xing Street, Taipei 11031, Taiwan.

\* Corresponding author: Clinical Big Data Research Center, Taipei Medical 
University Hospital, 250 Wu-Xing Street, Taipei 11031, Taiwan. Phone: 
+886-2-66382736 ext. 1515. Email address: emilysu@tmu.edu.tw.

The paper can be found here:
https://doi.org/10.1016/j.csbj.2022.08.011

Supplement and other files can be found in the paper.


## System requirements

We used R 4.0.2 programming language (R Foundation, Vienna, Austria) to conduct 
most steps of the data analysis. To reproduce our work, a set of hardware 
requirements may be needed. We used a single machine with 16 logical processors 
for the 2.10 GHz central processing unit (CPU) (Xeon® E5-2620 v4, Intel®, Santa 
Clara, CA, USA), 128 GB RAM, and 11 GB graphics processing unit (GPU) memory 
(GeForce GTX 1080 Ti, NVIDIA, Santa Clara, CA, USA). Parallelism was applied 
for CPU computing. However, we set the codes to be able for running only 
non-expensive and shorter computation. This only requires 3.6 GB RAM and 1 CPU.


## Installation guide

This analytical pipeline needs all files in this repository and Google Drive. Please manually dowload the latter files (17.47 GB) from this link:

https://drive.google.com/file/d/12sIGr0ys07WyMNDdZ55hmhTzF7PRoBA4/view?usp=sharing

Put the zip file from Google Drive in the same folder with pest.Rmd and pest.R; 
thus, the codes can unzip the files to the correct path. Docker installation is 
also provided, but one needs to manually download the zip file from Google 
Drive and upload it to the container. For the first time, it is suggested to 
run the codes sequentially from top to bottom. Installation approximately 
requires ~5 minutes.


## Demo

All codes require ~25 minutes to complete for non-expensive and shorter 
computation. We provided pre-existing files to substitute the expensive and 
longer computation. All complete computation may take weeks to 
complete. One can set a variable in "Set to run or not run very heavy 
computations" that defines whether the expensive and longer parts will be 
conducted:

Set run_heavy_computation=T to run the codes which require a computer with 
the system requirements above. Load data/log.rds to get details on how long 
each piece of codes need to be completed. Ensure to set other variables below 
to TRUE.

Set load_outlier_comp=T to find the steps for defining outliers.

Set load_trained_divnn=T to load trained DI-VNN models.

Set load_emulators=T to load the emulated tree models.

Set show_best_tree=T to show the best tree with different datasets.

## Instructions for use

The R Markdown (.Rmd) contains the programming codes for the data analysis. The 
codes for core steps in the analysis pipeline are also provided exclusively in 
an R Script (.R). The codes beyond the core steps were used for analytic 
decision or creating tables or figures. These are shown to provide details on 
how data are processed to construct all tables and figures in the main text.
