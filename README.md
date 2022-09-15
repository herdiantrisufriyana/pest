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

Supplement and other files can be found in any of above publications.


## System requirements

We used R 4.0.2 programming language (R Foundation, Vienna, Austria) to conduct 
most steps of the data analysis. To reproduce our work, a set of hardware 
requirements may be needed. We used a single machine with 16 logical processors 
for the 2.10 GHz central processing unit (CPU) (Xeon® E5-2620 v4, Intel®, Santa 
Clara, CA, USA), 128 GB RAM, and 11 GB graphics processing unit (GPU) memory 
(GeForce GTX 1080 Ti, NVIDIA, Santa Clara, CA, USA). Parallelism was applied 
for CPU computing. However, we set the codes to be able for running only 
non-expensive and shorter computation. This only requires 6 GB RAM and 1 CPU.


## Installation guide

Please follow through the R Markdown (pest.Rmd) or R Script (pest.R). 
Installation approximately requires ~5 minutes.


## Demo

All codes require ~20 minutes to complete for non-expensive and shorter 
computation. We provided pre-existing files to substitute the expensive and 
longer computation. One can set a variable that defines whether the expensive 
and longer parts will be conducted. All complete computation may take weeks to 
complete.


## Instructions for use

Briefly, all system requirements, installation guide, demo, and instructions 
for use are available in R Markdown (pest.Rmd) and R Script (pest.R), and 
other files in this repository.

The R Markdown (.Rmd) contains the programming codes for the data analysis. The 
codes for core steps in the analysis pipeline are also provided exclusively in 
an R Script (.R). The codes beyond the core steps were used for analytic 
decision or creating tables or figures. These are shown to provide details on 
how data are processed to construct all tables and figures in the main text.
