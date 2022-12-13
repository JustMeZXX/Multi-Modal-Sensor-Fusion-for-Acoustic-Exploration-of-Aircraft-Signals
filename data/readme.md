### Dataset
The dataset used in this paper can be downloaded from: [Dataset](https://smu.box.com/s/ka1oo3mrz1oamgie0fcdyvz7xm3puwif). 

The dataset consists of both acoustic and ads-b data from Julian day 217 to 247 of 2021, which are collected on the campus of Southern Methodist University (SMU). 

The labels for the acoustic aircraft detection are saved as .mat files in folder "./labels". The labeled aircraft are denoted as a 1x31 cell, where each cell records for one Julian day. At each cell, the label is a Mx4 matrix (M is the number of observations). For each observation, the first two columns are the start and end time, the third column is the source frequency, and the fourth column is the time at which the source frequency is measured. 
