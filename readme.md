# Direct and bisulfite-free 5-methylcytosine and 5-hydroxymethylcytosine sequencing at single-cell resolution with scTAPS and scCAPS+ 


Xiufei Chen<sup>1,2,3,10</sup>, Jingfei Cheng<sup>2,3,10</sup>, Linzhen Kong<sup>2,3,10</sup>, Xiao Shu<sup>2,3</sup>, Haiqi Xu<sup>2,3</sup>, Masato Inoue<sup>2,3</sup>, Marion Silvana Fernández Berroca<sup>l4</sup>, Dagny Sanden Døskeland<sup>4</sup>, Magnar Bjørås<sup>4,5,6</sup>, Shivan Sivakumar<sup>7,8</sup>, Yibin Liu<sup>9,* </sup>, Jing Ye<sup>4,* </sup>, Chun-Xiao Song<sup>2,3,* </sup> 

<sup>1</sup> Department of Central Laboratory, The First Affiliated Hospital of Ningbo University, Ningbo 315010, China 

<sup>2</sup> Ludwig Institute for Cancer Research, Nuffield Department of Medicine, University of Oxford, Oxford, UK. OX3 7FZ. 

<sup>3</sup> Target Discovery Institute, Nuffield Department of Medicine, University of Oxford, Oxford, UK. OX3 7FZ. 

<sup>4</sup> Department of Clinical and Molecular Medicine (IKOM), Norwegian University of Science and Technology (NTNU), 7491, Trondheim, Norway.  

<sup>5</sup> Department of Microbiology, Oslo University Hospital, University of Oslo, Oslo, 0424, Norway. 

<sup>6</sup> Centre for Embryology and Healthy Development, University of Oslo, Oslo, 0373, Norway. 

<sup>7</sup> The Institute of Immunology and Immunotherapy, College of Medical and Dental Sciences, University of Birmingham, Edgbaston, Birmingham, B15 2TT. 

<sup>8</sup> Kennedy Institute of Rheumatology, University of Oxford, Oxford, UK. OX3 7FY. 

<sup>9</sup> College of Chemistry and Molecular Sciences, and Taikang Center for Life and Medical Sciences, Wuhan University, Wuhan 430072, China. 

<sup>10</sup> These authors contributed equally. 

<sup>*</sup> Corresponding authors. Email: liuyibin@whu.edu.cn, jing.ye@ntnu.no, and chunxiao.song@ludwig.ox.ac.uk.   


## data preprocessing
preprocess/mESC/configure.yaml    
preprocess/mESC/lowinput_analysis.smk    
preprocess/mESC/mESC_sccaps.smk    

preprocess/t_cell/configure.yaml    
preprocess/t_cell/sctaps.smk    
preprocess/t_cell/standard_taps.smk    

preprocess/mouse_hippocampus/configure.yaml    
preprocess/mouse_hippocampus/nd_taps_convert_call.py    
preprocess/mouse_hippocampus/nd_taps_extract_parallel.py    
preprocess/mouse_hippocampus/sccaps.smk    

## Code for figures
Fig1:
figs/fig1_cd8t_cell_analysis.sh # Bash script to generate input data for plotting
figs/fig1_cd8t_cell.r   # R script to generate the figure
Fig2:
figs/fig2_mouse_hippocampus_analysis.sh # Bash script to generate input data for plotting
figs/fig2_mouse_hippocampus.r # R script to generate the figure
FigS1:
figs/lowinput_archive.r
FigS3:
figs/revision.r
FigS4:
figs/sc_archive.r
figs/revision.r
FigS5:
figs/figure2.r
figs/fig1_cd8t_cell.r
FigS7:
figs/fig2_mouse_hippocampus.r