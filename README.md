# MSEL_contour_extraction

Author of this Release Package: 
	Yuliang Guo (yuliang_guo@brown.edu)
This is multi-stage approach in extracting curve fragments features from image.
This package includes research code still under development. There are a lot redundent code included.
The evaluation result is not the same reported in the published papers.

Reference: 
	"A Multi-Stage Approach to Curve Extraction", Y.Guo, N.Kumar, M.Narayanan and B.Kimia, ECCV 2014
	"On Evaluating Methods for Recovering Image Curve Fragments", Y.Guo, B.Kimia, CVPRW 2012
    
Setup: 
	Currently, this package of research code only works on Matlab under 64-bit Ubuntu system . The main reason is due the dependency on the core function  'util/mex_unix64/mex_compute_curve_frags'.
    The cross-platform c++ version is still under development.

Prepare edge maps:

	Follow the directory structure in 'Data/TO_SEL_CFGD/'. Construct your own test folder 'Data/****/', prepare your edgemap in the format as our '.edg' file and put them in 'Data/****/edges/'. 
	
	We provide a few third-party edge detectors, but with more accurate orientation modification (third-order orientation correction), check:
    util/multispectral_TO/main_mex_TO_CFGD.m  (recommended) "No grouping left behind: From edges to curve fragments, Tamrakar and Kimia, ICCV 2007"
	util/SE_Dollar_detector/edgesDemo_SETO_CFGD.m    "Structured Forests for Fast Edge Detection, Dollar and Zitnick, ICCV 2013"
	util/gPb_detector/Main_CFGD_gPb_TO.m     "Contour detection and hierarchical image segmentation, Arbelaez etal, PAMI 2011"
	
Compute full set of curve fragments:
	Change the paths in 'demo_compute_SEL'.
	Change the 'Prefix' to 'TO_SEL' or 'gPb_SEL' or 'SE_SEL'depending on your edgemap. Prefix is to load in the trained classifier weights for the corresponding edgemap input.
	
    Tune parameters to control how broken the curve fragments are. The lower these thresholds are, the longer curve fragments generally be.
            params.merge_th = 0.4;
            params.merge_th_geom = 0.4;

Evaluation:
	check 'Eval/main_eval_CFGD_varying_num' with your dataset name (****), and your edge prefix name (####)

	
Evaluation contour maps directly:
	Need to prepare the contour maps in the format of '.cem' in 'Data/****/final_curves/'.
	look at 'util/io/write_cem.m' to understand the format of '.cem' file.
