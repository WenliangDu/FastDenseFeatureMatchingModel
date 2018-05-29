In this folder there are four functions:

	- GTM.m: This function is the implementaion of the GTM algorithm in Matlab. The function takes two sets of corresponding points P1 and P2 (initial
                 matches) and finds matches using GTM algorithm by finding two isomorphic graphs for them. The inputs and outputs are detailed in the following:
		 -Inputs:
			-  P1: an array with size N X 2 containing the first set of points.
			-  P2: an array with size N X 2 containing the corresponding points of the first set of points (the m-th point of P1 and P2 are initially matched).
			-  K: the number of closet neighbores of a point that are connected with point in the graph.
		  -Outputs:
			- inliers: the index of inlier matches in P1 and P2.
			- outliers: the index of outlier matches in P1 and P2.

	- GTM_C_EXE.m: The GTM algorithm is implemented in C++ and its execution file is used in this function for speeding up and utilizing in Matlab functions. 
		       All inputs and outputs are completely similar to those of the function GTM.m detailed above.

	- WGTM.m: This function is the implementaion of the WGTM algorithm in Matlab. The function takes two sets of corresponding points P1 and P2 (initial
                 matches) and finds matches using WGTM algorithm by the geometrical relation between points. The inputs and outputs are detailed in the following:
		 -Inputs:
			-  P1: an array with size N X 2 containing the first set of points.
			-  P2: an array with size N X 2 containing the corresponding points of the first set of points (the m-th point of P1 and P2 are initially matched).
			-  K: the number of closet neighbores of a point that are connected with point in the graph.
			-  emu: is the stopping threshold or the maximum change in the weighted graph after removing an outlier.
		  -Outputs:
			- inliers: the index of inlier matches in P1 and P2.
			- outliers: the index of outlier matches in P1 and P2.

	- WGTM_C_EXE.m: The WGTM algorithm is implemented in C++ and its execution file is used in this function for speeding up and utilizing in Matlab functions. 
		        All inputs and outputs are completely similar to those of the function WGTM.m detailed above.

***************
The source code of both functions in C++ exists in the folder "WGTM & GTM Source code in C++". In this folder, the project named "WGTM_C_Function" contains 
the c++ functions of both GTM and WGTM algorithm. This C++ source code is used to generate two execute file that are used in the function GTM_C_EXE.m and WGTM_C_EXE.m that
are explained above.