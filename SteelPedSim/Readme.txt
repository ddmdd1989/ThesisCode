"Non_CVX_iLL"
* This is a folder including the files to plot the nonconvex figure in the paper
* Please run "NonCVX_Conventional_MAC.m" to plot the figures shown in the nonconvex illustration of MAC value formulation
* "Non_CVX.mat" includes the local optimal result shown in the paper

"Results" 
* This is a folder including the model updating results of the steel pedestrian bridge and the code to plot the result figures in the paper
* Please run "PostProcess.m" to plot the result figures in the paper
* Followings are the descriptions of each result file:
  1. SteelPedBrdg_form1_JACon_LM.mat  - Case 1(a): applying Levenberg-Marquardt algorithm on MAC value formulation with analytical gradient 
  2. SteelPedBrdg_form2_JACon_LM.mat  - Case 2(a): applying Levenberg-Marquardt algorithm on eigenvector difference formulation with analytical gradient 
  3. SteelPedBrdg_form2_JACon_TRR.mat - Case 2(b): applying trust-region-reflective algorithm on eigenvector difference formulation with analytical gradient 

"SAP_Model"
* This is a folder including the sap2000 model of the steel pedestrian bridge

"SteelPedestrianBridge.m" 
* This is the main function to perform FE model updating of steel pedestrian bridge from 100 start points;
* Before run the code, please add "_Shared" folder into MATLAB path. 
* The code first loads the influence matrices of updating variables, the mass matrix and the measured DOF from "SteelPedestrianBridge.mat". Then the code calls the functions in the '"_Shared" to perform FE model updating

"LoadStructure.m"
* This is a function called by "SteelPedestrianBridge.m" to load the influence matrices of updating variables, mass matrix and build the nominal as well as actual stiffness matrices

"SteelPedestrianBridge.mat" 
* This is a file including the influence matrices of updating variables, mass matrix, and measured DOF