"Results" 
* This is a folder including the model updating results of the concrete building frame and the code to plot the result figures in the paper
* Please run "PostProcess.m" to plot the result figures in the paper
* Followings are the descriptions of each result file:
  1. ConcBuildFrm_form1_JACoff_LM.mat  - Case 1(a): applying Levenberg-Marquardt algorithm on MAC value formulation with numerical gradient 
  2. ConcBuildFrm_form1_JACon_LM.mat   - Case 1(a): applying Levenberg-Marquardt algorithm on MAC value formulation with analytical gradient 
  3. ConcBuildFrm_form2_JACoff_LM.mat  - Case 2(a): applying Levenberg-Marquardt algorithm on eigenvector difference formulation with numerical gradient 
  4. ConcBuildFrm_form2_JACoff_TRR.mat - Case 2(b): applying trust-region-reflective algorithm on eigenvector difference formulation with numerical gradient
  5. ConcBuildFrm_form2_JACon_LM.mat   - Case 2(a): applying Levenberg-Marquardt algorithm on eigenvector difference formulation with analytical gradient 
  6. ConcBuildFrm_form2_JACon_TRR.mat  - Case 2(b): applying trust-region-reflective algorithm on eigenvector difference formulation with analytical gradient 
  
"SAP_Model"
* This is a folder including the sap2000 model of the concrete building frame

"ConcreteBuildingFrame.m" 
* This is the main function to perform FE model updating of concrete building frame from 100 start points;
* Before run the code, please add "_Shared" folder into MATLAB path. 
* The code first loads the influence matrices of updating variables, the mass matrix and the measured DOF from "ConcreBuildingFrame.mat". Then the code calls the functions in the '"_Shared" to perform FE model updating

"LoadStructure.m"
* This is a function called by "ConcreteBuildingFrame.m" to load the influence matrices of updating variables, mass matrix and build the nominal as well as actual stiffness matrices

"ConcreBuildingFrame.mat" 
* This is a file including the influence matrices of updating variables, mass matrix, and measured DOF

