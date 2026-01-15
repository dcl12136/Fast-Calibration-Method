Title: Source Code for "Fast Calibration Method of Axial Distance Error in Ptychography Based on Linear Model"

Description:
This package contains the MATLAB implementation of the proposed "Linear Model-based Binary Search Autofocusing Strategy" (Algorithm 1 in the manuscript).

Files Included:
1. LinearModel.m
   - The core function for Linear Model Initialization (Line 3 in Algorithm 1).
   - Analytically reconstructs the initial complex object.

2. dPIE.m
   - The function for sharpness evaluation and axial distance update (Lines 6-13 in Algorithm 1).

3. Ptychography.m
   - Standard iterative update engine (ePIE) (Line 18 in Algorithm 1).

4. Propagate.m
   - Numerical implementation of angular spectrum propagation.

Usage Guide:
These functions are designed as modular components. Users can reconstruct the workflow of Algorithm 1 by integrating these functions into their main simulation loop as follows:

    % --- Pseudo-code for Integration ---
    
    % Step 1: Initialization
    % Use LinearModelePIE3 to get the initial object estimate
    [Obj_Init, ~] = LinearModelePIE3(N, Obj_Mov, delta_L0, Dect_Std);

    % Step 2: Autofocusing Loop
    % Feed Obj_Init into the iterative loop
    while not_converged
        % Update Object
        [Obj_Updated, ...] = Ptychography_1(..., Obj_Init, ...);
        
        % Update Distance (Sharpness check & Binary search)
        [new_dist, ...] = dPIE(Obj_Updated, current_dist, ...);
    end

Note:
This package provides the core algorithmic modules. Users are expected to provide their own main simulation scripts and ptychographic datasets.
