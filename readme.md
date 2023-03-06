<img src="assets/transmission_system_diagram.png">

# Design Optimization of CTR Transmission System

## Introduction
This repository contains code about how to design an optimal pyramid-shaped transmission system for multi-arm concentric-tube robots (CTR) . 

* Full paper PDF: [Manuscript submitted to TMECH.](assets/Manuscript.pdf)

* Authors: *Chao Zhang, Guangdu Cen, Xing Yang, et al.*


## File Tree

* assets (*PDF and figures, etc*)
* Transmission_System_Optimal_Design (*Matlab code*)
  * Pyramid-shaped_multi-arm
    * Case 1
      * CTR_params.m
      * determine_result.m
      * 
    * Case 2
  * Y-shaped_dual-arm

## Runding the Code
1. Add the path to Matlab commandd window

2. Run "optimal_design_framework.m" to get the Pareto front

<img src="assets/pareto_front.png">

3. Run "variables_varying.m" to see the variation of optimization variable **x** = [OD~s1~, OD~s2~, ..., OD~sm~, l~b~, \alpha, l~e~, \beta~1~, \beta~2~, ..., \beta~n-1~] ; 

<img src="assets/variable_variation.png">

4. Run "determine_result.m" to get the final optimal result






