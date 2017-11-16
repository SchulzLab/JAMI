

In adaptive partitioning we split a cube of $X, Y, Z$ recursively in eight smaller cubes until the small cubes are balanced. Subsequently, we sum up the CMI values of each sub-cube. Before starting with the recursive splitting, we can split off a sub-cube containing all the samples for which there are zeros on an axis and compute the CMI value for this sub-cube as a special case. The CMI value of cube $C = [X, Y, Z]$ is

$CMI_C = I(X; Y | Z)$.

We initialize with $CMI_C L:= 0$.

We consider three cases in succession:

### Case 1

In case $X \tilde= (0)^{m} *^{k}$ , i.e. m zeros on the X axis followed by k non-zero values (0s could also be missing values in the data (NAs)), we split of the sub-cube $C_{X_0} = [X_1^m = (0)^{\{m\}}, Y_1^m, Z_1^m]$

Since $X_1^m$ does not contribute any information, it follows that $CMI(C_{X_0}) = I(X_1^m;Y_1^m|Z_1^m) = 0$

We update
$CMI_C := CMI_C + 0$

In other words, in the sub-cube where we have only 0s on the $X$ axis, the CMI is zero, adding nothing to the overall CMI. The remaining cube, e.g.
$C := C \setminus C_{X_0}$ is considered further.

### Case 2
Is analogue to case 1 for swapping X and Y.

$CMI_C := CMI_C + 0$
$C := C \setminus C_{Y_0}$ 

### Case 3

If $Z \tilde= (0)^m*^k$

We split of the sub-cube $C_{Z_0} = [X_1^m, Y_1^m, Z_1^m  = {(0)^m}]$

$CMI(C_{Z_0}) = I(X_1^m;Y_1^m|Z_1^m) = I(X_1^m; Y_1^m)$

In other words, since the variable $Z$ does not yield any information, the CMI is the same as the mutual information of just X and Y. The difference of CMI and mutual information, otherwise known as interaction information, is zero.

$C := C \setminus C_{Z_0}$

The remaining cube $C$ is handled by the adaptive partitioning algorithm and is proven to be correct as long as we don't have duplicated expression values, which is unlikely in expression data except for the case of zeros that we just covered. 


