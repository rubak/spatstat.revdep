# CalSim 0.5.1
1. introduced multinomial tests

2. removed input vectors p_b,p_n,p_a for calibration_simplex (use p1,p2,p3 instead!)

3. fixed some bugs (mainly concerning very small samples):
 - samples can all be within a single bin
 - not all outcomes must occur
 - plot works for n = 1

4. reduced size of sample data from 50,000 rows to 10,000 rows

# CalSim 0.3.2
- added sample data and code examples illustrating the calibration simplex
- removed code example, which was using data from the scoring package
- changed input vector names for calibration_simplex (p_b -> p1, p_n -> p2, p_a -> p3). The old variable names are still supported, however they will trigger a warning and may be removed in future versions of this package.
- added (optional) parameter category_labels to the plot function and changed default labels to 1,2,3 instead of b,n,a.
- added README.md and NEWS.md

# CalSim 0.2.2
- initial package
