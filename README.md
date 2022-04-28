# trig_approx

This project approximates sin() using polynomials. It adopts the following constraints:

 - The derivative at 0 must be exactly 1
 - The derivative at pi/2 and -pi/2 must be exactly 0
 - The value at pi/2 and -pi/2 must be exactly 1 and -1, respectively.

Within these constraints, this program will search for polynomial coefficients that result in the smallest maximum error.

Inspired by [this project](https://gist.github.com/publik-void/067f7f2fef32dbe5c27d6e215f824c91) which applies the same "searching for smallest maximum error" concept, but without the constraints. As a result, the polynomials output by this program have a larger maximum error, but are more accurate in applications where the values and derivatives at extremes are important.

## Example Output

```
Degree 5 polynomial:
0.0074030612083868476x^5 + -0.16553878047471235x^3 + 1x
Max absolute error: 0.00039453431471303535
Max relative error: 0.00047212644532912407
Degree 7 polynomial:
Minimized absolute error: [max error: 0.000002074751765612781, -0.0001824394856250362x^7 + 0.008303363983515513x^5 + -0.1666494845036777x^3 + 1x]
Minimized relative error: [max error: 0.000002647680715894829, -0.00018276495338425956x^7 + 0.008304970102529945x^5 + -0.1666514659735894x^3 + 1x]
Degree 9 polynomial:
Minimized absolute error: [max error: 0.000000009030746150884283, 0.000002581565088764897x^9 + -0.00019797584650388015x^7 + 0.008332882618479833x^5 + -0.1666665119037877x^3 + 1x]
Minimized relative error: [max error: 0.000000011973297020517748, 0.0000025862706501651165x^9 + -0.00019800456241893334x^7 + 0.008332938382504606x^5 + -0.16666654535711545x^3 + 1x]
Degree 11 polynomial:
Minimized absolute error: [max error: 0.000000000030410229889810125, -0.00000002374615310900146x^11 + 0.000002751645455079729x^9 + -0.00019840672469706794x^7 + 0.008333329365468965x^5 + -0.16666666574820746x^3 + 1x]
Minimized relative error: [max error: 0.00000000004138178688606331, -0.00000002378379639286943x^11 + 0.000002751919536787619x^9 + -0.0001984074347988978x^7 + 0.008333330125662047x^5 + -0.16666666602270294x^3 + 1x]
Degree 13 polynomial:
Minimized absolute error: [max error: 0.000009199689431804714, 0.00004642420816971904x^13 + -0.0003870840543244714x^11 + 0.001224243813519676x^9 + -0.001991878023959648x^7 + 0.009527929928152047x^5 + -0.16694376085913895x^3 + 1x]
Minimized relative error: [max error: 0.000013219676566533778, 0.000049091568631976115x^13 + -0.0003973068334962563x^11 + 0.0012069221863960503x^9 + -0.0018681684699810976x^7 + 0.009353744929261113x^5 + -0.16687196402433616x^3 + 1x]
```
