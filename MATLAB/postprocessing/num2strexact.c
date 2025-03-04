/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    num2strexact
 * Filename:    num2strexact.c
 * Programmer:  James Tursa
 * Version:     2.0
 * Date:        June 10, 2020
 * Copyright:   (c) 2009-2020 by James Tursa, All Rights Reserved
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are 
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in 
 *       the documentation and/or other materials provided with the distribution
 *      
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  Building:
 *
 *  num2strexact requires that a mex routine be built (one time only). This
 *  process is typically self-building the first time you call the function
 *  as long as you have the files num2strexact.m and num2strexact.c in the same
 *  directory somewhere on the MATLAB path. If you need to manually build
 *  the mex function, here are the commands:
 *
 *  >> mex -setup
 *    (then follow instructions to select a C / C++ compiler of your choice)
 *  >> mex num2strexact.c
 *
 * num2strexact is a mex function that converts a double or single input to
 * the exact decimal string. The conversion is done with hundreds of digits
 * of precision to maintain the exact conversion. The conversion uses the
 * exact decimal value of each bit of the IEEE double precision floating
 * point format along with the exact application of 2^exponent. Inf and NaN
 * bit patterns are recognized, and denormalized numbers are handled also.
 *
 * Don't confuse the exact conversion with significance! Double numbers will
 * only be significant to about 15 decimal digits, and single numbers will
 * only be significant to about 7 decimal digits. For example,
 *
 * >> format hex
 * >> 1.2
 * ans =
 *     3ff3333333333333
 * >> num2strexact(1.2)
 * ans =
 *     1.1999999999999999555910790149937383830547332763671875
 *
 * >> 1.2 + eps(1.2)
 * ans =
 *     3ff3333333333334   <-- one bit different from 1.2
 * num2strexact(1.2 + eps(1.2))
 * ans =
 *     1.20000000000000017763568394002504646778106689453125
 *
 * >> num2strexact(eps(1.2))
 * ans =
 *     2.220446049250313080847263336181640625e-16
 *
 * You can see that 1.2 is not represented exactly in IEEE double format.
 * The difference shows up in the 18th digit for this example. Then note
 * that the very next number in the IEEE double format model is about 2e-16
 * bigger. The exact conversions are shown for each number, but they are
 * not significant beyond the 16th digit shown. There are no numbers in
 * between these two numbers that can be represented in IEEE double format.
 *
 * Syntax:
 *
 *   Y = num2strexact(X [,'fixed' or 'float'])
 *   [Y1 Y2] = num2strexact(X1,X2 [,'fixed' or 'float'])
 *   [Y1 Y2 Y3] = num2strexact(X1,X2,X3 [,'fixed' or 'float'])
 *       :           :
 *      etc         etc
 *
 * X = double or single or half
 *
 * NOTE: The half type can be input in one of two ways:
 *       1) A uint16 class variable containing the half bit patterns
 *       2) A half class variable. num2strexact must have access to the
 *          underlying bit patterns, so if you input a half class variable,
 *          then this will first be converted to the equivalent integer bit
 *          patterns with the storedInteger function (a temporary deep copy).
 *
 * The number of inputs must match the number of outputs, except in the
 * special case of 1 input and 0 outputs where the result will simply be
 * put into ans. If the input is a scalar, the output will be a char string.
 * If the input is any other size array, the result will be a cell array of
 * the same dimensions as the input, with each cell containing the char
 * string conversion of the corresponding element.
 *
 * The optional 'fixed' argument forces the result to be fixed point.
 *
 * The optional 'float' argument forces the result to be floating point. All
 * numbers will be printed with exponents, even if the exponent is 0.
 *
 * The default will be fixed or floating point depending on the size of the
 * decimal exponent.  Numbers with -1 <= (decimal exponent) <= 2 will be
 * displayed as fixed point, else the number will be displayed as floating
 * point.
 *
 * All NaN bit patterns, regardless of the sign bit and payload, will be
 * returned simply as the string 'NaN'.
 *
 * Infinity bit patterns will be returned as 'Inf' or '-Inf' as appropriate.
 *
 * Revision History:
 *   2009-Aug-09 , Initial Release
 *   2020-Jun-10 , Added optional 'fixed' or 'float' input
 *                 Added half support (input as half or uint16 type)
 *
 ****************************************************************************/

#include "mex.h"
#include <string.h>
#include <stddef.h>

#include "matlab_version.h"
#if MATLAB_VERSION <= 0x2006a
#define mwSize int
#endif

// Macros used for printing the digits, etc.

#define N  400          // Size of integer array to hold the decimal digits
#define P  100          // Index location of 1's digit, next to decimal point
#define E    4          // Number of decimal digits to store in each integer
#define F "%4d"         // Format to print each integer with
#define EN1 (E*N-1)     // Output character string length
#define MAXVAL 10000uL  // Maximum integer value based on E
#define MAXHALF 5000uL  // Half of MAXVAL

// Need various integer sizes

#define  INT16_TYPE          short
#define UINT16_TYPE unsigned short
#define  INT32_TYPE          int
#define UINT32_TYPE unsigned int
#define  INT64_TYPE          long long

// Prototypes

mxArray *mxArray2digits(const mxArray *);
mxArray *double2digits(void *vdouble);
mxArray *single2digits(void *vsingle);
mxArray *half2digits(void *vhalf);
void halfp2singles(void *target, void *source, ptrdiff_t numel);

// Global variables

int fixed, floating;
int check = 1;

// Gateway function --------------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i; // Loop index
    double d = 1.0;  // Used to check IEEE double bit pattern
    char *cfixed;

    if( check ) {
//\
// Check to ensure that integer sizes are as expected
///
    
        if( sizeof(INT64_TYPE) != 8 || sizeof(UINT32_TYPE) != 4 || sizeof(UINT16_TYPE) != 2 ) {
            mexErrMsgTxt("Sizes of integer types are not as expected");
        }
    
//\
// Check to ensure that floating point bit pattern is IEEE double
///
    
        if( *((INT64_TYPE *)&d) != 0x3ff0000000000000LL )  // Using unspecified behaviour, type punning
            mexErrMsgTxt("Floating point bit pattern is not IEEE double");

        check = 0;
    }
    
/* Check for fixed point output */    
    
    fixed = 0;
    floating = 0;
    if( nrhs && mxIsChar(prhs[nrhs-1]) ) {
        cfixed = mxArrayToString(prhs[nrhs-1]);
        if( strcmp(cfixed,"fixed") == 0 || strcmp(cfixed,"fix") == 0 ) {
            fixed = 1;
            nrhs--;
        }
        if( strcmp(cfixed,"float") == 0 || strcmp(cfixed,"floating") == 0 ) {
            floating = 1;
            nrhs--;
        }
        mxFree(cfixed);
    }
    
//\
// Require that the number of input and output arguments match, except for
// the special case of 1 input argument and 0 output arguments.
///
    
    if( (nrhs != nlhs) && ((nrhs != 1) || (nlhs != 0)) )
        mexErrMsgTxt("Mismatched number of input/output arguments\n");

//\
// Nothing to convert
///
    
    if( nrhs == 0 )
        return;

//\
// Check each input argument for correct type
///
    
    for( i=0; i<nrhs; i++ ) {
        if( mxIsComplex(prhs[i]) || mxIsSparse(prhs[i]) || 
            !( mxIsDouble(prhs[i]) || mxIsSingle(prhs[i]) || mxIsUint16(prhs[i]) || mxIsClass(prhs[i],"half") ) )
            mexErrMsgTxt("Input arguments must be full real double or single or half or uint16");
    }
    
//\
// Convert each input argument
///
    
    for( i=0; i<nrhs; i++ )
        plhs[i] = mxArray2digits(prhs[i]);
    
}

//----------------------------------------------------------------------------

mxArray *mxArray2digits(const mxArray *mx)
{
    mwSize i, numel, ndim;
    const mwSize *dims;
    mxArray *rx, *hx;
    double *dp;
    float *fp;
    UINT16_TYPE *hp;
    int ishalf;
    
    if( ishalf = mxIsClass(mx,"half") ) {
        if( mexCallMATLAB( 1, &hx, 1, &mx, "storedInteger" ) ) {
            mexErrMsgTxt("Cannot get stored integer for half variable");
        }
        mx = hx;
    }
    
    numel = mxGetNumberOfElements(mx);
    
    if( numel == 0 ) {  // Empty input, return empty char string
        return mxCreateString("");
        
    } else if( numel == 1 ) {  // Scalar input, return one char string
        if( mxIsDouble(mx) ) {
            return double2digits(mxGetData(mx));
        } else if( mxIsSingle(mx) ) {
            return single2digits(mxGetData(mx));
        } else {
            hx = half2digits(mxGetData(mx));
            if( ishalf ) mxDestroyArray(mx);
            return hx;
        }
        
    } else {  // Non-scalar array input, return cell array same dimensions
        ndim = mxGetNumberOfDimensions(mx);
        dims = mxGetDimensions(mx);
        rx = mxCreateCellArray(ndim, dims);
        if( mxIsDouble(mx) ) {
            dp = (double *) mxGetData(mx);
            for( i=0; i<numel; i++ ) {
                mxSetCell(rx, i, double2digits((void *)dp++));
            }
        } else if( mxIsSingle(mx) ) {
            fp = (float *) mxGetData(mx);
            for( i=0; i<numel; i++ ) {
                mxSetCell(rx, i, single2digits((void *)fp++));
            }
        } else {
            hp = (UINT16_TYPE *) mxGetData(mx);
            for( i=0; i<numel; i++ ) {
                mxSetCell(rx, i, half2digits((void *)hp++));
            }
            if( ishalf ) mxDestroyArray(mx);
        }
        return rx;
    }
}

//----------------------------------------------------------------------------

mxArray *double2digits(void *vdouble)
{
    UINT32_TYPE d[N] = {0uL};  // All values set to 0
    UINT32_TYPE t[N] = {0uL};  // All values set to 0
    UINT32_TYPE carry, dnew, bits, carrybits;
    UINT32_TYPE carryvalue[] = {    0uL,  625uL, 1250uL, 1875uL,
                            2500uL, 3125uL, 3750uL, 4375uL,
                            5000uL, 5625uL, 6250uL, 6875uL,
                            7500uL, 8125uL, 8750uL, 9375uL };
    INT64_TYPE usign, uexpo, umant, umbit;
    INT64_TYPE *udouble;
    int texpo, dexpo, k, i, j, jfirst, jlast, c, cfirst, clast, check0;
    char out[EN1];  // Character string output
    
//\
// Interpret the 8-byte double as an 8-byte integer. I should be using an
// unsigned 64-bit integer for this, but the lcc compiler that comes with
// MATLAB does not handle unsigned 64-bit integers properly, so I am forced
// to use signed 64-bit integers instead.
///
    
    udouble = (INT64_TYPE *) vdouble;
    
//\
// Is the input 0.0?
///
    
    if( *udouble == 0LL ) {
        if( floating ) {
            return mxCreateString("0e0");
        } else {
            return mxCreateString("0");
        }
    }
    if( *udouble == 0x8000000000000000LL ) {
        if( floating ) {
            return mxCreateString("-0e0");
        } else {
            return mxCreateString("-0");
        }
    }
    
//\
// Pick off the sign bit, exponent bits, and mantissa bits
///

    usign = 0x8000000000000000LL & *udouble;
    uexpo = 0x7FF0000000000000LL & *udouble;
    umant = 0x000FFFFFFFFFFFFFLL & *udouble;
    
//\
// Check for special bit patterns
///
    
    if( uexpo == 0uLL ) {  // Denormalized Number
        d[P] = 0uL;
        t[P] = 1uL;
        jfirst = P;
        jlast  = P;
    } else if( uexpo == 0x7FF0000000000000LL ) {
        if( umant ) {
            return mxCreateString("NaN");
        } else if( usign ) {
            return mxCreateString("-Inf");
        } else {
            return mxCreateString("Inf");
        }
    } else {  // Normalized Number
        d[P] = 1uL;
        t[P+1] = MAXHALF;
        jfirst = P;
        jlast  = P + 1;
    }
    
//\
// Calculate actual power of 2 exponent
///
    
    texpo = ((int)(uexpo >> 52)) - 1023;
    
//\
// Convert mantissa bits into decimal integer values one at a time.
// The d array is the accumulated mantissa decimal digits. The t array
// is the exact power of 2 for the particular mantissa bit being converted
// prior to the application of the exponent. At the end of this section d
// will be 1.xxxxxxxxxx...
///
    
    umbit = 0x0010000000000000LL;
    while( umbit >>= 1 ) {  // For each mantissa bit
        if( umant & umbit ) {  // If mantissa bit is set
            carry = 0uL;
            for( i=jlast; i>=P; i-- ) {  // Loop to add d = d + t
                d[i] += t[i] + carry;
                dnew = d[i] % MAXVAL;
                if( dnew != d[i] ) {
                    carry = d[i] / MAXVAL;
                    d[i] = dnew;
                } else {
                    carry = 0uL;
                }
            }
        }
        carry = 0uL;
        for( i=P; i<=jlast; i++ ) {  // Loop to divide t = t / 2
            if( (t[i] += carry) & 1uL ) {
                carry = MAXVAL;
            } else {
                carry = 0uL;
            }
            t[i] >>= 1;  // t[i] = t[i] / 2
        }
        if( carry ) {
            t[++jlast] = MAXHALF;
        }
    }
    
//\
// Clean up ending zeros
///
    
    while( d[jlast] == 0uL ) {
        jlast--;
    }
    
//\
// Apply the positive exponent
///
    
    k = 16;  // Start with 2^16
    while( texpo > 0 ) {
        if( texpo >= k ) {
            check0 = 1;
            carry = 0uL;
            for( i=jlast; i>=jfirst; i-- ) {  // Loop to multiply d = d * 2^k
                d[i] = (d[i] << k) + carry;
                dnew = d[i] % MAXVAL;
                if( dnew != d[i] ) {
                    carry = d[i] / MAXVAL;
                    d[i] = dnew;
                } else {
                    carry = 0uL;
                }
                if( check0 && dnew == 0uL ) {
                    --jlast;
                } else {
                    check0 = 0;
                }
            }
            while( carry ) {  // Loop to add left-over carry value to front end
                d[--jfirst] = carry;
                dnew = d[jfirst] % MAXVAL;
                if( dnew != d[jfirst] ) {
                    carry = d[jfirst] / MAXVAL;
                    d[jfirst] = dnew;
                } else {
                    break;
                }
            }
            texpo -= k;  // Decrease exponent texpo = texpo - k
        } else {
            k >>= 1;  // Divide k = k / 2
        }
    }
    
//\
// Apply the negative exponent
///
    
    k = 4;  // Start with 2^4
    bits = 0xFuL;  // Check the 4 trailing bits
    while( texpo < 0 ) {
        if( texpo <= -k ) {
            check0 = 1;
            carry = 0uL;
            for( i=jfirst; i<=jlast; i++ ) {  // Loop to divide d = d / 2^k
                carrybits = d[i] & bits;
                d[i] >>= k;  // Divide d[i] = d[i] / 2^k
                d[i] += carry; // Add in the carry value from the previous pass
                carry = carryvalue[carrybits << (4 - k)];
                if( check0 && d[i] == 0uL ) {
                    ++jfirst;
                } else {
                    check0 = 0;
                }
            }
            if( carry ) {  // Add left-over carry value to back end
                d[++jlast] = carry;
            }
            texpo += k;  // Increase exponent texpo = texpo + k
        } else {
            --k;  // Decrement k = k - 1
            bits >>= 1;  // Shift bits right by 1 (i.e., check 1 fewer trailing bit)
        }
    }
    
//\
// Print digits to string
///
    
    for( i=jfirst; i<=jlast; i++ ) {
        sprintf(out+i*E, F, d[i]);
    }
    
//\
// Indexes of first and last digits to print
///
    
    cfirst = jfirst * E;
    clast  = jlast  * E + E;
    
//\
// Replace blanks with zeros
///
    
    for( c=cfirst; c<clast; c++ ) {
        if( out[c] == ' ' ) out[c] = '0';
    }
    
//\
// Clean up trailing zeros
///
    
    out[clast--] = '\0';
    while( clast > cfirst && out[clast] == '0' ) {
        out[clast--] = '\0';
    }
    
//\
// Fix up exponent depending on how big the first integer is
///
    
    dexpo = (P - jfirst) * E;
    if( d[jfirst] < 10 ) {
        cfirst += 3;
    } else if( d[jfirst] < 100 ) {
        cfirst += 2;
        dexpo++;
    } else if( d[jfirst] < 1000 ) {
        cfirst++;
        dexpo += 2;
    } else {
        dexpo += 3;
    }
    
//\
// Place the decimal point, sign for fixed point
///

    if( fixed || (dexpo >= -1 && dexpo <= 2 && !floating) ) {
        
        if( dexpo == 0) {
            out[cfirst-1] = out[cfirst];
            out[cfirst--] = '.';
            
        } else if( dexpo > 0 ) {
            i = cfirst--;
            dexpo++;
            while( dexpo-- ) {
                if( i <= clast ) {
                    out[i-1] = out[i];
                } else {
                    out[i-1] = '0';
                }
                i++;
            }
            out[i-1] = '.';
            if( i > clast ) {
                clast = i - 1;
                out[i] = '\0';
            }
            
        } else {
            i = cfirst - 1;
            dexpo++;
            while( dexpo++ ) {
                out[i--] = '0';
            }
            out[i--] = '.';
            out[i] = '0';
            cfirst = i;
            
        }
        if( usign ) {
            out[--cfirst] = '-';
        }
        if( out[clast] == '.' ) {
            out[clast--] = '\0';
        }
        
    } else {
    
    
//\
// Place the decimal point, sign, and exponent for floating point
///
    
        if( dexpo == -1 && !floating ) {
            out[--cfirst] = '.';
            out[--cfirst] = '0';
            dexpo = 0;
        } else {
            out[cfirst-1] = out[cfirst];
            out[cfirst--] = '.';
        }
        if( usign ) {
            out[--cfirst] = '-';
        }
        if( out[clast] == '.' ) {
            out[clast--] = '\0';
        }
        if( dexpo != 0 || floating ) {
            out[++clast] = 'e';
            if( dexpo < 0 ) {
                out[++clast] = '-';
                dexpo = -dexpo;
            }
            clast++;
            if( dexpo < 10 ) {
                sprintf(out+clast, "%1d", dexpo);
            } else if( dexpo < 100 ) {
                sprintf(out+clast, "%2d", dexpo);
            } else {
                sprintf(out+clast, "%3d", dexpo);
            }
        }
        
    }
    
    return mxCreateString(out+cfirst);
}

//----------------------------------------------------------------------------

mxArray *single2digits(void *vsingle)
{
    double d;
    
    d = (double)(*((float *)vsingle));  // No custom single code. Just convert to
    return double2digits((void *)&d);   // double and call the double code.
}

//----------------------------------------------------------------------------

mxArray *half2digits(void *vhalf)
{
    double d;
    float f;
    
    halfp2singles(&f, vhalf, 1);       // No custom half code. Convert to single
    d = f;                             // then to double
    return double2digits((void *)&d);  // then call double code.
}

//-----------------------------------------------------------------------------

// This routine borrowed from halfprecision.c code
// https://www.mathworks.com/matlabcentral/fileexchange/23173-ieee-754r-half-precision-floating-point-converter

void halfp2singles(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) target; // Type pun output as an unsigned 32-bit int
    UINT16_TYPE h, hs, he, hm;
    UINT32_TYPE xs, xe, xm;
    INT32_TYPE xes;
    ptrdiff_t i;
    int e;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
    
    for( i=0; i<numel; i++ ) {
        h = hp[i];
        if( (h & 0x7FFFu) == 0 ) {  // Signed zero
            xp[i] = ((UINT32_TYPE) h) << 16;  // Return the signed zero
        } else { // Not zero
            hs = h & 0x8000u;  // Pick off sign bit
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( he == 0 ) {  // Denormal will convert to normalized
                e = -1; // The following loop figures out how much extra to adjust the exponent
                do {
                    e++;
                    hm <<= 1;
                } while( (hm & 0x0400u) == 0 ); // Shift until leading bit overflows into exponent bit
                xs = ((UINT32_TYPE) hs) << 16; // Sign bit
                xes = ((INT32_TYPE) (he >> 10)) - 15 + 127 - e; // Exponent unbias the halfp, then bias the single
                xe = (UINT32_TYPE) (xes << 23); // Exponent
                xm = ((UINT32_TYPE) (hm & 0x03FFu)) << 13; // Mantissa
                xp[i] = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            } else if( he == 0x7C00u ) {  // Inf or NaN (all the exponent bits are set)
                xp[i] = (((UINT32_TYPE) hs) << 16) | ((UINT32_TYPE) 0x7F800000u) | (((UINT32_TYPE) hm) << 13); // Signed Inf or NaN
            } else { // Normalized number
                xs = ((UINT32_TYPE) hs) << 16; // Sign bit
                xes = ((INT32_TYPE) (he >> 10)) - 15 + 127; // Exponent unbias the halfp, then bias the single
                xe = (UINT32_TYPE) (xes << 23); // Exponent
                xm = ((UINT32_TYPE) hm) << 13; // Mantissa
                xp[i] = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            }
        }
    }
}
