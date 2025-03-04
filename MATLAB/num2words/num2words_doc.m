%% NUM2WORDS Examples
% The function <https://www.mathworks.com/matlabcentral/fileexchange/47221
% |NUM2WORDS|> converts a numeric scalar into a string with the number value
% given in English words, e.g. |1024| -> |'one thousand and twenty-four'|.
% Optional arguments control the handling of the numeric input, as well as
% many string formatting and dialect options.
% The options are explained in this document, together with examples.
%
% The number format is based on <http://www.blackwasp.co.uk/NumberToWords.aspx>
%
% In the |NUM2WORDS| documentation the word _string_ refers to a 1xN
% character vector, and not to the MATLAB string class.
%% Basic Usage
% For integer values |NUM2WORDS| can be called without any options:
num2words(0)
num2words(+1024)
num2words(-1024)
%% Decimal Digits
% |NUM2WORDS| rounds the input number value to the digit specified by the
% either the |order| or the |sigfig| option:
%
% * |sigfig| rounds the value to the requested number of significant figures,
% * |order| rounds the value to the specified order of magnitude.
%
% The default |order| is zero, which means that any fractional part is
% rounded to give an integer. Fractional digits can be included by selecting
% an appropriate value for either the |order| or |sigfig| options:
num2words(1.234) % default 'order' is zero
num2words(1.234, 'order',-2)
num2words(1.234, 'sigfig',3)
%% Trailing Zeros
% By default trailing fractional zeros are not included in the string.
% Setting the |trz| option to |true| will keep
% the trailing zeros up to the requested |order| or |sigfig|:
num2words(1, 'sigfig',3, 'trz',false) % default
num2words(1, 'sigfig',3, 'trz',true)
%% Floating-Point Precision
% |NUM2WORDS| has internal limits on the significant figures for all
% <http://www.mathworks.com/help/matlab/matlab_prog/floating-point-numbers.html
% floating-point numbers>: 15 for double, and 6 for single. These limits
% ensure that the least-unexpected output is returned. As |NUM2WORDS| is
% based around |SPRINTF|, an |SPRINTF| example shows why this is required:
sprintf('%#.15g',1e23) % fifteen significant figures
sprintf('%#.16g',1e23) % sixteen significant figures
%% Integer Precision
% In contrast, all <http://www.mathworks.com/help/matlab/matlab_prog/integers.html
% integer class numbers> are parsed at their full precision (even |UINT64|,
% and positive and negative |INT64| values):
num2words(intmax('uint64'))
num2words(intmin('int64'), 'pos',true)
num2words(intmax('int64'), 'pos',true)
%% Type: Decimal
% The default type: it can be used for all integer and decimal fraction
% values. All of the above examples use this number type.
%% Type: Ordinal
% The last number word is changed to have an ordinal ending:
num2words(1, 'type','ordinal')
num2words(12, 'type','ordinal')
num2words(123, 'type','ordinal')
%% Type: Highest
% This uses the highest magnitude word together with a significand that
% includes fractional digits when required. In most cases it would be
% useful to also specify the number of significant figures:
num2words(1234567.89, 'type','highest')
num2words(1234567.89, 'type','highest', 'sigfig',2)
%% Type: Money
% Treats the number as being a value of currency, and returns a string that
% contains the currency unit and subunit names. The currency subunit is
% 1/100 of the unit, which suits almost all currencies in the world today.
% The currency unit and subunit names can be supplied as string arguments:
% the string format also selects whether the name has a regular, irregular
% or invariant plural form (see the |NUM2WORDS| help for details).
%
% Note that the default |order| is changed to -2, to match the subunit.
num2words(23.5, 'type','money')
num2words(23.5, 'type','money', 'unit','Pound|', 'subunit','Penny|Pence')
num2words(101,  'type','money', 'unit','Dalmatian|', 'case','title')
num2words(1001, 'type','money', 'unit','Night|', 'case','title')
%% Type: Cheque
% Similar to |money|, but follows the style used in many countries: if no
% trailing subunits then the units are followed by the word |'only'|, and
% leading units are always indicated, even if they are zero (preventing fraud):
num2words(5.0,'type','cheque')
num2words(0.5,'type','cheque')
%% String: Capitalisation
% The default returns a lower-case string. It is also possible to select
% upper-case, title-case (all words except |'and'| have an initial capital
% letter) and sentence-case (only the first word has an initial capital):
num2words(-1023,'case','lower') % default
num2words(-1023,'case','upper')
num2words(-1023,'case','title')
num2words(-1023,'case','sentence')
%% String: Grammar
% There are some string formatting options that allow for different English
% dialects and style guidelines. The option defaults can be changed inside
% the mfile, if you wish to make another English dialect the default.
num2words(1234,'and',true,'comma',true, 'hyphen',true) % default
num2words(1234,'and',false,'comma',false, 'hyphen',false)
%% String: 'Positive' Prefix
% |NUM2WORDS| always prepends |'negative'| for negative number values
% (including negative zero). The |'positive'| prefix is optional:
num2words(1, 'pos',false) % default
num2words(1, 'pos',true)
%% Number Scales
% |NUM2WORDS| supports a number of common and not-so-common
% <https://en.wikipedia.org/wiki/Names_of_large_numbers number scales>:
%
% * <https://en.wikipedia.org/wiki/Short_and_long_scale |short| and
% |long| scales> are explained on many websites. Most contemporary English
% dialects use the |short| scale (and this is the |NUM2WORDS| default).
% * <https://en.wikipedia.org/wiki/Indian_numbering_system
% |indian| number system> is commonly used in South Asia (the Indian
% subcontinent). Defaults back to |short| for values greater than 1e21.
% * <https://en.wikipedia.org/wiki/Jacques_Pelletier_du_Mans |peletier|
% scale> is used in many non-English speaking European countries.
% * <http://www.unc.edu/~rowlett/units/large.html |rowlett| scale> was
% designed to avoid the ambiguity of the |short| and |long| scales.
% * <https://sites.google.com/site/pointlesslargenumberstuff/home/1/knuthyllions
% |knuth| scale> (aka _-yllion_) uses a logarithmic naming system to use
% very few names to cover a very wide range of values.
num2words(1e9, 'scale','short') % default
num2words(1e9, 'scale','long')
num2words(1e9, 'scale','indian')
num2words(1e9, 'scale','peletier')
num2words(1e9, 'scale','rowlett')
num2words(1e9, 'scale','knuth')
%% Faster Conversion: |NUM2WORDSQ|
% The bonus function |NUM2WORDSQ| does not support any options and runs
% around twice as fast as |NUM2WORDS|:
num2wordsq(1024)
%% Reverse Conversion: |WORDS2NUM|
% The function <https://www.mathworks.com/matlabcentral/fileexchange/52925
% |WORDS2NUM|> converts a number string into a numeric value:
words2num('one thousand and twenty-four')