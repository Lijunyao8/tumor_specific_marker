package Statistics::Outliers;

################################################################################
##                                                                            ##
##         I N T R O D U C T O R Y   P O D   D O C U M E N T A T I O N        ##
##                                                                            ##
################################################################################

=head1 NAME

Statistics::Outliers - identify outliers in a list of numerical data

=head1 SYNOPSIS

	use Statistics::Outliers ':all';
	use Statistics::Outliers 'tukey_iqr';
	use Statistics::Outliers 'peirce_criterion';
	use Statistics::Outliers 'median_absolute_dev';

=head1 DESCRIPTION

Data L<outliers|https://en.wikipedia.org/wiki/Outlier>
can arise from various error processes in experimental
work, including instrument or human error, contamination,
etc.
Such elements can pose enormous problems for statistical data
analysis.
Reliably identifying and eliminating outliers is
therefore a chronic problem for an enormously wide array of
applications in science, engineering, economics,
etc.
This package implements several well-known methods from
the mathematical statistics literature toward solving this
problem, although there are many others that are not
included.
Since there is really no formal mathematical specification
of what constitutes an outlier, there is likewise no
universal, optimal technique for identifying and removing
them.

=head1 SCOPE

This package has a limited scope of
methods.
First, it focuses strictly on outlier I<removal> and has no procedures for
clipping, aka L<"winsorizing"|https://en.wikipedia.org/wiki/Winsorizing>.
(This is where outlier values are I<replaced> by more reasonable values,
rather than being
removed.)
Second, there are several traditional methods that are omitted here for reasons
noted:

=over

=item *

Dixon's Outlier Test [Dixon50] is not included because it is
restricted to small data sets, generally taken to be no more than 25
elements [Sokal95], and is intended to find only a single
outlier.

=item *

Also not implemented here are casual methods, like the
L<3-sigma rule|https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule>,
e.g. where the Guassian
Z > 3.
They are trivial to construct using
L<existing packages|Statistics::Descriptive>.
More importantly, they are not particularly robust because they have
no provision to compensate for the effect of data outliers on the
calculated parameters of the model itself, i.e. the sample mean and standard
deviation.

=back

There are many other worthwhile methods that might
be incorporated into future versions of this
package.
The Wikipedia page on L<outliers|https://en.wikipedia.org/wiki/Outlier>
has a partial survey of some of
them.
[Rousseeuw93] also discusses several very promising
methods.

=head1 METHODS

The subtleties and esoterica of various statistical methods for identifying
and removing outliers require some familiarity with the mathematical statistics
literature.
We give rules-of-thumb for specific approaches, where they exist, but it may
be advisable to try several different ones to check the consistency of the
result.
For example, some methods depend to lesser or greater extents
on whether the data can be reasonably assumed to have come from a
(L<normally-distributed|http://en.wikipedia.org/wiki/Normal_distribution>)
population (testable using
L<Statistics::Normality>).
There are aspects of mathematical assumptions, the
subjectivity of the process, and the user's level of due
diligence that often render outlier removal a very tricky
proposition.
It is probably a good idea to examine the data
graphically.

Each function in this package accepts a reference
to an unordered list of floating-point values, e.g.

	$list = [2.5, 3.1, 9.5, 1.334, ...]

and additional options (in some cases) and returns an ordered list with outliers
removed.

=cut

#__STANDARD PERL PACKAGES
   use strict;
   use Carp;
   use Statistics::Descriptive;

#__SET UP EXPORTING
   use vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION );
   require Exporter;
   @ISA = qw(Exporter);
   @EXPORT = qw();
   @EXPORT_OK = qw/
      tukey_iqr
      peirce_criterion
      median_absolute_dev
   /;
   %EXPORT_TAGS = (all => [qw/
      tukey_iqr
      peirce_criterion
      median_absolute_dev
   /]);
   my $pkg = 'Statistics::Outliers';

#__VARIABLES
   our ($VERSION, $errtxt);
   $VERSION = '0.1';

################################################################################
##                                                                            ##
##                        P U B L I C   R O U T I N E S                       ##
##                                                                            ##
################################################################################

=head2 Peirce's criterion

Peirce's iterative criterion is one of the original
outlier methods [Peirce52,Gould55], but is nevertheless quite
sophisticated and still considered to be an extremely robust method
[Ross03].

	$list = peirce_criterion ($list);

but a 2nd argument can also be passed optionally to skip internal sorting (to
save CPU cycles) if the input data are already sorted lowest-to-highest.

	$list = peirce_criterion ($list, 1);

Called in a list context, the method returns the reference to the
trimmed list, as well as references to lists of low-end and high-end
values that were filtered:

	($list, $lo_rejects, $hi_rejects) = peirce_criterion ($list);

Because there is an iteration process for determining how many
elements are outliers, as well as an inner iteration for each of these
outer iterations to find the range ratio, this function requires
somewhat more computational effort than non-iterative methods,
which may be a consideration if there are many datasets that require
processing.

=cut

#  LOCAL PROGRAMMER NOTES:
#
#  (1) CERTIFICATION: THIS IMPLEMENTATION GETS THE EXACT SAME RESULT AS THE
#      EXAMPLE GIVEN IN [Ross03] (MCW, 8-DEC-2018)
#
#  (2) THIS DOES *NOT* GET BOTH OUTLIERS IN THE EXAMPLE SHOWN AT
#      eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers
#
#      my $list = [qw/1 2 3 3 4 4 4 5 5.5 6 6 6.5 7 7 7.5 8 9 12 52 90/];
#                                          .----------------------+  |
#           but misses this one (barely)--^    only finds this one---+
#
#  (3) ROUTINE USES AN INDICATOR LIST WHICH IS DETERMINED DURING THE COURSE
#      OF THE CALCULATION (0 = NOT OUTLIER, 1 = OUTLIER), TO INITIALIZE IT WE
#      USE THE REPETITION ("x") OPERATOR IN A LIST CONTEXT BY ENCLOSING IT IN
#      PARENTHESES, E.G.
#
#      @ones = (1) x 80;          # a list of 80 1's
#
#      SEE: http://perldoc.perl.org/perlop.html#Multiplicative-Operators
#
# hardwired for m = 1: Dardis paper pp 4 -- m is degrees of freedom
# and m>2 not of practical importance

sub peirce_criterion {
   my ($list, $no_sort) = @_;
   $no_sort = 0 unless $no_sort;

#__ARGUMENT MUST BE REFERENCE TO LIST OF NUMBERS
   my $error = _check_args_ ($list);
   croak $error if defined $error;

#__ORDER THE LIST SMALL TO LARGE UNLESS INPUT IS ALREADY SORTED
   @{$list} = sort _numerical_ @{$list} unless $no_sort;
   my $num_vals = scalar @{$list};

#__INITIALIZE CORRESPONDING LIST OF INDICATORS (0 = NOT OUTLIER, 1 = OUTLIER)
   my @indicators = (0) x $num_vals;

#__CALCULATE MEAN AND STANDARD DEVIATION
   my ($mean, $standev) = _mean_and_standev_ ($list);

#__ITERATER OVER SERIES OF DECREASING MAXIMUM DEVIATIONS
   my $j = 1;
   while ($j <= $num_vals) {

   #__CALCULATE GOULD'S X (AKA PEIRCE'S RATIO R)
      my $x = _goulds_x_ ($num_vals, $j);

   #__MAXIMUM ALLOWABLE DEVIATION
      my $max_dev = $standev * $x;

   #__CALCULATE ACTUAL DEVIATIONS AND COMPARE TO MAXIMUM ALLOWABLE DEVIATION
      my $num_outliers_on_this_round = 0;
      for (my $i = 0; $i <= $#{$list}; $i++) {

      #__SKIP ANY ALREADY-KNOWN OUTLIERS
         next if $indicators[$i];

      #__DEVIATION OF CURRENT ELEMENT
         my $dev = abs ($mean - $list->[$i]);

      #__MARK AS AN OUTLIER IF IT EXCEEDS THE MAXIMUM
         if ($dev > $max_dev) {
            $indicators[$i] = 1;
            $num_outliers_on_this_round++;
         }
      }

   #__ADVANCE THE GOULD COUNTER BY THE NUMBER OF OUTLIERS ELIMINATED THIS ROUND
      $j += $num_outliers_on_this_round;

   #__END THE IDENTIFICATION PROCESS IF NO ELEMENTS ELIMINATED ON THIS ROUND
      last unless $num_outliers_on_this_round;

   #__SANITY CHECK: IS THIS A STRANGE CASE WHERE MUCH OF THE DATA IS FLAGGED
      croak "funny business: too much data eliminated" if $j > $num_vals / 2;
   }

#__ELIMINATE OUTLIERS FROM THE DATA ACCORDING TO ACCOMPANYING INDICATOR ARRAY
   my ($num_low_outs, $num_high_outs) = _read_indicators_ (@indicators);
   my ($low_end_rejects, $high_end_rejects) = ([], []);
   @{$low_end_rejects} = splice (@{$list}, 0, $num_low_outs) if $num_low_outs;
   @{$high_end_rejects} = splice (@{$list}, -$num_high_outs) if $num_high_outs;

#__RETURN CLEANED LIST
   if (wantarray) {
      return ($list, $low_end_rejects, $high_end_rejects);
   } else {
      return $list;
   }
}

# ------------------------------------------------------------------------------

=head2 Tukey's Inter-Quartile Method

This non-parametric approach, frequently referred to as the interquartile
(IQR) method, is based on identifiying low/hi outliers based on distance
of values below/above (respectively) the 1st and 3rd quartiles
[Tukey77,McGill78].
By default this distance is 1.5 times the inter-quartile distance:

	$list = tukey_iqr ($list);

but the method takes an optional second argument to change the
coefficient:

	$list = tukey_iqr ($list, 2.5);

Note that the method will L<carp|Carp> unless
the coefficient is within the range of 1 to 3
(inclusive).
A somewhat smaller range of 1.5 to 3 (inclusive) is recommended by
[Tukey77].
Called in a list context, the method returns the reference to the
trimmed list, as well as references to lists of low-end and high-end
values that were filtered:

	($list, $lo_rejects, $hi_rejects) = tukey_iqr ($list);

This method may fail if much of the data set consists of identical
values, such that quartiles 1 and 3, i.e. Q1 and Q3, are
equal.
Here, the IQR will necessarily be zero, so all elements of the
dataset I<not> equal to the median will be identified as outliers
and subsequently removed, regardless of the value of the IQR
coefficient.
The function will L<carp|Carp> in this case, but will nevertheless execute the
calculation.

=cut

#  LOCAL PROGRAMMER NOTES:
#
#  (1) CERTIFICATION: THIS IMPLEMENTATION GETS THE EXACT SAME RESULT AS SEVERAL
#      EXAMPLE PROBLEMS
#
#      * https://www.purplemath.com/modules/boxwhisk3.htm
#
#  (2) WE USE A METHOD CALL TO Statistics::Descriptive TO GET THE QUARTILES
#      AND THE PACKAGE DOES ITS OWN INTERNAL SORTING. WE HAVE TO REPEAT SORTING
#      HERE MANUALLY TO CLEAN THE LIST, WHICH IS NOT VERY EFFICIENT --- CHANGE
#      THIS WHEN THERE IS TIME (10-DEC-2018)

sub tukey_iqr {
   my ($list, $coefficient) = @_;
   $coefficient = 1.5 unless defined $coefficient;

#__ARGUMENT MUST BE REFERENCE TO LIST OF NUMBERS
   my $error = _check_args_ ($list);
   croak $error if defined $error;

#__CHECK COEFFICIENT TOO
   $error = _check_args_ ([$coefficient]);
   croak $error if defined $error;

#__WARN USER IF COEFFICIENT CHOICE IS QUESTIONABLE
   if ($coefficient < 1) {
      carp "coefficient < 1 may result in too many non-outliers being removed";
   } elsif ($coefficient > 3) {
      carp "coefficient > 3 may be too restrictive for removing outliers";
   }

#__SORT LIST (EVEN THOUGH THE QUARTILE SUBROUTINE WILL DO THAT AGAIN: SEE NOTES)
   @{$list} = sort _numerical_ @{$list};

#__GET QUARTILES
   my $stat = Statistics::Descriptive::Full->new();
   $stat->add_data (@{$list});
   my $q1 = $stat->quantile (1);
   my $q3 = $stat->quantile (3);

#__INTER-QUARTILE RANGE SCALED BY THE TUKEY COEFFICIENT
   my $scaled_iqr = $coefficient * ($q3 - $q1);
   if ($scaled_iqr == 0) {
      carp "scaled IQR is zero: too many identical elements, try diff method?";
   }

#__TALLY OUTLIERS ON THE LOW AND HIGH ENDS RELATIVE TO THE TUKEY "HINGES"
   my ($thresh_low, $thresh_high) = ($q1 - $scaled_iqr, $q3 + $scaled_iqr);
   my ($num_low_outs, $num_high_outs) = (0, 0);

#__TALLY THE LOW END
   foreach my $element (@{$list}) {
      if ($element < $thresh_low) {
         $num_low_outs++;
      } else {
         last;
      }
   }

#__TALLY THE HIGH END
   foreach my $element (reverse @{$list}) {
      if ($element > $thresh_high) {
         $num_high_outs++;
      } else {
         last;
      }
   }

#__ELIMINATE OUTLIERS FROM THE DATA ACCORDING TO TALLIES
   my ($low_end_rejects, $high_end_rejects) = ([], []);
   @{$low_end_rejects} = splice (@{$list}, 0, $num_low_outs) if $num_low_outs;
   @{$high_end_rejects} = splice (@{$list}, -$num_high_outs) if $num_high_outs;

#__RETURN CLEANED LIST
   if (wantarray) {
      return ($list, $low_end_rejects, $high_end_rejects);
   } else {
      return $list;
   }
}

# ------------------------------------------------------------------------------

=head2 Median Absolute Deviation (MAD) Method

This non-parametric approach is based on the median of an auxiliary
set of numbers that are themselves the absolute values of the deviation
between the median of the original set and each original datum
[Rousseeuw93].

	$list = median_absolute_dev ($list);

Internally, the method is scaled to be consistent
with what a normal distribution of the data would
imply.
Without getting into technicalities, it means the threshold for deciding
whether a datum is an outlier can be set in a similar context as to how
it is done in terms of standard deviations of the mean for actual normal
data.
The default threshold is 3, but the method takes an optional second
argument to change this threshold:

	$list = median_absolute_dev ($list, 2);

Note that the method will L<carp|Carp> unless the threshold is within
the range of 2 to 3 (inclusive), which is the typical range of usage
[Rousseeuw93].
Called in a list context, the method returns the reference to the
trimmed list, as well as references to lists of low-end and high-end
values that were filtered:

	($list, $lo_rejects, $hi_rejects) = median_absolute_dev ($list);

As with Tukey's method, this procedure may fail
if much of the data set consists of identical
values.
The function will likewise L<carp|Carp> in
this case, but will nevertheless execute the
calculation.

=cut

#  LOCAL PROGRAMMER NOTES:
#
#  (1) CERTIFICATION: SOLVES THE EXAMPLE AT
#      eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers

sub median_absolute_dev {
   my ($list, $threshold) = @_;
   $threshold = 3 unless defined $threshold;

#__ARGUMENT MUST BE REFERENCE TO LIST OF NUMBERS
   my $error = _check_args_ ($list);
   croak $error if defined $error;

#__CHECK THRESHOLD TOO
   $error = _check_args_ ([$threshold]);
   croak $error if defined $error;

#__WARN USER IF THRESHOLD CHOICE IS QUESTIONABLE
   if ($threshold < 2) {
      carp "threshold < 2 may result in too many non-outliers being removed";
   } elsif ($threshold > 3) {
      carp "threshold > 3 may be too restrictive for removing outliers";
   }

#__CONSISTENCY SCALING COEFFICIENT (E.G. SEE [Rousseeuw93])
#  SEE ALSO: https://en.wikipedia.org/wiki/Median_absolute_deviation
   my $coefficient = 1.4826;

#__SORT SAMPLE LIST
   @{$list} = sort _numerical_ @{$list};

#__GET SAMPLE MEDIAN
   my $median_sample = _median_ ($list);

#__GET LIST OF ABSOLUTE DEVIATIONS FROM THE SAMPLE MEDIAN
   my $devs_from_median = [];
   foreach my $val (@{$list}) {
      push @{$devs_from_median}, abs ($val - $median_sample);
   }

#__SORT DEVIATIONS LIST
   @{$devs_from_median} = sort _numerical_ @{$devs_from_median};

#__GET MEDIAN OF DEVIATIONS I.E. THE "MEDIAN ABSOLUTE DEVIATION" (MAD)
   my $mad = _median_ ($devs_from_median);

#__RENDER A CONSISTENT ESTIMATOR BY PROPER SCALAING (AGAIN, SEE [Rousseeuw93])
   $mad *= $coefficient;

#__CARP IF MAD=0 (TOO MANY IDENTICAL VALS) AND SET TO EPSILON TO AVOID DIV-BY-0
   unless ($mad) {
      carp "MAD is zero: too many identical elements, try diff method?";
      $mad = 0.000000000000001;
   }

#__DETERMINE OUTLIERS BY COMPARING ABSOLUTE DEVIATION OF SAMPLE DATA FROM
#  SAMPLE MEDIAN AS SCALED BY MAD VERSUS THE THRESHOLD --- STORE THIS IN AN
#  INDICATOR LIST (0 = NOT OUTLIER, 1 = OUTLIER) 
   my @indicators = ();
   foreach my $val (@{$devs_from_median}) {
      my $mad_scaled_distance = $val / $mad;
      if ($mad_scaled_distance > $threshold) {
         push @indicators, 1;
      } else {
         push @indicators, 0;
      }
   }

#__ELIMINATE OUTLIERS FROM THE DATA ACCORDING TO ACCOMPANYING INDICATOR ARRAY
   my ($num_low_outs, $num_high_outs) = _read_indicators_ (@indicators);
   my ($low_end_rejects, $high_end_rejects) = ([], []);
   @{$low_end_rejects} = splice (@{$list}, 0, $num_low_outs) if $num_low_outs;
   @{$high_end_rejects} = splice (@{$list}, -$num_high_outs) if $num_high_outs;

#__RETURN CLEANED LIST
   if (wantarray) {
      return ($list, $low_end_rejects, $high_end_rejects);
   } else {
      return $list;
   }
}

################################################################################
##                                                                            ##
##             P R I V A T E   R O U T I N E S   A R E   H E R E              ##
##                                                                            ##
################################################################################

#  ---------
#  GOULD'S X
#  ---------
#
#  LOCAL PROGRAMMER NOTES:
#
#  THIS ROUTINE IS PATTERNED DIRECTLY AFTER A PYTHON SCRIPT IN THE WIKIPEDIA
#  PAGE ON PEIRCE'S CRITERION en.wikipedia.org/wiki/Peirce%27s_criterion
#  (IN FACT, WE EVEN DUPLICATE THEIR CODING COMMENTS, WITH SOME AUGMENTATION),
#  WHICH IS ITSELF PATTERNED DIRECTLY AFTER GOULD'S 1855 RE-STATEMENT [Gould55]
#  OF PEIRCE'S ORIGINAL METHOD [Peirce52] TO CALCULATE THE COEFFICIENT OF
#  MAXIMUM ALLOWABLE DEVIATION ("X" AND "R" SEEM TO BOTH BE USED TO REPRESENT
#  THIS VALUE) WHICH IS USED IN THE LARGER OUTLIER IDENTIFICATION PROCESS.

sub _goulds_x_ {
   my ($num_obs, $num_outs) = @_;
   my $x_sqr;

#__THIS IMPLEMENTATION ASSUMES M=1 PROCESS THAT CREATES NOISE THEREFORE OUTLIERS
   my $m = 1;

#__CALCULATE SQUARED THRESHOLD ERROR DEVIATION FOR COHORT SIZE AND # OUTLIERS
   if ($num_obs > 1) {

   #__CALCULATE Q (THE N-TH ROOT OF GOULD'S EQUATION B):
      my $q = (
         $num_outs**($num_outs / $num_obs) *
         ($num_obs - $num_outs)**(($num_obs - $num_outs)/$num_obs)
      ) / $num_obs;

   #__INITIALIZE R VALUES
      my ($r_new, $r_old) = (1, 0);

   #__ITERATE TO CONVERGE ON R
      while (abs($r_new - $r_old) > 0.0000000000000002 * $num_obs) {

      #__CALCULATE LAMBDA (THE 1/(N-n)-TH ROOT OF GOULD'S EQUATION A')
         my $l_div = $r_new**$num_outs;
         $l_div = 0.000001 if $l_div ==0;
         my $lambda = ($q**$num_obs/$l_div)**(1/($num_obs - $num_outs));

      #__CALCULATE X-SQUARED (GOULD'S EQUATION C)
         $x_sqr = 1 + ($num_obs - $m - $num_outs)*(1 - $lambda**2)/$num_outs;

      #__NUMERICAL ARTIFACT: DO NOT ALLOW NEGATIVE X-SQUARED
         if ($x_sqr < 0) {
            $x_sqr = 0;
            $r_old = $r_new;

      #__OTHERIWSE USE X-SQUARED TO UPDATE R (GOULD'S EQUATION D)
         } else {
            $r_old = $r_new;
            my $argument = sqrt ($x_sqr / 2);
            $r_new = exp (($x_sqr - 1)/2) * _erfc_ ($argument);
         }
      }
   } else {
      return 0;
   }

#__RETURN FINAL GOULD X
   return sqrt($x_sqr);
}

#  ---------------------------------
#  ERROR FUNCTION AND ITS COMPLEMENT
#  ---------------------------------
#
#  LOCAL PROGRAMMER NOTES:
#
#  THE METHOD USED HERE IS A CURVE FIT FROM ABRAMOWITZ AND STEGUN 7.1.26
#  (BIBTEX KEY abramowitz:1972::math:text)
#
#  MODELED AFTER PUBLIC DOMAIN PYTHON CODE AT THE LINK BELOW (WITH SOME MODS)
#  https://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/

sub _erf_ {
   my ($x) = @_;

#__FIRST TERM OF TAYLOR SERIES IE 2*x/sqrt(pi) IF ARGUMENT IS VERY SMALL
   return 1.1283791670955126 * $x if abs ($x) < 0.001;

#__CURVE FIT CONSTANTS
   my ($a1, $a2, $a3, $a4, $a5, $p) = qw/0.254829592 -0.284496736 1.421413741
                                         -1.453152027 1.061405429 0.3275911/;
#__NOTE THE SIGN OF THE ARGUMENT
   my $sign = 1;
   $sign = -1 if $x <0;
   $x = abs ($x);

#__ABRAMOWITZ AND STEGUN FIT
   my $t = 1/(1 + $p * $x);
   my $y = 1 - ((((($a5*$t + $a4)*$t) + $a3)*$t + $a2)*$t + $a1)*$t*exp(-$x*$x);
   return $sign * $y;
}

sub _erfc_ {
   my ($x) = @_;
   return 1 - _erf_ ($x);
}

#  -------------
#  SIMPLE MEDIAN
#  -------------
#
#  ASSUMES THE ARGUMENT IS A REFERENCE TO A ** SORTED ** LIST

sub _median_ {
   my ($list) = @_;
   my $length = scalar @{$list};
   my $index = int ($length/2);

#__IF LIST HAS ODD NUMBER OF ELEMENTS THEN THE MIDDLE IS THE MEDIAN
   if ($length % 2) {
       return $list->[$index];

#__ELSE FOR EVEN RETURN THE AVERAGE OF THE 2-EQUALLY-MIDDLE-VALUES
   } else {
       return ($list->[$index - 1] + $list->[$index]) / 2;
   }
}

#  -----------------
#  NUMERICAL SORTING
#  -----------------

sub _numerical_ {$a <=> $b}

#  ----------------------
#  SIMPLE ARITHMETIC MEAN
#  ----------------------

sub _mean_ {
   my ($list) = @_;
   my $mean = 0;
   my $num_vals = scalar @{$list};
   foreach my $number (@{$list}) {
      $mean += $number;
   }
   $mean /= $num_vals;
   return $mean;
}

#  ---------------------------
#  MEAN AND STANDARD DEVIATION
#  ---------------------------

sub _mean_and_standev_ {
   my ($list) = @_;

#__ARGUMENT MUST BE REFERENCES TO LISTS OF NUMBERS
   my $error = _check_args_ ($list);
   croak $error if defined $error;

#__CENTRAL LOCATION IS THE MEAN
   my $x_1_bar = _mean_ ($list);

#__SUM OF SQUARES OF DIFFEFERENCES FROM MEAN
   my $sum_sq_of_diffs_1 = _sum_square_of_diffs_ ($list, $x_1_bar);

#__UNBIASED ESTIMATOR IS SIZE MINUS 1
   my $denominator = scalar @{$list} - 1;

#__SAMPLE VARIANCE
   my $s_sqr;
   if ($denominator) {
      $s_sqr = $sum_sq_of_diffs_1 / $denominator;
   } else {
      $s_sqr = 0;
   }

#__RETURN MEAN AND STANDARD DEVIATION
   return ($x_1_bar, sqrt ($s_sqr));
}

#  -----------------------------
#  SUM OF SQUARES OF DIFFERENCES
#  -----------------------------

sub _sum_square_of_diffs_ {
   my ($list, $central_location) = @_;
   my $sum_square_of_diffs = 0;
   foreach my $number (@{$list}) {
      my $sq_diff = ($number - $central_location)*($number - $central_location);
      $sum_square_of_diffs += $sq_diff;
   }
   return $sum_square_of_diffs;
}

#  ----------------------------------------------
#  CONVERT INDICATOR ARRAY TO TALLIES OF OUTLIERS
#  ----------------------------------------------

sub _read_indicators_ {
   my (@indicators) = @_;
   my ($num_low, $num_high) = (0, 0);

#__HOW MANY OUTLIERS ARE ON THE LOW END
   foreach my $indicator (@indicators) {
      if ($indicator) {
         $num_low++;
      } else {
         last;
      }
   }

#__HOW MANY OUTLIERS ARE ON THE HIGH END
   foreach my $indicator (reverse @indicators) {
      if ($indicator) {
         $num_high++;
      } else {
         last;
      }
   }

#__RETURN TALLIES
   return ($num_low, $num_high);
}

#  ------------------
#  VALIDATE ARGUMENTS
#  ------------------

sub _check_args_ {
   my ($list) = @_;

#__ARGUMENT MUST BE REFERENCE TO LIST OF NUMBERS (NOTE REGEXP::COMMON
#  DOES NOT SEEM TO PROVIDE A GENERAL REGEXP FOR FLOATING POINTS OF VARIOUS
#  FORMS (AS ITS TITLE WOULD SEEM TO INDICATE) --- FOUND A PRETTY GOOD REGEXP
#  AT http://perl.active-venture.com/pod/perlretut-regexp.html WHICH SEEMS TO
#  BE PART OF THE PERL REGULAR EXPRESSIONS TUTORIAL
   return "argument must be list ref" unless ref $list eq "ARRAY";
   return "argument list must have >= 1 element" unless scalar @{$list};
   foreach my $string (@{$list}) {
      return "'$string' not a number" unless
         $string =~ /^[+-]?\ *(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$/;
   }

#__RETURN NO ERROR MESSAGE IF EVERYTHING IS OK
   return undef;
}

################################################################################
##                                                                            ##
##             T R A I L I N G   P O D   D O C U M E N T A T I O N            ##
##                                                                            ##
################################################################################

=head1 REFERENCES

=over

=item *

[Dixon50] WJ Dixon (1950)
I<Analysis of Extreme Values>,
Annals of Mathematical Statistics B<21>(4), 488-506.

=item *

[Gould55] BA Gould (1855)
I<On Peirce's Criterion for the Rejection of Doubtful Observations, With
Tables for Facilitating Its Application>,
Astronomical Journal B<4>(11), 81-87.

=item *

[McGill78] R McGill, JW Tukey, and WA Larsen (1978)
I<Variation of Box Plots>, The American Statistician
B<32>(1), 12-16.

=item *

[Peirce52] B Peirce (1852)
I<Criterion for the Rejection of Doubtful Observations>,
Astronomical Journal B<2>(21), 161-163.

=item *

[Ross03] SM Ross (2003)
I<Peirce's Criterion for the Elimination of Suspect Experimental Data>,
Journal of Engineering Technology B<20>(2), 38-41.

=item *

[Rousseeuw93] PJ Rousseeuw and C Croux (1993)
I<Alternatives to the Median Absolute Deviation>,
Journal of the American Statistical Association B<88>(424), 1273-1283.

=item *

[Sokal95] RR Sokal and FJ Rohlf (1995) I<Biometry>, WH Freeman and Co.

=item *

[Tukey77] JW Tukey (1977) I<Exploratory Data Analysis>, Addison-Wesley.

=back

=head1 AUTHOR

Michael C. Wendl

S<mwendl@wustl.edu>

Copyright (C) 2018 Washington University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

=cut

################################################################################
##                                                                            ##
##                                -  E N D  -                                 ##
##                                                                            ##
################################################################################

1;
