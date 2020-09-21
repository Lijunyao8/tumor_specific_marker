#!/usr/bin/perl
#
#  THIS SCRIPT IS FOR FINDING TISSUE-SPECIFIC EXPRESSION OF GENES FROM A
#  GENE LIST FROM LIJUN IN TERMS OF RNA EXPRESSION DATA FROM THE HUMAN PROTEIN
#  ATLAS
#
#  ######################
#  #  PROGRAMMER NOTES  #
#  ######################
#
#  -------
#  GENERAL
#  -------
#
#  -----------------------
#  FILES AND THEIR FORMATS
#  -----------------------
#
############
#  SET UP  #
############

#__STANDARD PRAGMAS
   use strict;
   use warnings;

my ($infile, $hpafile) = @ARGV;

if (not defined $infile) {
  die "Need path of input gene list\n";
}

if (not defined $hpafile) {
  die "Need path of expression matrix downloaded from HPA\n";
}
#__STANDARD PERL/C-PAN LIBRARIES
   use Statistics::Distributions;
   use Math::SigFigs;

#__SPECIAL LIBRARIES
   use lib "/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/GTEX";
   use Statistics::Outliers;
   use HPA_Read_Tools;

###################
#  CONFIGURATION  #
###################

#__GENE LIST DATA FILE
#  my $infile = "DE_genes_filtered_surface_3DB_GTEX_HPA_updated.txt";
#   my $infile = "DE_genes_filtered_surface_3DB_GTEX_HPA_IPA_review_p0.05.txt";

#__HUMAN PROTEIN ATLAS DATA FILE
#   my $hpafile = "rna_tissue_consensus.tsv";

#__SET GLOBAL FORMATTING VARIABLES FOR OUTPUT
   $^ = "OUTPUT_TOP";
   $~ = "OUTPUT";
   $- = 0; # closes out any previous report so current "TOP" will be used

####################
#  PRE-PROCESSING  #
####################

#__GET GENE LIST
   my $gene_list = HPA_Read_Tools::read_gene_list ($infile);

#__GET EXPRESSION DATA FOR THESE GENES FROM HUMAN PROTEIN ATLAS DATA
   my ($data, $count) = HPA_Read_Tools::read_hpa_data ($hpafile, $gene_list);

#__DIAGNOSTICS: GENES HAVING NO HUMAN PROTEIN ATLAS DATA
#  foreach my $gene (sort keys %{$gene_list}) {
#     print "$gene    $gene_list->{$gene}\n" unless $gene_list->{$gene};
#     print "$gene, " unless $gene_list->{$gene};
#  }
#  print "\n";

###############
#  MAIN CODE  #
###############
#
#  DATA STRUCTURE KEYS ON EXPRESSION VAL WHICH MAY NOT BE UNIQUE WITHIN A GENE
#  push @{$data->{$gene}->{$expression_val}}, $tissue;
#  $data = {
#     CCND1 = {
#        9.7  => [stomach, testis],
#        6.0  => [tonsil],
#        28.0 => [amygdala],
#         :   :     :
#     },
#     :
#  };

#__DISCERN GOOD TEST CANDIDATES USING T-TEST ANALYSIS
   my $highs = 0;
   my ($num_tests, $pval_results) = (0, {});
   foreach my $gene (sort keys %{$data}) {

   #__FORMULATE LIST OF EXPRESSION VALUES FOR THIS GENE: ACCOUNTS FOR DUPLICATES
      my ($values, $tissues) = ([], []);
      foreach my $expr_val (reverse sort _numerical_ keys %{$data->{$gene}}) {
         foreach my $tissue (@{$data->{$gene}->{$expr_val}}) {
            push @{$values}, $expr_val;
            push @{$tissues}, $tissue;
         }
      }

   #__TEST EVERY ELEMENT IN THE LIST OF EXPRESSIONS FOR THIS GENE
   #
   #  Here we use something unusual, which is to march through the list of
   #  numerical values and, for each, pop out that value from the list using
   #  'splice', test it against the remainders in the list, then pop it back in
   #  at the same original position. This obviates the need to make additional
   #  data structures.
####  print "there are $num values\n";

# @{$values} = qw/clam1 clam2 clam3 clam4/;
# @{$values} = qw/3.02 4.02 3.88 3.34 3.87 3.18/;
# @{$tissues} = qw/eye hand toe nose face boner/;

      for (my $i = 0; $i <= $#{$values}; $i++) {

      #__POP-OUT THE NEXT ELEMENT FOR TESTING
         my $element = splice (@{$values}, $i, 1);
         my $tissue = splice (@{$tissues}, $i, 1);

      #__DESCRIPTIVE STATISTICS OF THE BASELINE (REMAINING) LIST
         my ($mean, $variance) = _mean_and_variance_ ($values);
         my $num = scalar @{$values};

      #__T-STATISTIC FOR 1-SAMPLE TEST (SEE E.G. sokal:1981::math:stats PP 230)
         my $t;
#        if ($element >= $mean) {
            $t = $element - $mean;
            $t /= sqrt ($variance * ($num + 1) / $num);
#        } else {
      ######print "    $gene   $tissue   TEST ($avg_test) < BACKGROUND " .
      ######      "($avg_bgrnd) : NO TEST PERFORMED\n";
#           next;
#        }

      #__DEGREES OF FREEDOM FOR 1-SAMPLE TEST (SEE sokal:1981::math:stats AGAIN)
         my $dof = $num - 1;

      #__HYPOTHESIS TEST
         my $p_value = Statistics::Distributions::tprob ($dof, $t);

      #__SAVE TEST RESULTS: SAVE FOR EACH TISSUE WITH THIS EXPRESSION VALUE
         push @{$pval_results->{$p_value}}, [$gene, $tissue, $element];
         $num_tests++;

# YOUAREHERE
#  print "$i  $tissue\n";
#  print "   element-->$element<----  LIST IS ", join (',', @{$values}), "\n";
#  print "   mean = $mean    variance = $variance\n";
#  print "   T-statistic = $t    DOF = $dof    P = $p_value\n";

      #__POP-IN ELEMENT AND TISSUE TO THEIR ORIGINAL POSITIONS
         splice (@{$values}, $i, 0, $element);
         splice (@{$tissues}, $i, 0, $tissue);
      }
# die "done";
   }

#__MULTIPLE TEST CORRECTION
   my ($rank, $fdr_previous, $min_non_zero, $fdr_non_zero) = (1, 0, 99999, 0);
   foreach my $p_value (sort _numerical_ keys %{$pval_results}) {
      foreach my $list_reference (@{$pval_results->{$p_value}}) {
         my ($gene, $tissue, $element) = @{$list_reference};

      #__ACTUAL FDR CALCULATION
         my $fdr = $p_value * $num_tests / $rank;

      #__MODIFIER 1: FDR DOES NOT EXCEED UNITY
         $fdr = 1 if $fdr > 1;

      #__MODIFIER 2: FDR IS MONOTONIC
         $fdr = $fdr_previous if $fdr < $fdr_previous;
         $fdr_previous = $fdr;

      #__CHOP IT AT THE CUT-OFF
         if ($fdr <= 0.05) {
            $highs++;
         } else {
            last;
         }

      #__RIGMAROLE ROUNDING FOR 2 DIGITS BEHIND THE DECIMAL IN SCI-NOTATION
         my $fdr_round = FormatSigFigs ($fdr, 5);
         my $fdr_sci_notation = sprintf ("%.4e", $fdr_round);

      #__RIGMAROLE ROUNDING FOR 2 DIGITS BEHIND THE DECIMAL IN SCI-NOTATION
         my $p_value_round = FormatSigFigs ($p_value, 5);
         my $p_value_sci_notation = sprintf ("%.4e", $p_value_round);

#__FORMATTING
format OUTPUT_TOP =
  ----------------------------------------------------------------------------
  GENE           TISSUE               P-VALUE        FDR            EXPRESSION
  ----------------------------------------------------------------------------
.
format OUTPUT =
  @<<<<<<<<<<<   @<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<  @<<<<<<<<<<<<  @<<<<<<
  $gene, $tissue, $p_value_round, $fdr_round, $element
.
         write;

      #__OUTPUT AND INCREMENT
      #  print "RANK $rank:  $list_reference->[0]   $list_reference->[1]    " .
      #        "PVAL = $p_value   FDR = $fdr_sci_notation\n";
         $rank++;
      }
   }


#####################
#  POST-PROCESSING  #
#####################

#__ALSO OUTPUT SUMMARY STATISTICS
   print "total number of expression values examined: $count\n";
   print "number deemed to be high outliers: $highs\n";

################################################################################
#                                                                              #
#                                  SUBROUTINES                                 #
#                                                                              #
################################################################################

#__STANDARD PERL NUMERICAL SORT
   sub _numerical_ {$a <=> $b}

#############################
#  SUROUTINES:  STATISTICS  #
#############################

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
#  MEAN AND VARIANCE
#  ---------------------------

sub _mean_and_variance_ {
   my ($list) = @_;

#__CENTRAL LOCATION IS THE MEAN
   my $x_1_bar = _mean_ ($list);

#__SUM OF SQUARES OF DIFFEFERENCES FROM MEAN
   my $sum_sq_of_diffs_1 = _sum_square_of_diffs_ ($list, $x_1_bar);

#__UNBIASED ESTIMATOR IS SIZE MINUS 1
   my $denominator = scalar @{$list} - 1;

#__SAMPLE VARIANCE
   my $variance;
   if ($denominator) {
      $variance = $sum_sq_of_diffs_1 / $denominator;
   } else {
      $variance = 0;
   }

#__RETURN MEAN AND SAMPLE VARIANCE
   return ($x_1_bar, $variance);
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

################################################################################
#                                                                              #
#                                     DATA                                     #
#                                                                              #
################################################################################
