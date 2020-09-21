#!/usr/bin/perl
#
#  SINGLE USE SCRIPT TO GET META-DATA FOR MM PAPER GENE EXPRESSION
#  FROM GTEX: SPECIFICALLY LIST OF TISSUE-TYPES OF THE SAMPLES
 
############
#  SET UP  #
############

#__REGULAR STUFF
   use strict;
   use warnings;

###################
#  CONFIGURATION  #
###################

#__GTEX HOME
   my $gtex_dir = "./DOWNLOAD_DB";

#__GTEX DATA FILE
   my $gtex_data_file = join (
      '/', $gtex_dir, "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"
   );

#__GTEX ANNOTATION FILE
   my $gtex_annot_file = join (
      '/', $gtex_dir, "GTEx_v7_Annotations_SampleAttributesDS.txt"
   );

####################
#  PRE-PROCESSING  #
####################

###############
#  MAIN CODE  #
###############

#__GET META-DATA
   my ($samples, $types) = ({}, {});
   open (F, $gtex_annot_file) || die "cant open file '$gtex_annot_file'";
   while (<F>) {
      chomp;
      next unless /^GTEX/;
      my @fields = split /\t/;
#foreach my $c (@fields) {
#  print "-->$c<--\n";
#}
      my ($sample, $tissue_type) = ($fields[0], $fields[6]);
      $samples->{$sample} = $tissue_type;
      $types->{$tissue_type} = 1;
   }
   close (F);

#__OUTPUT
   print "#  TISSUE TYPES IN THESE DATA\n";
   print "#  ================================\n";
   foreach my $type (sort keys %{$types}) {
      print "#  $type\n";
   }
   print "#  ================================\n";
   print "#  col 1: GTEx sample name\n";
   print "#  col 2: tissue type\n";
   print "#  ================================\n";
   foreach my $sample (sort keys %{$samples}) {
      print "$sample\t$samples->{$sample}\n";
   }

#####################
#  POST-PROCESSING  #
#####################

################################################################################
#                                                                              #
#                                  SUBROUTINES                                 #
#                                                                              #
################################################################################
 
################################################################################
#                                                                              #
#                                     DATA                                     #
#                                                                              #
################################################################################
