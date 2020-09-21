#!/usr/bin/perl
#
#  SINGLE-USE SCRIPT TO GET BULK EXPRESSION DATA OVER GTEX SAMPLES
#  FOR JUST THE GENES IN TABLE 7A OF THE MM PAPER
 
#use Getopt::Long;

# save arguments following -h or --host in the scalar $host
# # the '=s' means that an argument follows the option
# # they can follow by a space or '=' ( --host=127.0.0.1 )
#GetOptions( 'infile=i' => \my $infile
#            ,'gtexdir=g' => \my $gtexdir  # same for --user or -u
#                               );
#
############
#  SET UP  #
############

#__REGULAR STUFF
   use strict;
   use warnings;

 
# save arguments following -h or --host in the scalar $host
# # # the '=s' means that an argument follows the option
# # # they can follow by a space or '=' ( --host=127.0.0.1 )
#   use Getopt::Long;
#GetOptions('infile=i' => \my $infile,
#           'gtexdir=g' => \my $gtexdir);
my ($infile) = @ARGV;

if (not defined $infile) {
  die "Need path of input gene list\n";
}

###################
#  CONFIGURATION  #
###################
#__GENE LIST Folder
#   my $marker_dir = "/diskmnt/Datasets/mmy_scratch/lyao/MMY/Analysis/cell_surface_markers/Scripts/V8/automate_test/";
#
#__FILE FROM WHICH WE GET GENE NAMES
   #my $file_genes = "Gene_list_of_DE_genes_filtered_surface_3DB.txt";
   #my $file_genes = join (
   #   '/', $indir, "Gene_list_of_DE_genes_filtered_surface_3DB.txt"
   #);
#__GTEX HOME
   #my $gtex_dir = "/diskmnt/Datasets/GTEX_tpm";

#__GTEX DATA FILE
   my $gtex_data_file = join (
      '/', "DOWNLOAD_DB", "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"
   );

#__GTEX ANNOTATION FILE
   my $gtex_annot_file = join (
      '/', "DOWNLOAD_DB", "GTEx_v7_Annotations_SampleAttributesDS.txt"
   );

####################
#  PRE-PROCESSING  #
####################

#__GET THE GENE LIST
   my $genes = &genes ($infile);
   my $num_genes = scalar keys %{$genes};

###############
#  MAIN CODE  #
###############
 
#__GET RUIYANG GENES
   my $count = 0;
   open (F, $gtex_data_file) || die "cant open file '$gtex_data_file'";
   while (<F>) {
      print if /^Name/;
      next unless /^ENSG/; # THESE ARE DATA LINES IN THE ANALYSIS FILE
      my @fields = split /\t/;
      my $gene = $fields[1];
      if (defined $genes->{$gene}) {
         $count++;
         $genes->{$gene} = 0;
         print;
      }
   }
   print "# CAPTURED $count OF $num_genes DE GENES\n";
   close (F);
   foreach my $gene (keys %{$genes}) {
      print "# not captured: $gene\n" if $genes->{$gene};
   }

#####################
#  POST-PROCESSING  #
#####################

################################################################################
#                                                                              #
#                                  SUBROUTINES                                 #
#                                                                              #
################################################################################
 
#__READ GENE LIST
   sub genes {
      my ($file) = @_;
#my @local_genes = qw/
#PMEL CD63 S100A13 GAPDH GPNMB GSTP1 SFRP1 FXYD3 ATP1A1 NDRG1 BSG HSPB1 BAALC PRAME VIM RAB13 MARCKSL1 CAV1 PLP1 BACE2 CCT3 APP EDNRB RAC1 SYNGR1 CTNNB1 FXYD1 HSP90AB1 PMP22 ATP5B YWHAE BRI3 PERP SCCPDH ERBB3 C1QBP GPM6B BAMBI FLOT1 CADM1 RHOBTB3 TNFRSF12A AP2S1 ALDH1A3 GPR56 PLEKHA4 GNG5 GNG11 LAMTOR2 CD151 STOML2 PHB HSPD1 PRKAR1A PRNP ETV5 VDAC1 NRP2 GPI GAS7 EMP2 PHLDA3 HDLBP ADI1 RDX RAB7A NEDD4L GNB2 NDUFC2 PDIA6 TRIP6 SDHB SDC3 ADRM1 SOD1 PON2 HAX1 LY6K PPFIBP1 HM13 CNP SEPT2 NUDT1 COMT PTPRZ1 WLS ATP6V1D SCARB1 IFNGR2 NELFE METAP2 ATP1B3 PDAP1 CAPN2 COA6 SERBP1 APH1A ADSS THEM4 ANXA4 ATRAID CCDC124 MLEC PARVB SRI MESDC2 LTBR RAB5C DCXR WASL STXBP6 EPN1 STX3 PPP2R1A ABHD12 VAPA TMEM179B CTNND1 TNFRSF1A AAMP /;

   #__READ GENES FROM FILE
      my $genes = {};
      open (F, $file) || die "cant open file $file";
      while (<F>) {
         chomp;
         $genes->{$_} = 1;         
         print "GENE IS -->$_<--\n";
      }
      close (F);
#     while (<DATA>) {}
#     foreach my $gene (@local_genes) {
#        chomp;
#        $genes->{$gene} = 1;
#     }
      return $genes;
   }

################################################################################
#                                                                              #
#                                     DATA                                     #
#                                                                              #
################################################################################
#
#  THIS IS THE GENE LIST FOR COMPARISON FROM FIGURE 7A OF THE MM PAPER

