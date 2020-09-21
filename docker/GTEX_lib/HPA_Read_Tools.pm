package HPA_Read_Tools;

################################################################################
##                                                                            ##
##         I N T R O D U C T O R Y   P O D   D O C U M E N T A T I O N        ##
##                                                                            ##
################################################################################

=head1 NAME

HPA_Read_Tools - File-reading tools for Human Protein Atlas flat files

=head1 SYNOPSIS

	use HPA_Read_Tools;

=cut

#__STANDARD PERL PACKAGES
   use strict;

#__SET UP EXPORTING
   use vars qw( @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION );
   require Exporter;
   @ISA = qw(Exporter);
   @EXPORT = qw();
   @EXPORT_OK = qw/read_gene_list read_hpa_data/;
   %EXPORT_TAGS = (all => [qw/read_gene_list read_hpa_data/]);
   my $pkg = 'HPA_Read_Tools';

#__VARIABLES
   our ($VERSION, $errtxt);
   $VERSION = '0.1';

################################################################################
##                                                                            ##
##                        P U B L I C   R O U T I N E S                       ##
##                                                                            ##
################################################################################

#  --------------------------------
#  READ GENE LIST
#  --------------------------------

sub read_gene_list {
   my ($file) = @_;
   my $gene_list = {};

#__READ FILE CONTAINING GENE LIST
   open (F, $file) || die "cant open file '$file'";
   while (<F>) {
      chomp;
      next if /^#/;

   #__PARSE LINE AND GET GENE
      my @fields = split /\t/;
      my $gene = shift @fields;

   #__SAVE GENE: HASH_VAL=0 (DEFAULT) CHANGES TO 1 ONCE DATA ARE READ
      $gene_list->{$gene} = 0;
#     push @{$gene_list}, $gene;
   }
   close (F);
   return $gene_list;
}

#  --------------------------------
#  READ HUMAN PROTEIN ATLAS DATA
#  --------------------------------
#
#  THIS DATA STRUCTURE KEYS ON THE EXPRESSION VALUE WHICH MAY NOT BE
#  UNIQUE WITHIN A GENE
#  $data = {
#     CCND1 = {
#        9.7  => [stomach, testis],
#        6.0  => [tonsil],
#        28.0 => [amygdala],
#         :   :     :
#     },
#     :
#  };

sub read_hpa_data {
   my ($file, $gene_list) = @_;
   my ($data, $count) = ({}, 0);

#__READ HUMAN PROTEIN ATLAS FILE
   open (F, $file) || die "cant open file '$file'";
   while (<F>) {

   #__READ A DATA LINE
      chomp;
      next if /^#/;
      my ($code, $gene, $tissue, $expression_val) = split /\t/;

   #__SAVE TISSUE AND EXPRESSION VALUE IF THIS GENE IS ON THE GENE LIST
      if (defined $gene_list->{$gene}) {

      #__SAVE THE DATA
######   $data->{$gene}->{$tissue} = $expression_val;
         push @{$data->{$gene}->{$expression_val}}, $tissue;
         $count++;

      #__FLAG THIS GENE AS HAVING HUMAN PROTEIN ATLAS DATA
         $gene_list->{$gene} = 1;
      }
   }
   close (F);
   return ($data, $count);
}

################################################################################
##                                                                            ##
##             P R I V A T E   R O U T I N E S   A R E   H E R E              ##
##                                                                            ##
################################################################################

################################################################################
##                                                                            ##
##             T R A I L I N G   P O D   D O C U M E N T A T I O N            ##
##                                                                            ##
################################################################################

################################################################################
##                                                                            ##
##                                -  E N D  -                                 ##
##                                                                            ##
################################################################################

1;
