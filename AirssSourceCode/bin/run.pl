#!/usr/bin/env perl 

#==================================================================================#
#                                   run.pl                                         #
#==================================================================================#
#                                                                                  #
# This file is part of the AIRSS structure prediction package.                     #
#                                                                                  #
# AIRSS is free software; you can redistribute it and/or modify it under the terms #
# of the GNU General Public License version 2 as published by the Free Software    #
# Foundation.                                                                      #
#                                                                                  #
# This program is distributed in the hope that it will be useful, but WITHOUT ANY  #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  #
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.        #           
#                                                                                  #
# You should have received a copy of the GNU General Public License along with this#
# program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street,#                   
# Fifth Floor, Boston, MA  02110-1301, USA.                                        #
#                                                                                  #
#----------------------------------------------------------------------------------#
# This script relaxes a bunch of structure files using CASTEP                      #
#----------------------------------------------------------------------------------#
# Written by Chris Pickard, Copyright (c) 2005-2018                                #
#----------------------------------------------------------------------------------#
#                                                                                  #
#==================================================================================#

use strict;
use Getopt::Long;
use File::Copy;
use List::Util qw ( shuffle );
use Time::HiRes qw ( sleep );

sub usage {
  printf STDERR "Usage: run.pl [-mpinp] [-conventional] [-keep] [-check] Run many Castep computations\n";
  printf STDERR "       -mpinp          Number of cores per mpi Castep (0)\n";
  printf STDERR "       -conventional   Use conventional cell - cif only (0)\n";
  printf STDERR "       -keep           Keep all output files (0)\n";
  printf STDERR "       -check          Check crystal structure files (0)\n";
  exit();
}

my ($opt_mpinp,$opt_conventional,$opt_keep,$opt_check,$opt_help) = (0,0,0,0,0);

&GetOptions("mpinp=n"      => \$opt_mpinp,
            "conventional" => \$opt_conventional,
	    "keep"         => \$opt_keep,
	    "check"        => \$opt_check,
            "h|help"       => \$opt_help)|| &usage();

if ($opt_help) {
  &usage();
}
;

# Loop over the structure files, use SHLX if available, otherwise CIF

my @files = <*-*.{res,cif}>;

foreach (@files) {
  
  my $file = $_;
  my @tmp = split('\.',$_);
  my $seed=$tmp[0];
  my $type=$tmp[1];
  @tmp = split('-',$seed);
  my $root=$tmp[0];

  my $duplicate=0;
  
  if (-e "jobs.txt") {
    open(JOBFILE, "+<jobs.txt") || die;
    while (<JOBFILE>) {
      if ($_ eq $seed."\n") {
  	$duplicate=1;
      }
    }
    if ($duplicate == 0) {
      print JOBFILE $seed."\n";
    }
    close JOBFILE || die;
  } else {
    open(JOBFILE, ">jobs.txt") || die;
    print JOBFILE $seed."\n";	
    close JOBFILE || die;
  }  
  
  if ($duplicate == 0) {

    if ($opt_check) {

      if ($type eq "cif") {
	chomp(my $status = `cif2cell $file -p castep -o /dev/null 2>&1 | grep Error | awk '{print \$1}'`);
	if (index($status,"Error") > -1) {
	  system("mkdir -p bad_cif ; mv $file bad_cif");
	} 
      }

    } else {
      
      if ($type eq "res") {
	system("(cabal res cell < $file ; sed -e '/^%BLOCK [Ll][Aa][Tt]*/, /^%ENDBLOCK [Ll][Aa][Tt]*/d' $root.cell | sed -e '/^%BLOCK [Pp][Oo][Ss]*/, /^%ENDBLOCK [Pp][Oo][Ss]*/d') > $seed.cell");
      }
     
      my $reduce='';
      if ($opt_conventional == 1) {
	$reduce = '--no-reduce';
      }

      if ($type eq "cif") {
	system("cif2cell $file -p castep -o $seed.tmp.cell $reduce &> /dev/null ; (cat $seed.tmp.cell | sed -n -e '/^%BLOCK [Ll][Aa][Tt]*/, /^%ENDBLOCK [Pp][Oo][Ss]*/p' | grep -v angstrom | sed 's/D /H /g'; sed -e '/^%BLOCK [Ll][Aa][Tt]*/, /^%ENDBLOCK [Ll][Aa][Tt]*/d' $root.cell | sed -e '/^%BLOCK [Pp][Oo][Ss]*/, /^%ENDBLOCK [Pp][Oo][Ss]*/d') > $seed.cell ; rm $seed.tmp.cell");
      }

      my $executable= 'castep'; if ($opt_mpinp > 0 ) {$executable= '"mpirun -np '.$opt_mpinp.' castep"'}
            
      #system("ln -sf $root.param $seed.param ; eval $executable $seed");
      system("cp $root.param $seed.param");

      # Ensure that the spin in the cell and param are consistent

      my $spintot=`grep SPIN= $seed.cell| awk 'BEGIN {FS="SPIN="};{ sum += \$2 } END {printf "%10.3f",sum}'`;

      open  PARAMFILE, "$seed.param" or die $!;
      my @paramdata = <PARAMFILE>;
      close PARAMFILE;
      
      open  PARAMFILE, ">$seed.param" or die $!;
      foreach (@paramdata) {
	my @vec = split(' ',$_);
	if ( lc $vec[0] ne "spin" ) {
	  print PARAMFILE $_;
	}
	;
      }
      print PARAMFILE "spin : ".$spintot;
      close PARAMFILE;
      
      if ($opt_mpinp > 0 ) {
	system("eval $executable $seed");
      } else {
	system("$executable $seed");
      }
      
      if (-e $seed."-out.cell") {
	
	if (-e $seed.".dome_bin") {
	  system("(echo 'task : dos';echo 'broadening : linear'; echo 'compute_band_gap : true')>$seed.odi;optados $seed;rm $seed.linear.dat;sed 's/legend 0.85, 0.8/legend 0.2, 0.8/g' $seed.linear.agr | sed 's/Electronic Density of States/$seed/g' > $seed.dos.agr ; (echo '\@target G0.s1'; echo '\@type xy' ; echo 'EfD 0'; echo 'EfD 5' ; echo '&' ;echo ' @    xaxis ticklabel font 4';echo ' @    yaxis ticklabel font 4') >> $seed.dos.agr");
	}
	
	# Construct the res file

	system("castep2res $seed > $seed.res ; rm $seed-out.cell");
	
	# Add data to DOS plot
	
	if (-e $seed.".dos.agr") {
	  my $efd='';
	  if (-e $seed.".odo") {
	    chomp($efd = `grep EfD $seed.odo | tail -1 | awk '{print \$7}'`);
	  }
	  system("mv $seed.dos.agr $seed.dos.tmp.agr ; sed 's/EfD/$efd/g' $seed.dos.tmp.agr | sed 's/Generated by OptaDOS//g' > $seed.dos.agr ; rm $seed.dos.tmp.agr ")
	}
	  
      } else {
	system("mkdir -p bad_castep ; mv $seed.* bad_castep")
      }
	
      if (! $opt_keep) {
	my @delete_files = <$seed.*>;
	foreach (@delete_files) {
	  my @tmp = split('\.',$_); my $type=$tmp[1];
	  if (($type ne "res")&&($type ne "cif")&&($type ne "magres")&&($type ne "castep")&&($type ne "odo")&&($type ne "dos")) {
	    unlink($_);
	  }	
	}
      }

      #system("mkdir -p good_castep ; mv $seed.* good_castep")
	
    }
  }
    
}
  
  
