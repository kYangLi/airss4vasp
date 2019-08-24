#!/usr/bin/env perl

use Getopt::Long;
use strict;

sub usage {
   printf STDERR "Usage: tidy.pl Tidy up the result of a AIRSS run\n";
   exit();
}

my ($opt_help) = (0);

&GetOptions("h|help" => \$opt_help)|| &usage();

if ($opt_help) {&usage()};

print "Files will be removed - <ENTER> to continue";
my $ctemp = <STDIN> ;
if (-e "jobs.txt") {print "jobs.txt detected - no files deleted\n";exit()}

my @files = <*-*.cell>;

foreach (@files) {
    my @tmp = split('\.',$_);
    my $seed=$tmp[0];
    
    if (-e "$seed.res") {
	open RESFILE, "$seed.res";
	while (<RESFILE>) {
	    my @tmp = split(' ',$_);
	    if ($tmp[1] eq "Converted") {unlink glob("$seed.*")}
	}
    } else {
	unlink glob("$seed.*");
    }

}
unlink glob("*.conv");
system "rm -fr trash symmetry* ToMatServer.txt stdout.txt StartUp* MatStudioServer.log findsym.log *.temp gulptmp_* *.minsep";
