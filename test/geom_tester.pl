#!/usr/bin/perl

# geom_tester.pl
# Some driver code to test the Geometric.pm package

use strict;
use warnings;
use Geometric;

my @ps = qw(0.2 0.5 0.8);
my @ks = qw(0 1 2 3 4 5 6 7 8 9 10);

foreach my $p ( @ps )	{
	print "p = $p\n";
	my $Geom = Geometric->new( "p" => $p );
	foreach my $k ( @ks )	{
		my $prob = $Geom->probability($k);
		print "k = $k    --- GeomPr(X = k) = $prob\n";
	}
	print "\n\n";
}

# Estimate K from a GeomProb and a given p
print "\n\nEstimate K --> GeomPr = 0.16, p=0.8\n";		# K should equal 1
my $Geom = Geometric->new( "p" => 0.8 );
my $estk = $Geom->estimate_k(0.16);
print " === K = $estk === \n\n";

# Maximum Likelihood Estimation of p from a list of k-sampled trials
print "\n\nMLE estimation of p from a list of samples, which are the number of k failures before the first success\n";	
my @samples = qw(1 1 1 1 1 2);
my $p = Geometric->geom_mle( \@samples );
print " === P =~ $p === \n\n";

# Draw random samples from a Geometric distribution
$Geom->{_p} = 0.8;
my $random_sample = $Geom->rgeom(75);
print "\n\nRandomyl sample a Geometric distribution with parameter p=0.8\n";
print join(" ", @{$random_sample}), "\n\n";

exit;