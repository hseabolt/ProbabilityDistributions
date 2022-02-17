#!/usr/bin/perl

# package Binomial.pm
# An HS custom Perl package for Binomial distributions, properties, and sampling methods.
# There is no object here, just a methods class.

# Author: MH Seabolt
# Last Updated: 5-26-2020	

package Binomial; 

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(max min sum lg log10 pow round);			 #Import from other packages

use strict;
use warnings;
use List::Util qw( max min sum first shuffle uniq );
use Scalar::Util;
use Storable qw(dclone);
use Data::Dumper;
use Carp;
our $AUTOLOAD;

# @INC libraries for my PC and Biolinux or HS custom classes 
use lib '/media/hunter/Data/scripts';
use lib 'D:\scripts';
use lib '/scicomp/home/ngr8/Biolinux/scripts';


####################################################################################################################
# Class data and methods 
# Attributes 
{
	my %_attribute_properties = (
		_n => '',						# The number of independent Bernoulli trials WITH REPLACEMENT (without replacement is the hypergeometric distribution!)
		_p => '',  						# The probability of a success, given as a positive real number between 0 and 1 inclusive.
		_q => '',						# The probability of a failure, also given as a positive real number between 0 and 1 inclusive.
	);
	
	# Global variable counter
	my $_count = 0;
	
	# Return a list of all attributes
	sub _all_attributes	{
		keys %_attribute_properties;
	}
	
	# Return the default value for a given attribute
    	sub _attribute_default 	{
     	my( $self, $attribute ) = @_;
        	$_attribute_properties{$attribute};
    	}
    
	# Manage the count of existing objects
	sub get_count	{
		$_count;
	}
	sub _incr_count	{
		++$_count;
	}
	sub _decr_count	{
		--$_count;
	}	
}

############################################################
#                       CONSTRUCTORS                       #
############################################################

# The contructor method
# Construct a new graph (my $node = Markov->new() );
# Returns a scalar reference to a
sub new				{
	my ( $class, %arg ) = @_;
	
	# Create the new object
	my $self = bless {}, $class;

	foreach my $attribute ( $self->_all_attributes() ) {
		# E.g. attribute = "_name",  argument = "name"
		my ($argument) = ( $attribute =~ /^_(.*)/ );
		# If explicitly given
		if (exists $arg{$argument}) 	{
			$self->{$attribute} = $arg{$argument};
		}
		else	{
			$self->{$attribute} = $self->_attribute_default($attribute);
		}
   	}
	
	# Construction code specific to this class
	######################################################################################
 
	# Check that we have either P or Q
	if ( not $self->{_p} && not $self->{_q} )	{
		return;
	}
	elsif ( $self->{_p} && not $self->{_q} )	{
		# P must be a positive real number ge 0 and le 1.
		return if ( $self->{_p} < 0 && $self->{_p} > 1.0 ); 
		$self->{_q} = 1 - $self->{_p};
	}
	elsif ( $self->{_q} && not $self->{_p} )	{
		# Q must be a positive real number ge 0 and le 1.
		return if ( $self->{_q} < 0 && $self->{_q} > 1.0 ); 
		$self->{_p} = 1 - $self->{_q};
	}
	else	{
		# Here we have both p and q.  Check that they meet the required criteria and that they sum to 1.0
		return if ( $self->{_p} < 0 && $self->{_p} > 1.0 ); 
		return if ( $self->{_q} < 0 && $self->{_q} > 1.0 ); 
		return if ( $self->{_p} + $self->{_q} != 1 );
	}

	# N must be a positive integer ge 0
	return if ( $self->{_n} < 0 || $self->{_n} != int $self->{_n} );
	
	######################################################################################
	
    $class->_incr_count();
	return $self;
}

# The clone method
# All attributes are copied from the calling object, unless specifically overriden
# Called from an existing object ( Syntax: $cloned_obj = $obj->clone(); )
sub clone	{
	my ( $caller, %arg ) = @_;
	# Extract the class name from the calling object
	my $class = ref($caller);
		
	# Create a new object
	my $self = bless {}, $class;
		
	foreach my $attribute ( $self->_all_attributes() )	{
		my ($argument) = ( $attribute =~ /^_(.*)/ );
			
		# If explicitly given
		if ( exists $arg{$argument} )	{
			$self->{$attribute} = $arg{$argument};
		}
			
		# Otherwise, copy attribute of new object from the calling object
		else	{
			$self->{$attribute} = $caller->{$attribute};
		}
	}
	$self->_incr_count();
	return $self;
}

# When an object is no longer being used, garbage collect it and adjust count of existing objects
sub DESTROY	{
	my ( $self ) = @_;
	$self->_decr_count();
}

# Get and set the value of p and q. Note that if you use the set_p() function, it will recalculate 1 as 1 - p.
sub get_p	{	return $_[0]->{_p};	}
sub p 		{	return $_[0]->{_p};	}
sub set_p	{	
	return if ( $_[1] < 0 || $_[1] > 1.0 );
	$_[0]->{_p} = $_[1];
	$_[0]->{_q} = 1 - $_[1];
}

# Get the value of q, given p.  Note that if you use the set_q() function, it will recalculate p as 1 - q.
sub get_q	{	1 - $_[0]->{_p};	}
sub q 		{	1 - $_[0]->{_p};	}
sub set_q	{	
	return if ( $_[1] < 0 || $_[1] > 1.0 );
	$_[0]->{_q} = $_[1];
	$_[0]->{p} = 1 - $_[1];
}

################################################################################## 

############################################################
#         Properties of the Binomial Distribution      	   #
############################################################

# Some readability has been sacrificed here in favor of computational speed, but in general the equations are given and the order of args is the same that they appear in the given equation.

# Expected value (the expected mean of the distribution, given our parameters)
# I've included a few aliases here, but be very careful that we don't accidentally run into name collisions.
# Expected value (mean) == Np
sub expected	{		$_[0]->{_n} * $_[0]->{_p};		}		

# Variance == Npq
sub variance	{	$_[0]->{_n} * $_[0]->{_p} * $_[0]->{_q};		}

# Median == floor(Np) --> we are just using the int function here to get the floor
sub median		{		int($_[0]->{_n} * $_[0]->{_p});		} 

# Mode == floor((N+1)p)
sub mode 		{  		int($_[0]->{_n} + 1) * $_[0]->{_p};		} 

# Skewness == (q-p) / sqrt(Npq)
sub skewness 	{		($_[0]->{_q} - $_[0]->{_p}) / sqrt($_[0]->{_n} * $_[0]->{_p} * $_[0]->{_q});		}

# Kurtosis == (1 - 6pq) / (Npq)
sub kurtosis	{		(1 - 6*($_[0]->{_p} * $_[0]->{_q})) / ($_[0]->{_n} * $_[0]->{_p} * $_[0]->{_q});		}

# Fisher information for a fixed number of trials N, and the probability of success, p.
# Gn(p) = n/pq
sub fisher_information		{		$_[0]->{_n} / ($_[0]->{_p} * $_[0]->{_q});	}

############################################################
#     			    Distribution Functions    		 	   #
############################################################

# Unlike the above subroutines, here we are using explicit variable names in the more complex equations.

# Probability mass function --> the probability of getting exactly k successes, given N trials and a probability of success p.
# 'pmf' and 'dbinom' are aliases here, they does the exact same thing.
sub binomial	{
	my ( $self, $k, $log ) = @_;
	my $n = $self->{_n};
	my $p = $self->{_p};
	my $q = 1 - $p;
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $binom = $self->_factorial($n) / $self->_factorial($k) / $self->_factorial($n-$k) * ( $p ** $k ) * ($q ** ($n - $k));
	( $log && int($log) == 1 )? return log($binom) : return $binom;
}
sub pmf	{
	my ( $self, $k, $log ) = @_;
	my $n = $self->{_n};
	my $p = $self->{_p};
	my $q = 1 - $p;
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $binom = $self->_factorial($n) / $self->_factorial($k) / $self->_factorial($n-$k) * ( $p ** $k ) * ($q ** ($n - $k));
	( $log && int($log) == 1 )? return log($binom) : return $binom;
}
sub dbinom	{
	my ( $self, $k, $log ) = @_;
	my $n = $self->{_n};
	my $p = $self->{_p};
	my $q = 1 - $p;
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $binom = $self->_factorial($n) / $self->_factorial($k) / $self->_factorial($n-$k) * ( $p ** $k ) * ($q ** ($n - $k));
	( $log && int($log) == 1 )? return log($binom) : return $binom;
}
sub probability	{
	my ( $self, $k, $log ) = @_;
	my $n = $self->{_n};
	my $p = $self->{_p};
	my $q = 1 - $p;
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $binom = $self->_factorial($n) / $self->_factorial($k) / $self->_factorial($n-$k) * ( $p ** $k ) * ($q ** ($n - $k));
	( $log && int($log) == 1 )? return log($binom) : return $binom;
}
sub binom	{
	my ( $self, $k, $log ) = @_;
	my $n = $self->{_n};
	my $p = $self->{_p};
	my $q = 1 - $p;
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $binom = $self->_factorial($n) / $self->_factorial($k) / $self->_factorial($n-$k) * ( $p ** $k ) * ($q ** ($n - $k));
	( $log && int($log) == 1 )? return log($binom) : return $binom;
}

# The binomial coefficient == (N choose K) --> the number of possible ways (combinations)
# of achieving a given numbers of successes, given the probability of success p and a number of possible outcomes N.
sub coefficient	{
	my ( $self, $k ) = @_;
	my $n = $self->{_n};
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) );
	return $self->_choose($n, $k);
}

# Cumulative distribution (density) function, CDF --> the sum of consecutive probabilities of a range of outcomes
# Calculate the cumulative distribution function given the number of independent trials, k.
# CDF == Pr(X <= k) == SUM(i=0, floor(k)) of (N choose i)p**i(1-p)**(n-i)
# Eg. Pr( Y <= 2 ) = pbinom(2,5,0.1) = 0.99144
# 'pbinom' is an alias for this function
sub cdf		{
	my ( $self, $k, $lower_tail, $log ) = @_;
	my $n = $self->{_n};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $cmf = 0.0;
	$p = 1 - $p if ( $lower_tail && int($lower_tail) == 1 );
	my $q = 1 - $p;
	for ( my $i=0; $i <= $k; $i++ )		{
		$cmf += $self->_choose($n, $i) * ( $p**$i ) * ($q**($n-$i));
	}
	( $log && int($log) == 1 )? return log($cmf) : return $cmf;
}
sub pbinom	{
	my ( $self, $k, $lower_tail, $log ) = @_;
	my $n = $self->{_n};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $cmf = 0.0;
	$p = 1 - $p if ( $lower_tail && int($lower_tail) == 1 );
	my $q = 1 - $p;
	for ( my $i=0; $i <= $k; $i++ )		{
		$cmf += $self->_choose($n, $i) * ( $p**$i ) * ($q**($n-$i));
	}
	( $log && int($log) == 1 )? return log($cmf) : return $cmf;
}
sub cumulative_distribution		{
	my ( $self, $k, $lower_tail, $log ) = @_;
	my $n = $self->{_n};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $cmf = 0.0;
	$p = 1 - $p if ( $lower_tail && int($lower_tail) == 1 );
	my $q = 1 - $p;
	for ( my $i=0; $i <= $k; $i++ )		{
		$cmf += $self->_choose($n, $i) * ( $p**$i ) * ($q**($n-$i));
	}
	( $log && int($log) == 1 )? return log($cmf) : return $cmf;
}

########################################
#  		 ESTIMATE P FROM A SAMPLE      #
########################################

# Maximum likelihood estimation of P, given a sample size N and the number of successes (and alias).
# We assume that the list of samples is encoded as 1 for a success and 0 for a failure.
sub binomial_mle	{
	my @sample = @{$_[1]};
	my $bigsum;
	foreach ( @sample )	{
		$bigsum += $_ if ( $_ == 1 );
	}
	return $bigsum / scalar @sample;
}
sub estimate_p_from_sample		{
	my @sample = @{$_[1]};
	my $bigsum;
	foreach ( @sample )	{
		$bigsum += $_ if ( $_ == 1 );
	}
	return $bigsum / scalar @sample;
}


############################################################
#     			    Sampling Functions    	   		 	   #
############################################################

# Here we want to randomly sample from a binomial distribution, which can
# be done by first computing the probability of observing each possible outcome,
# followed by using a pseudo-random number generator between 0 and 1.
# 'rbinom' is an alias here, just to be similar to R functionality.
# The variables here are slightly different:
# N possible outcomes, K trials, Probability of success P
sub sample	{
	my ( $self, $k, ) = @_;
	return unless ( $k >= 0 && $k == int $k );
	# We return an index to a possible outcome (ie. a slot in the Galton board) by choosing the number of RIGHTWARD moves
	my @outcomes;
	for ( my $i=0; $i < $k; $i++ )	{
		my $slot = 0;
		foreach ( my $j=0; $j < $self->{_n}; $j++ )		{
			my $rand = rand(1);
			$slot++ if ( $rand >= $self->{_p} );
		}
		push @outcomes, $slot;	
	}
	( scalar @outcomes == 1 )? return $outcomes[0] : return \@outcomes;
}

# Randomly sample the Binomial distribution K times
# Think of this as a Galton board with N pegs, K balls(trials), and a leftward probability P.
# This is just an alias for 'sample'
sub rbinom		{
	my ( $self, $k, ) = @_;
	return unless ( $k >= 0 && $k == int $k );
	# We return an index to a possible outcome (ie. a slot in the Galton board) by choosing the number of RIGHTWARD moves
	my @outcomes;
	for ( my $i=0; $i < $k; $i++ )	{
		my $slot = 0;
		foreach ( my $j=0; $j < $self->{_n}; $j++ )		{
			my $rand = rand(1);
			$slot++ if ( $rand >= $self->{_p} );
		}
		push @outcomes, $slot;	
	}
	( scalar @outcomes == 1 )? return $outcomes[0] : return \@outcomes;
}

####################################
#    FACTORIALS AND COMBINATIONS   #
####################################

# Computes the factorial of an integer value.
# Use the Gamma function for non-integers!
sub _factorial	{
	my ( $self, $n ) = @_;
	my $res = 1;
	return undef unless $n >= 0 and $n == int($n);		# Note: Non-integers require the Gamma function!
	$res *= $n-- while ( $n > 1 );
	return $res;
}

# choose(n,k) is the number of ways to choose k elements from a set of n elements,
# when the order of selection is irrelevant.
# Given by the formula for n choose k =
#   (n - k)!
# -------------
#  k!(n - k)!
sub _choose	{
	my ( $self, $n, $k ) = @_;
	my ( $result, $j ) = ( 1, 1 );
	
	return 0 if $k > $n || $k < 0;
	$k = ($n-$k) if ( $n-$k ) < $k;
	
	while ( $j <= $k )	{
		$result *= $n--;
		$result /= $j++;
	}
	return $result;
}

# Normalizes a list of numerical values	-- here probabilities	
sub _normalize	{
	my ( $self, @list ) = @_;
	my $probabilities;
	my $total = sum @list;	
	foreach my $probability ( @list )	{
		push @{$probabilities}, ($probability/$total);
	}
	return $probabilities;
}
####################################



# Probability of success here: 1.0 :-)
1;