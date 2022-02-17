#!/usr/bin/perl

# package NegativeBinomial.pm
# An HS custom Perl package for Negative Binomial distributions, properties, and sampling methods.

# Author: MH Seabolt
# Last Updated: 6-10-2021

package NegativeBinomial; 

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
		_r => '',						# A positive integer for the number of failures until the experiment is stopped.
		_p => '',  						# The probability of a success, given as a positive real number between 0 and 1 inclusive.
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

	# Check that we have parameter p and that falls within the range of 0 to 1 inclusive.
	return if ( $self->{_p} < 0 && $self->{_p} > 1.0 ); 

	# r must be a positive integer ge 0
	return if ( $self->{_r} < 0 || $self->{_r} != int $self->{_r} );
	
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
sub get_r	{	return $_[0]->{_r};	}
sub r 		{	return $_[0]->{_r};	}
sub set_r	{	
	return if ( $_[1] < 0 || $_[1] != int $_[1] );
	$_[0]->{_r} = $_[1];
}

################################################################################## 

############################################################
#         Properties of the Binomial Distribution      	   #
############################################################

# Some readability has been sacrificed here in favor of computational speed, but in general the equations are given and the order of args is the same that they appear in the given equation.

# Expected value (the expected mean of the distribution, given our parameters)
# I've included a few aliases here, but be very careful that we don't accidentally run into name collisions.
# Expected value (mean) 
sub expected	{		($_[0]->{_r} * $_[0]->{_p})/(1 - $_[0]->{_p});		}		
sub mean 		{		($_[0]->{_r} * $_[0]->{_p})/(1 - $_[0]->{_p});		}	

# Variance
sub variance	{	($_[0]->{_r} * $_[0]->{_p})/(1 - $_[0]->{_p})**2;		}

# Mode
sub mode 		{  		
	return 0 if ( $_[0]->{_r} >= 1 );
	return int(($_[0]->{_r} - 1) * $_[0]->{_p})/(1 - $_[0]->{_p});
} 

# Skewness
sub skewness 	{	(1 + $_[0]->{_p}) / sqrt($_[0]->{_p} * $_[0]->{_r});		}

# Kurtosis
sub kurtosis	{	(6/$_[0]->{_r}) + ((1 - $_[0]->{_p})**2 / ($_[0]->{_r} * $_[0]->{_p}));		}

# Fisher information for a fixed number of trials N, and the probability of success, p.
sub fisher_information		{		$_[0]->{_r} / ($_[0]->{_p}*((1-$_[0]->{_p})**2));	}

############################################################
#     			    Distribution Functions    		 	   #
############################################################

# Unlike the above subroutines, here we are using explicit variable names in the more complex equations.

# Probability mass function --> where k is the number of successes, r is the number of failures, and p is the probability of success. 
# 'pmf' and 'dnegbinom' are aliases here, they does the exact same thing.
sub negative_binomial	{
	my ( $self, $k, $log ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $r > 0 && $r == int($r) && $p > 0 && $p < 1 );
	my $negbinomial = $self->_choose(($k+$r-1),($r-1))*((1-$p)**$k)*($p**$r);
	( $log && int($log) == 1 )? return log($negbinomial) : return $negbinomial;
}
sub pmf	{
	my ( $self, $k, $log ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $r > 0 && $r == int($r) && $p > 0 && $p < 1 );
	my $negbinomial = $self->_choose(($k+$r-1),($r-1))*((1-$p)**$k)*($p**$r);
	( $log && int($log) == 1 )? return log($negbinomial) : return $negbinomial;
}
sub dnegbinom	{
	my ( $self, $k, $log ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $r > 0 && $r == int($r) && $p > 0 && $p < 1 );
	my $negbinomial = $self->_choose(($k+$r-1),($r-1))*((1-$p)**$k)*($p**$r);
	( $log && int($log) == 1 )? return log($negbinomial) : return $negbinomial;
}
sub probability	{
	my ( $self, $k, $log ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $r > 0 && $r == int($r) && $p > 0 && $p < 1 );
	my $negbinomial = $self->_choose(($k+$r-1),($r-1))*((1-$p)**$k)*($p**$r);
	( $log && int($log) == 1 )? return log($negbinomial) : return $negbinomial;
}
sub negbinom	{
	my ( $self, $k, $log ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $r > 0 && $r == int($r) && $p > 0 && $p < 1 );
	my $negbinomial = $self->_choose(($k+$r-1),($r-1))*((1-$p)**$k)*($p**$r);
	( $log && int($log) == 1 )? return log($negbinomial) : return $negbinomial;
}

# The negative binomial coefficient == (N choose K) --> the number of possible ways (combinations)
# of achieving a given numbers of successes, given the probability of success p and a number of possible outcomes N.
sub coefficient	{
	my ( $self, $k ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	return unless ( $k >= 0 && $k == int($k) && $r > 0 && $r == int($r) && $p > 0 && $p < 1 );
	return $self->_choose(($k+$r-1),($r-1));
}

# Cumulative distribution (density) function, CDF --> the sum of consecutive probabilities of a range of outcomes
# Calculate the cumulative distribution function given the number of independent trials, k.
# 'pnegbinom' is an alias for this function
sub cdf		{
	my ( $self, $k, $log ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	my $n = $r + $k;
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $cmf = 0.0;
	my $q = 1 - $p;
	for ( my $i=0; $i <= $k; $i++ )		{
		$cmf += $self->_choose($n, $i) * ( $p**$i ) * ($q**($n-$i));
	}
	( $log && int($log) == 1 )? return log($cmf) : return $cmf;
}
sub pnegbinom	{
	my ( $self, $k, $log ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	my $n = $r + $k;
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $cmf = 0.0;
	my $q = 1 - $p;
	for ( my $i=0; $i <= $k; $i++ )		{
		$cmf += $self->_choose($n, $i) * ( $p**$i ) * ($q**($n-$i));
	}
	( $log && int($log) == 1 )? return log($cmf) : return $cmf;
}
sub cumulative_distribution		{
	my ( $self, $k, $log ) = @_;
	my $r = $self->{_r};
	my $p = $self->{_p};
	my $n = $r + $k;
	return unless ( $k >= 0 && $k == int($k) && $n > 0 && $n == int($n) && $p > 0 && $p < 1 );
	my $cmf = 0.0;
	my $q = 1 - $p;
	for ( my $i=0; $i <= $k; $i++ )		{
		$cmf += $self->_choose($n, $i) * ( $p**$i ) * ($q**($n-$i));
	}
	( $log && int($log) == 1 )? return log($cmf) : return $cmf;
}

############################################################
#     			    Sampling Functions    	   		 	   #
############################################################

# Here we want to randomly sample from a negative binomial distribution, which can
# be done by first computing the probability of observing each possible outcome,
# followed by using a pseudo-random number generator between 0 and 1.
# 'rnegbinom' is an alias here, just to be similar to R functionality.
sub sample	{
	my ( $self, $k ) = @_;
	return unless ( $k >= 0 && $k == int $k );
	# We return an index to a possible outcome
	my @outcomes;
	for ( my $i=0; $i < $k; $i++ )	{
		my $slot = 0;
		foreach ( my $j=0; $j < ($self->{_r} + $k); $j++ )		{
			my $rand = rand(1);
			$slot++ if ( $rand >= $self->{_p} );
		}
		push @outcomes, $slot;	
	}
	( scalar @outcomes == 1 )? return $outcomes[0] : return \@outcomes;
}

# This is just an alias for 'sample'
sub rnegbinom		{
	my ( $self, $k, ) = @_;
	return unless ( $k >= 0 && $k == int $k );
	# We return an index to a possible outcome (ie. a slot in the Galton board) by choosing the number of RIGHTWARD moves
	my @outcomes;
	for ( my $i=0; $i < $k; $i++ )	{
		my $slot = 0;
		foreach ( my $j=0; $j < ($self->{_r} + $k); $j++ )		{
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