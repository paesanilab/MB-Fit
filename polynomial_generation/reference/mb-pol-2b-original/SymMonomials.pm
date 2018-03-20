# (c) vbabin@zoo.astana.kz 01/2013

use strict;
use warnings;

#------------------------------------------------------------------------------#

BEGIN {
    use Carp;
    use Exporter ();

    our (@ISA, @EXPORT);
    @ISA = qw(Exporter);
    @EXPORT = qw(&integer_permutations &cartesian_product
                 &expand_molecule &make_permutations);
}

#------------------------------------------------------------------------------#

sub expand_molecule {
    my $short_name = shift;

    croak "bad molecule '$short_name'"
        unless $short_name =~ /^([a-zA-Z]+\d+)+$/;

    my $atoms = [];
    while ($short_name =~ /([a-zA-Z]+)(\d+)/g) {
        for (my $n = 0; $n < $2; ++$n) {
            push @{$atoms}, "$1$n";
        }

    }

    return $atoms;
}

#------------------------------------------------------------------------------#

sub integer_permutations {

    my @X = @_;
    my $nX = scalar(@X);

    my $result = [];

    {
        my @Xcopy = @X;
        push @{$result}, \@Xcopy;
    }

    while ($nX > 1) {

        # 1. Find the largest index k such that X[k] < X[k + 1].
        #    If no such index exists, the permutation is the last
        #    permutation.

        my $k = $nX - 2;

        while ($k >= 0) {
            last if $X[$k] < $X[$k + 1];
            --$k;
        }

        last if $k < 0;

        # 2. Find the largest index l such that X[k] < X[l]. Since k + 1
        #    is such an index, l is well defined and satisfies k < l.

        my $l = $nX - 1;

        while ($l > $k) {
            last if $X[$k] < $X[$l];
            --$l;
        }

        die unless $l > $k;

        # 3. Swap X[k] with X[l].

        {
            my $tmp = $X[$k];
            $X[$k] = $X[$l];
            $X[$l] = $tmp;
        }

        # 4. Reverse the sequence from X[k + 1] up to and including
        #    the final element X[n].

        {
            ++$k;
            my @tmp = @X[$k .. $nX - 1];
            while ($k < $nX) {
                $X[$k] = pop @tmp;
                ++$k;
            }
        }

        my @Xcopy = @X;
        push @{$result}, \@Xcopy;
    }

    return $result;
}

#------------------------------------------------------------------------------#

# stollen from Math::Cartesian::Product

sub cartesian_product {

    my @C = @_; # Lists to be multiplied
    my @c = (); # Current element of cartesian product
    my @P = (); # Cartesian product
    my $n = 0;  # Number of elements in product

    return 0 if @C == 0; # Empty product

    @C == grep {ref eq 'ARRAY'} @C
      or croak("arrays of things required by _cartesian");

    my $p;

    $p = sub {
        if (@c < @C) {
            for (@{$C[@c]}) {
                push @c, $_;
                &$p();
                pop @c;
            }
        } else {
            my $p = [@c];
            push @P, $p;
        }
    };

    &$p();

    return \@P;
}

#------------------------------------------------------------------------------#

#
# returns a list of the permutations of the like atoms
#

sub make_permutations {
    my $short_name = shift;

    croak "bad molecule '$short_name'"
        unless $short_name =~ /^([a-zA-Z]+\d+)+$/;

    my $idx = 0;
    my @perm1 = ();
    while ($short_name =~ /([a-zA-Z]+)(\d+)/g) {
        croak "bad molecule '$short_name'" unless $2 > 0;
        push @perm1, integer_permutations($idx .. $idx + $2 - 1);
        $idx += $2;
    }

    my $perms = cartesian_product(@perm1);

    my @result = ();

    foreach my $cc (@{$perms}) {
        my @curr = ();
        foreach (@{$cc}) {
            push @curr, @{$_};
        }
        push @result, \@curr;
    }

    return \@result;
}

#------------------------------------------------------------------------------#

1;

#------------------------------------------------------------------------------#
