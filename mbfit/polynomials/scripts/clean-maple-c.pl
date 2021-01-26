#!/usr/bin/perl

use strict;
use warnings;

my @entries = ();
my @this_entry = ();

# pre-process the input file. Maple inserts linebreaks between variable 
# definitions if they are > 80 characters. To keep from breaking this when
# renaming, we will read each variable into the entries array. Each element
# in entries is an array reference corresponding to a 'this_entry' array. 
# Each @this_entry element is a line from the maple-generated file.
#
# The size of @this_entry corresponds to the number of lines occupied by that
# variable in the maple-generated file

while (<STDIN>){
    if ((not $_ =~ m/{/) && (not $_ =~ m/}/) && (not $_ =~ m/double/)
      &&(not $_ =~ m/void/)){
        my $line = $_;
        chomp ($line);
        $line =~ s/\s+$//;

        push @this_entry, $line;

        #if this line ends in a semicolon, it terminates the entry
        if($line =~ m/;$/){
            my @arr = @this_entry;
            push @entries, \@arr;
            
            @this_entry = ();
        }
    }
}

# check to see if this code calculates gradients

my $doing_gradients = 0;
foreach my $r_entry (@entries){ # loop over every variable
    my $first_line = $r_entry->[0];
    if($first_line =~ m/crea_par/){
        $doing_gradients = 1;
    }
}

open (POLY_H, "poly-model.h") or die "poly-model.h must be in this directory\n";
my $n_vars = 0;
my $size = 0;
while(<POLY_H>){
    if($_ =~ m/n_vars = (\d+)/){
        $n_vars = $1;
    }
    if($_ =~ m/size = (\d+)/){
        $size = $1;
    }
}

print "#include \"poly-model.h\"

namespace mb_system {
\n";
if($doing_gradients == 1){
    print "double poly_model::eval(const double a[$size], const double x[$n_vars],
                        double g[$n_vars])\n";
}else{
    print "double poly_model::eval(const double a[$size], const double x[$n_vars])\n";
}
print "{\n";

# with the file repackaged into @entries, loop through each variable 
# and reformat it for c++ code

my $energy;
my $count = 0;
my @grd = ();

foreach my $r_entry (@entries){ # loop over every variable
    my $first_line = $r_entry->[0];

    if($first_line =~ m/^\s*t(\d+)\s*=/){
        $first_line =~ s/^\s*t(\d+)\s*=/    const double t$1 =/g;
        print $first_line, "\n";

        for(my $l = 1; $l < @{$r_entry}; ++$l){
            print $r_entry->[$l], "\n";
        }

    }elsif($first_line =~ m/^\s+crea_par/){
        # first crea_par[0] is the energy. all others are gradients
        my $line = $first_line;
        $line =~ s/^\s+//;
        chomp($line);
        my @toks = split(/\s+/, $line);

        # if the value spans multiple lines, add them here
        my $value = $toks[2];
        for(my $l = 1; $l < @{$r_entry}; ++$l){
            $value = $value."\n".$r_entry->[$l];
        }

        if($count == 0){
            $energy = $value;
        }else{
            push @grd, $value;
        }
        ++$count;
    }elsif($first_line =~ m/^\s+return/){
        if(defined $energy){
            for(my $i = 0; $i < scalar @grd; ++$i){
                print "    g[$i] = ", $grd[$i], "\n";
            }
            print "    return $energy\n";
        }else{
            print $first_line;
            for(my $l = 1; $l < @{$r_entry}; ++$l){
                print $r_entry->[$l], "\n";
            }
        }
    }
}

print "
}

} // namespace mb_system";
