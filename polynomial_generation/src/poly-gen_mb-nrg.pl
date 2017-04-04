#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;

use SymMonomials;

# check that command line parameters are defined
if(not @ARGV or @ARGV != 2 or $ARGV[0] !~ m/\d/){
    print STDERR "\nusage: poly-gen_mb-nrg.pl order input_file > log\n";
    &print_input_file_info;
    die ;
}

my $poly_order = shift;
my $input_file = shift;

#####
# Check for maximum polynomial order currently supported
#####
if ($poly_order > 6){
    die "Polynomial order > 6 not implemented.\n";
}

#####
# Open the output files for the polynomials
#####

#my $out_type = &read_output_type_from_input($input_file);

open (OUT_VARS, ">vars.cpp")
    or die "cannot open vars.cpp for write\n";
open (OUT_CPP, ">poly-direct.cpp")
    or die "cannot open poly-direct.cpp for write\n";
open (OUT_MAPLE_NOGRD, ">poly-nogrd.maple")
    or die "cannot open poly-nogrd.maple for write\n";
open (OUT_MAPLE_GRD, ">poly-grd.maple")
    or die "cannot open poly-grd.maple for write\n";

#####
# Read the molecule names that define the system from the input file (see usage)
# Each molecule name must be of the format:
#    atomName.numberOfAtoms.atomName.numberOfAtoms (no period)
# for example, O1H2 for water
#              H1O1Cl1 for HOCl
#####

my $r_molecules = &read_molecules_from_input($input_file);
my @molecules = @{$r_molecules};

# assume that the N-body order is the same as the number of molecules
my $NB_order = scalar@molecules;

print "<> molecules(", scalar@molecules, ") <>\n\n";
for(my $i = 0; $i < scalar@molecules; ++$i){
    print $molecules[$i], "\n";
}

if(scalar@molecules > 3){
    die "Currently interactions are only supported up to the 3B level\n";
}

# Parse the molecule names to determine and name the atoms in the system

my $atom_names = [];
my @natm_per_molec = ();
my @index_this_molec = ();
my $molname = 'a';
for(my $i = 0; $i < scalar@molecules; ++$i){
    push @{$atom_names}, 
                     (map {$_ . $molname;} @{&expand_molecule($molecules[$i])});
    my $ind = 0;
    for(my $j = 0; $j < scalar@natm_per_molec; ++$j){
        $ind += $natm_per_molec[$j];
    }
    push @index_this_molec, $ind;
    push @natm_per_molec, scalar @{&expand_molecule($molecules[$i])};
    ++$molname;
}

# remove number from oxygen
$atom_names = [map {$_ =~ s/-//; $_} @{$atom_names}];
# exchanges molecule index (a..z) with atom index (1,2,etc.)
$atom_names = [map {$_ =~ s/(\d)(\w)$/${2}${1}/; $_} @{$atom_names}];

print  "\n<> atoms (", scalar @{$atom_names}, 
              ") <>\n\n", join(':', @{$atom_names}), "\n";

my $natom = scalar(@{$atom_names});

# assigns an index to each atom name: e.g.,
#                            if $atom_names->[1] = 'Ha1', $atom_idx->{'Ha1'} = 1
my $atom_idx = {};
for (my $n = 0; $n < $natom; ++$n) {
    $atom_idx->{$atom_names->[$n]} = $n;
}

#####
# Now build the total permutation group
#####

# Step 1: need permutation group for each molecule
my @molec_permutations = ();
for(my $i = 0; $i < scalar@molecules; ++$i){
    my $mp = make_permutations($molecules[$i]);

    # shift all permutations by the atom_idx of the first atom of this molecule
    foreach (@{$mp}){
        @{$_} = (map {$_ + $index_this_molec[$i];} @{$_});
    }
    push @molec_permutations, $mp;
}

my $atom_permutations = [];

# Step 2: construct the total permutation group of the N-body system,
#         including permutation of like fragments
foreach my $p1 (@{$molec_permutations[0]}) {
    if(scalar@molecules > 1){
        foreach my $p2 (@{$molec_permutations[1]}) {
            if(scalar@molecules > 2){ # 3-body
                foreach my $p3 (@{$molec_permutations[2]}) {
                    my @xp1 = @{$p1};
                    my @xp2 = @{$p2};
                    my @xp3 = @{$p3};

                    push @{$atom_permutations}, [@xp1, @xp2, @xp3];
                    # if equivalent, swap fragments
                    if($molecules[0] eq $molecules[1]){# A == B
                        push @{$atom_permutations}, [@xp2, @xp1, @xp3];
                    }
                    if($molecules[0] eq $molecules[2]){# A == C
                        push @{$atom_permutations}, [@xp3, @xp2, @xp1];
                    }
                    if($molecules[1] eq $molecules[2]){# B == C
                        push @{$atom_permutations}, [@xp1, @xp3, @xp2];
                    }
                }
            }else{ # 2-body
                my @xp1 = @{$p1};
                my @xp2 = @{$p2};

                push @{$atom_permutations}, [@xp1, @xp2];
                # if equivalent, swap fragments
                if($molecules[0] eq $molecules[1]){# A == B
                    push @{$atom_permutations}, [@xp2, @xp1];
                }
            }
        }
    }else{ # 1-body
        my @xp1 = @{$p1};

        push @{$atom_permutations}, [@xp1];
    }
}

print "\n<> permutations (", scalar(@{$atom_permutations}), ") <>\n\n";
foreach (@{$atom_permutations}) {
    print join(':', @{$_}), "\n";
}

#####
# Read distance variables from input_file
#####

my $r_variables = &read_distance_variable_from_input($input_file);
my @variables = @{$r_variables};

my $nvar = scalar(@variables);

print "\n<> variables ($nvar) <>\n\n";
for (my $v = 0; $v < $nvar; ++$v) {
#    print $variables[$v]->{ATOMS}->[0]->{ID}, " ", 
#          $variables[$v]->{ATOMS}->[1]->{ID}, "\n"; 
    printf "%2d : %3s(%1s) <===> %3s(%1s) : %s\n", $v,
        $variables[$v]->{ATOMS}->[0]->{NAME},
	$variables[$v]->{ATOMS}->[0]->{MOLEC}, 
        $variables[$v]->{ATOMS}->[1]->{NAME},
	$variables[$v]->{ATOMS}->[1]->{MOLEC},
	$variables[$v]->{CLASS};
}

# map variables to integers
my $vtr_to_var = {}; # "vtr" = "variable's triple" (atom, atom, class)
for (my $v = 0; $v < $nvar; ++$v) {
    my @this_var = ($variables[$v]->{ATOMS}->[0]->{ID},
                    $variables[$v]->{ATOMS}->[1]->{ID},
                    $variables[$v]->{CLASS});
    my $key = join(':', sort @this_var);
    $vtr_to_var->{$key} = $v;
}

#
# monomial is a list of length $nvar of the powers of the variables;
# for example, (1, 2, 0, ..., 0) <=> x1*x2*x2; total degree is the
# sum of the elements
#
# permutation initially acts on the atoms and so is represented by a list
# of length $natom : (0, 1, 2, 3, ..., $natom) is "identity", and so on
#

my $i0 = 0; #counter for polynomial array index
my $tot_terms = 0; #should equal i0 at the end

# loop over 1st degree monomials
my $mono1_orb = undef;
if($poly_order > 0){

    print "\n<> 1st degree <>\n\n";

    my %mono_all;
    for (my $k1 = 0; $k1 < $nvar; ++$k1) {
        my @mono = (0) x $nvar;
        $mono[$k1] += 1;
        $mono_all{join(':', @mono)} = undef;
    }

    print scalar keys %mono_all, " possible 1st degree monomials\n";
    $mono1_orb = find_orbits(\%mono_all);
    die unless scalar(keys %mono_all) == 0;
    $tot_terms += scalar(@{$mono1_orb});

    print scalar @{$mono1_orb}, " <<== accepted 1st degree terms\n";

    print OUT_VARS "\n// <-> variables <->\n\n";
    for (my $v = 0; $v < $nvar; ++$v) {
        print OUT_VARS "    x[$v] = \@VAR\@-\|$variables[$v]->{CLASS}\|";

	my $r_atms = $variables[$v]->{ATOMS};

        print OUT_VARS "(";
	for (my $i = 0; $i < scalar@{$r_atms}; ++$i){
	    print OUT_VARS $r_atms->[$i]->{NAME},
	                   $r_atms->[$i]->{MOLEC};
	}
        print OUT_VARS ";\n";
    }

}

# loop over 2nd degree monomials
my $mono2_orb = undef;
if($poly_order > 1){

    print "\n<> 2nd degree <>\n\n";

    my %mono_all;
    for (my $k1 = 0; $k1 < $nvar; ++$k1) {
        for (my $k2 = 0; $k2 < $nvar; ++$k2) {
            my @mono = (0) x $nvar;
            $mono[$k1] += 1;
            $mono[$k2] += 1;
            $mono_all{join(':', @mono)} = undef;
        }
    }

    &check_for_filtered_monomials_and_do_it($input_file, \%mono_all, 2);

    print scalar keys %mono_all, " possible (filtered) 2nd degree monomials\n";
    $mono2_orb = find_orbits(\%mono_all);
    die unless scalar(keys %mono_all) == 0;
    $tot_terms += scalar(@{$mono2_orb});

    print  scalar @{$mono2_orb}, " <<== accepted 2nd degree terms\n";

}

# loop over 3rd degree monomials
my $mono3_orb = undef;
if($poly_order > 2){

    print "\n<> 3rd degree <>\n\n";

    my %mono_all;
    for (my $k1 = 0; $k1 < $nvar; ++$k1) {
        for (my $k2 = 0; $k2 < $nvar; ++$k2) {
            for (my $k3 = 0; $k3 < $nvar; ++$k3) {
                my @mono = (0) x $nvar;
                $mono[$k1] += 1;
                $mono[$k2] += 1;
                $mono[$k3] += 1;
                $mono_all{join(':', @mono)} = undef;
            }
        }
    }

    &check_for_filtered_monomials_and_do_it($input_file, \%mono_all, 3);

    print scalar keys %mono_all, " possible (filtered) 3rd degree monomials\n";
    $mono3_orb = find_orbits(\%mono_all);
    die unless scalar(keys %mono_all) == 0;
    $tot_terms += scalar(@{$mono3_orb});

    print scalar @{$mono3_orb}, " <<== accepted 3rd degree terms\n";

}

# loop over 4th degree monomials
my $mono4_orb = undef;
if($poly_order > 3){

    print "\n<> 4th degree <>\n\n";

    my %mono_all;
    for (my $k1 = 0; $k1 < $nvar; ++$k1) {
        for (my $k2 = 0; $k2 < $nvar; ++$k2) {
            for (my $k3 = 0; $k3 < $nvar; ++$k3) {
                for (my $k4 = 0; $k4 < $nvar; ++$k4) {
                    my @mono = (0) x $nvar;
                    $mono[$k1] += 1;
                    $mono[$k2] += 1;
                    $mono[$k3] += 1;
                    $mono[$k4] += 1;
                    $mono_all{join(':', @mono)} = undef;
                }
            }
        }
    }

    &check_for_filtered_monomials_and_do_it($input_file, \%mono_all, 4);

    print scalar keys %mono_all, " possible (filtered) 4th degree monomials\n";
    $mono4_orb = find_orbits(\%mono_all);
    die unless scalar(keys %mono_all) == 0;
    $tot_terms += scalar(@{$mono4_orb});

    print scalar @{$mono4_orb}, " <<== accepted 4th degree terms\n";

}

# loop over 5th degree monomials
my $mono5_orb = undef;
if($poly_order > 4){

    print "\n<> 5th degree <>\n\n";

    my %mono_all;
    for (my $k1 = 0; $k1 < $nvar; ++$k1) {
        for (my $k2 = 0; $k2 < $nvar; ++$k2) {
            for (my $k3 = 0; $k3 < $nvar; ++$k3) {
                for (my $k4 = 0; $k4 < $nvar; ++$k4) {
                    for (my $k5 = 0; $k5 < $nvar; ++$k5) {
                        my @mono = (0) x $nvar;
                        $mono[$k1] += 1;
                        $mono[$k2] += 1;
                        $mono[$k3] += 1;
                        $mono[$k4] += 1;
                        $mono[$k5] += 1;
                        $mono_all{join(':', @mono)} = undef;
                    }
                }
            }
        }
    }

    &check_for_filtered_monomials_and_do_it($input_file, \%mono_all, 5);

    print scalar keys %mono_all, " possible (filtered) 5th degree monomials\n";
    $mono5_orb = find_orbits(\%mono_all);
    die unless scalar(keys %mono_all) == 0;
    $tot_terms += scalar(@{$mono5_orb});

    print scalar @{$mono5_orb}, " <<== accepted 5th degree terms\n";

}

# loop over 6th degree monomials
my $mono6_orb = undef;
if($poly_order > 5){

    print "\n<> 6th degree <>\n\n";

    my %mono_all;
    for (my $k1 = 0; $k1 < $nvar; ++$k1) {
        for (my $k2 = 0; $k2 < $nvar; ++$k2) {
            for (my $k3 = 0; $k3 < $nvar; ++$k3) {
                for (my $k4 = 0; $k4 < $nvar; ++$k4) {
                    for (my $k5 = 0; $k5 < $nvar; ++$k5) {
                        for (my $k6 = 0; $k6 < $nvar; ++$k6) {
                            my @mono = (0) x $nvar;
                            $mono[$k1] += 1;
                            $mono[$k2] += 1;
                            $mono[$k3] += 1;
                            $mono[$k4] += 1;
                            $mono[$k5] += 1;
                            $mono[$k6] += 1;
                            $mono_all{join(':', @mono)} = undef;
                        }
                    }
                }
            }
        }
    }

    &check_for_filtered_monomials_and_do_it($input_file, \%mono_all, 6);

    print scalar keys %mono_all, " possible (filtered) 6th degree monomials\n";
    $mono6_orb = find_orbits(\%mono_all);
    die unless scalar(keys %mono_all) == 0;
    $tot_terms += scalar(@{$mono6_orb});

    print scalar @{$mono6_orb}, " <<== accepted 6th degree terms\n";

}

print "\n Total number of terms: $tot_terms\n";

&print_cpp_header($tot_terms, scalar@variables);
&print_cpp_direct_opening($tot_terms, scalar@variables);

# print the polynomials

if($poly_order > 0){
    &print_orbits($i0, $mono1_orb);
    $i0 += scalar(@{$mono1_orb});
    print OUT_CPP "\n";
    print OUT_MAPLE_NOGRD "\n";
    print OUT_MAPLE_GRD "\n";
}

if($poly_order > 1){
    &print_orbits($i0, $mono2_orb);
    $i0 += scalar(@{$mono2_orb});
    print OUT_CPP "\n";
    print OUT_MAPLE_NOGRD "\n";
    print OUT_MAPLE_GRD "\n";
}

if($poly_order > 2){
    &print_orbits($i0, $mono3_orb);
    $i0 += scalar(@{$mono3_orb});
    print OUT_CPP "\n";
    print OUT_MAPLE_NOGRD "\n";
    print OUT_MAPLE_GRD "\n";
}

if($poly_order > 3){
    &print_orbits($i0, $mono4_orb);
    $i0 += scalar(@{$mono4_orb});
    print OUT_CPP "\n";
    print OUT_MAPLE_NOGRD "\n";
    print OUT_MAPLE_GRD "\n";
}

if($poly_order > 4){
    &print_orbits($i0, $mono5_orb);
    $i0 += scalar(@{$mono5_orb});
    print OUT_CPP "\n";
    print OUT_MAPLE_NOGRD "\n";
    print OUT_MAPLE_GRD "\n";
}

if($poly_order > 5){
    &print_orbits($i0, $mono6_orb);
    $i0 += scalar(@{$mono6_orb});
    print OUT_CPP "\n";
    print OUT_MAPLE_NOGRD "\n";
    print OUT_MAPLE_GRD "\n";
}

die "Something is wrong: total number of term inconsistent\n"
    unless $i0 == $tot_terms;

# print the maple instructions
&print_maple_energy($i0, scalar@variables);
&print_maple_gradient($i0, scalar@variables);
&print_cpp_direct_closing($i0, scalar@variables);

close(OUT_VARS);
close(OUT_CPP);
close(OUT_MAPLE_NOGRD);
close(OUT_MAPLE_GRD);

################################################################################

sub find_orbits {
    my $all = shift;

    my @orbits;

    my $maxtries = scalar(keys %{$all});

    while (1) {
        my @still_there = keys %{$all};
        last if scalar @still_there == 0 or $maxtries-- == 0;
        my $curr = shift @still_there;
        push @orbits, $curr if is_NB($curr);
        for my $p (@{$atom_permutations}) {
            my $permuted = permute_monomial($curr, $p);
            delete ${$all}{$permuted};
        }
    }

    die "\n == something is wrong ==\n" if $maxtries < 0;

    return \@orbits;
}

################################################################################

sub permute_monomial {
    my ($src_str, $p) = @_;

    my @src = split(':', $src_str);
    my $nvar = scalar @src;
    my @dst = (0) x $nvar;

    for (my $v = 0; $v < $nvar; ++$v) {
        next unless $src[$v] > 0;

        my $vtr = $variables[$v];
	my $atm = $vtr->{ATOMS};
        my ($a1, $a2) = ($atom_idx->{$atm->[0]->{ID}},
	                 $atom_idx->{$atm->[1]->{ID}});
        my @permuted = ($atom_names->[$p->[$a1]], $atom_names->[$p->[$a2]]);
        push @permuted, $vtr->{CLASS}; # class stays the same
        #push @permuted, $vtr->[2]; # class stays the same
        my $permuted_vtr = join(':', sort @permuted);
        
        #       print "Permuted vtr: $permuted_vtr\n";

        if(not exists($vtr_to_var->{$permuted_vtr})){
            print STDERR "failed in permute_monomial: \n",
                         '     $vtr_to_var->{$permuted_vtr} does not exist\n',
                         '               $permuted_vtr  = ', $permuted_vtr, "\n",
                         ' $vtr_to_var->{$permuted_vtr} = ', $vtr_to_var->{$permuted_vtr}, "\n";
            die;
        }
        $dst[$vtr_to_var->{$permuted_vtr}] = $src[$v];
    }
    return join(':', @dst);
}

################################################################################

sub is_NB {
    my $mono_str = shift;
    my @mono = split(':', $mono_str);

    die unless $nvar == scalar @mono;

    my %mols;

    for (my $v = 0; $v < $nvar; ++$v) {

	#if we are testing this variable
        next unless $mono[$v] > 0;

        my $r_atms = $variables[$v]->{ATOMS};

        # if NB_order > 1, exclude purely intramolecular variables
	if($NB_order > 1){
    	    next if ($r_atms->[0]->{MOLEC} eq $r_atms->[1]->{MOLEC});
	}

	for(my $i = 0; $i < scalar@{$r_atms}; ++$i){
	    $mols{$r_atms->[$i]->{MOLEC}} = undef;
	}
    }

    my $nmol = scalar keys %mols;
    return ($nmol == $NB_order ? 1 : 0);
}

################################################################################

sub filter_monomials {
    my ($monomials, $class, $max_degree) = @_;

    my @todo = keys(%{$monomials});
    foreach my $m_str (@todo) {
        my $kd = get_class_degree($m_str, $class);
        delete ${$monomials}{$m_str} if $kd > $max_degree;
    }
}

################################################################################

sub filter_monomials_at_least {
    my ($monomials, $class, $min_degree) = @_;

    my @todo = keys(%{$monomials});
    foreach my $m_str (@todo) {
        my $kd = get_class_degree($m_str, $class);
        delete ${$monomials}{$m_str} if $kd < $min_degree;
    }
}

################################################################################

sub get_class_degree {
    my ($m_str, $class) = @_;

    my @m = split(':', $m_str);
    my $nvar = scalar @m;

    my $class_degree = 0;
    for (my $v = 0; $v < $nvar; ++$v) {
        next unless $m[$v] > 0;
        my $vtr = $variables[$v];
        $class_degree += $m[$v] if ($vtr->{CLASS} eq $class);
    }

    return $class_degree;
}

################################################################################

sub print_one {
    my ($i, $m) = @_;
    my %terms;
    for my $p (@{$atom_permutations}) {
        my $xx = permute_monomial($m, $p);
        my @yy = split(':', $xx);
        my @pp;
        for (my $n = 0; $n < scalar @yy; ++$n) {
            next unless $yy[$n] > 0;
            for (my $p = 0; $p < $yy[$n]; ++$p) {
                push @pp, "x[$n]";
            }
        }
        my $curr = join('*', sort @pp);
        $terms{$curr} = undef;
    }
    print OUT_CPP "    p[$i] = ", join(' + ', keys %terms), ";\n";
}

################################################################################

# maple cannot take derivatives w.r.t. arrays,
# therefore xNN are used in place of x[XX];
# p[] is indexed starting from 1; codegen[C] will
# change that back to the 0-based indexing

sub print_one_maple {
    my ($i, $m) = @_;
    my %terms;
    for my $p (@{$atom_permutations}) {
        my $xx = permute_monomial($m, $p);
        my @yy = split(':', $xx);
        my @pp;
        for (my $n = 0; $n < scalar @yy; ++$n) {
            next unless $yy[$n] > 0;
            for (my $p = 0; $p < $yy[$n]; ++$p) {
                my $nn = sprintf "x%02d", $n + 1;
                push @pp, $nn;
            }
        }
        my $curr = join('*', sort @pp);
        $terms{$curr} = undef;
    }
    $i++;
    print OUT_MAPLE_GRD   "    p[$i] := ", join('+', keys %terms), ":\n";
    print OUT_MAPLE_NOGRD "    p[$i] := ", join('+', keys %terms), ":\n";
}

################################################################################

sub print_orbits {
    my ($i0, $orbs) = @_;

    foreach (@{$orbs}) {
        print_one_maple($i0, $_);
        print_one($i0++, $_);
    }
}

################################################################################

sub read_distance_variable_from_input {
    my ($input_file) = @_;

    open (INP, $input_file) or die "Cannot open $input_file for read\n";

    my @variables = ();

    while(<INP>){
        if($_ =~ m/add_variable/i && $_ !~ m/^#/i && $_ !~ m/&!/i){
            my $line = $_;
            # remove junk, leaving space separated values
            $line =~ s/add_variable//i;
            $line =~ s/'//g;
            $line =~ s/\[//g;
            $line =~ s/\]//g;
            $line =~ s/;//g;
            $line =~ s/,//g;

            # remove whitespace from beginning and end
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            chomp($line);

            my @toks = split(/\s+/, $line);
            die "error reading distance variable from input, expected format:
            add_variable['H1', 'a', 'H2', 'a',  'x-intra-HH']\n"
                unless scalar @toks == 5;

	    my $variable = {};
	    my @atoms = ();
	    for(my $i = 0; $i < 2; ++$i){# for each atom
		my $atom = {};
		$atom->{NAME} = $toks[2*$i + 0]; # 'H1' for first atom
		$atom->{MOLEC} = $toks[2*$i + 1];# 'a'  for first atom

		my $id = $atom->{NAME}.$atom->{MOLEC};
		if($id =~ /(\w+)(\d+)(\w+)/){
		    $id = $1.$3.$2;
		}
		$atom->{ID} = $id;

		push @atoms, $atom;
	    }
	    $variable->{ATOMS} = \@atoms;
	    $variable->{CLASS} = $toks[4]; # 'x-intra-HH'
            push @variables, $variable;
        }
    }
    close(INP);
    return \@variables;
}

################################################################################

sub check_for_filtered_monomials_and_do_it {
    my ($input_file, $r_mono_all, $poly_order) = @_;

    open (INP, $input_file) or die "Cannot open $input_file for read\n";

    while(<INP>){
        if($_ =~ m/add_monomial_filter/i && $_ !~ m/^#/i && $_ !~ m/&!/i){
            my $line = $_;
            # remove junk, leaving space separated values
            $line =~ s/add_monomial_filter//i;
            $line =~ s/'//g;
            $line =~ s/\[//g;
            $line =~ s/\]//g;
            $line =~ s/;//g;
            $line =~ s/,//g;

            # remove whitespace from beginning and end
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            chomp($line);

            my @toks = split(/\s+/, $line);

            if( scalar@toks != 3){
                print STDERR "error reading monomial filter request from input, expected format:\n";
                print STDERR "  >>> add_monomial_filter[4, 'x-intra-HH', 1]\n";
                print STDERR " [ * 4 is polynomial order to filter, \n";
                print STDERR "   * 'x-intra-HH' is the class being filtered,\n";
                print STDERR "   * 1 is the maximum allowed order of this class at this poly. order]\n";

                die;
            }

            if($toks[0] == $poly_order){
                &filter_monomials($r_mono_all, $toks[1], $toks[2])
            }
        }
    }
    close(INP);
}

################################################################################

sub read_molecules_from_input{
    my ($input_file) = @_;

    open (INP, $input_file) or die "Cannot open $input_file for read\n";

    my @molecules = ();

    while(<INP>){
        if($_ =~ m/add_molecule/i && $_ !~ m/^#/i && $_ !~ m/&!/i){
            my $line = $_;
            # remove junk, leaving space separated values
            $line =~ s/add_molecule//i;
            $line =~ s/'//g;
            $line =~ s/\[//g;
            $line =~ s/\]//g;
            $line =~ s/;//g;
            $line =~ s/,//g;

            # remove whitespace from beginning and end
            $line =~ s/^\s+//;
            $line =~ s/\s+$//;
            chomp($line);

            my @toks = split(/\s+/, $line);

            if( scalar@toks != 1){
                print STDERR "error reading molecule, expected format:\n";
                print STDERR "  >>> add_molecule['O1H2']\n";
                die;
            }

            push @molecules, $toks[0];
        }
    }
    close(INP);
    return \@molecules;
}

################################################################################

sub read_output_type_from_input{
    my ($input_file) = @_;

    open (INP, $input_file) or die "Cannot open $input_file for read\n";

    while(<INP>){
        if($_ =~ m/maple/i && $_ !~ m/^#/i && $_ !~ m/&!/i){
            # output requested as maple input file
            if($_ =~ m/grad/i){
                return 2; # request e+grd
            }else{
                return 1; # request e
            }
        }
    }
    close(INP);
    return 0; # normal C output requested
}

################################################################################

sub print_input_file_info {

    print STDERR "\n >>>> Example options for input file\n\n";
    print STDERR "add_molecule['O1H2']\n\n";
    print STDERR "add_variable['H1', 'a', 'H2', 'a', 'x-intra-HH']\n\n";
    print STDERR "add_monomial_filter[4, 'x-intra-HH', 1]\n";
    print STDERR " >   4 is polynomial order to filter, \n";
    print STDERR " >   'x-intra-HH' is the class being filtered,\n";
    print STDERR " >   1 is the maximum allowed order of this class at this poly. order\n";
    print STDERR "\n\n";

}

################################################################################

sub print_maple_energy {

    my ($nterms, $nvar) = @_;

print OUT_MAPLE_NOGRD "\nenergy := 0;
for k from 1 by 1 to $nterms do
    energy := energy + a[k]*p[k]:
od:

args := [";

for(my $i = 1; $i < $nvar; ++$i){
    if(($i)%10 == 0){
        printf OUT_MAPLE_NOGRD "x%02d,\n         ", $i;
    }else{
        printf OUT_MAPLE_NOGRD "x%02d, ", $i;
    }
}
printf OUT_MAPLE_NOGRD "x%02d]:\n\n", $nvar;
print OUT_MAPLE_NOGRD "energy := convert(energy, 'horner', args):

energy_proc := codegen[makeproc](energy, parameters = args):
codegen[cost](energy_proc);

xxx := codegen[optimize](energy_proc):
codegen[cost](xxx);

xxx := codegen[packargs](xxx, args, x):
xxx := codegen[optimize](xxx):

codegen[C](xxx, optimized, filename=\"poly-nogrd.c\"):";

}

################################################################################

sub print_maple_gradient {

    my ($nterms, $nvar) = @_;

print OUT_MAPLE_GRD "\nenergy := 0;
for k from 1 by 1 to $nterms do
    energy := energy + a[k]*p[k]:
od:

args := [";

for(my $i = 1; $i < $nvar; ++$i){
    if(($i)%10 == 0){
        printf OUT_MAPLE_GRD "x%02d,\n         ", $i;
    }else{
        printf OUT_MAPLE_GRD "x%02d, ", $i;
    }
}
printf OUT_MAPLE_GRD "x%02d]:\n\n", $nvar;
print OUT_MAPLE_GRD "energy := convert(energy, 'horner', args):

energy_proc := codegen[makeproc](energy, parameters = args):
codegen[cost](energy_proc);

xxx := codegen[GRADIENT](energy_proc, args, function_value = true):
codegen[cost](xxx);
xxx := codegen[optimize](xxx):
codegen[cost](xxx);

xxx := codegen[packargs](xxx, args, x):
xxx := codegen[optimize](xxx):

codegen[C](xxx, optimized, filename=\"poly-grd.c\"):";

}

################################################################################

sub print_cpp_header {

    my ($nterms, $nvar) = @_;

    open (POLY_H, ">poly-model.h") or die "Cannot open poly-model.h for write\n";

    print POLY_H "#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = $nvar;
    static const unsigned size = $nterms;

    static double eval(const double a[$nterms],
                       const double x[$nvar]);

    static double eval(const double a[$nterms],
                       const double x[$nvar],
                             double g[$nvar]);

    static double eval_direct(const double a[$nterms],
                              const double x[$nvar]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H";

    close (POLY_H);
}

################################################################################

sub print_cpp_direct_opening {

    my ($nterms, $nvar) = @_;

    print OUT_CPP "#include \"poly-model.h\"

namespace mb_system {

double poly_model::eval_direct(const double a[$nterms], const double x[$nvar])
{
    double p[$nterms];\n"

}

################################################################################

sub print_cpp_direct_closing {

    my ($nterms, $nvar) = @_;

    print OUT_CPP "    double energy(0);
    for(int i = 0; i < $nterms; ++i)
        energy += p[i]*a[i];

    return energy;

}
} // namespace mb_system";

}
