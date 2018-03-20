#!/usr/bin/perl

use strict;
use warnings;

use FindBin;
use lib $FindBin::Bin;

use SymMonomials;

my $monomer = 'O1H2';

my $atom_names = [];

push @{$atom_names}, (map {$_ . 'a';} @{&expand_molecule($monomer)});
push @{$atom_names}, (map {$_ . 'b';} @{&expand_molecule($monomer)});
push @{$atom_names}, (map {$_ . 'c';} @{&expand_molecule($monomer)});

# make it look more familiar (second character is also
# used to distinguish the molecules)

$atom_names = [map {
    $_ =~ s/O0/O/g;
    $_ =~ s/2(a|b|c)/${1}3/g;
    $_ =~ s/1(a|b|c)/${1}2/g;
    $_ =~ s/0(a|b|c)/${1}1/g;
    $_;} @{$atom_names}];

print  "\n<> atoms <>\n\n", join(':', @{$atom_names}), "\n";

my $natom = scalar(@{$atom_names});
my $natom_monomer = $natom/3;

my $atom_idx = {};
for (my $n = 0; $n < $natom; ++$n) {
    $atom_idx->{$atom_names->[$n]} = $n;
}

# build the permutation group

my $mp = make_permutations($monomer);
my $atom_permutations = [];

foreach my $p1 (@{$mp}) {
    foreach my $p2 (@{$mp}) {
        foreach my $p3 (@{$mp}) {
            my @xp1 = @{$p1};
            my @xp2 = (map {$_ + 1*$natom_monomer;} @{$p2});
            my @xp3 = (map {$_ + 2*$natom_monomer;} @{$p3});
            # swap the molecules explicitly
            push @{$atom_permutations}, [@xp1, @xp2, @xp3];
            push @{$atom_permutations}, [@xp1, @xp3, @xp2];
            push @{$atom_permutations}, [@xp2, @xp1, @xp3];
            push @{$atom_permutations}, [@xp2, @xp3, @xp1];
            push @{$atom_permutations}, [@xp3, @xp1, @xp2];
            push @{$atom_permutations}, [@xp3, @xp2, @xp1];
        }
    }
}

print "\n<> permutations (", scalar(@{$atom_permutations}), ") <>\n\n";
foreach (@{$atom_permutations}) {
    print join(':', @{$_}), "\n";
}

# atom, atom, class (class names should appear last once sorted)
my @variables = (['Ha1', 'Ha2', 'x-intra-HH'],
                 ['Hb1', 'Hb2', 'x-intra-HH'],
                 ['Hc1', 'Hc2', 'x-intra-HH'],
                 ['Oa',  'Ha1', 'x-intra-OH'],
                 ['Oa',  'Ha2', 'x-intra-OH'],
                 ['Ob',  'Hb1', 'x-intra-OH'],
                 ['Ob',  'Hb2', 'x-intra-OH'],
                 ['Oc',  'Hc1', 'x-intra-OH'],
                 ['Oc',  'Hc2', 'x-intra-OH'],
                 ['Ha1', 'Hb1', 'x-HH'],
                 ['Ha1', 'Hb2', 'x-HH'],
                 ['Ha1', 'Hc1', 'x-HH'],
                 ['Ha1', 'Hc2', 'x-HH'],
                 ['Ha2', 'Hb1', 'x-HH'],
                 ['Ha2', 'Hb2', 'x-HH'],
                 ['Ha2', 'Hc1', 'x-HH'],
                 ['Ha2', 'Hc2', 'x-HH'],
                 ['Hb1', 'Hc1', 'x-HH'],
                 ['Hb1', 'Hc2', 'x-HH'],
                 ['Hb2', 'Hc1', 'x-HH'],
                 ['Hb2', 'Hc2', 'x-HH'],
                 ['Oa',  'Hb1', 'x-OH'],
                 ['Oa',  'Hb2', 'x-OH'],
                 ['Oa',  'Hc1', 'x-OH'],
                 ['Oa',  'Hc2', 'x-OH'],
                 ['Ob',  'Ha1', 'x-OH'],
                 ['Ob',  'Ha2', 'x-OH'],
                 ['Ob',  'Hc1', 'x-OH'],
                 ['Ob',  'Hc2', 'x-OH'],
                 ['Oc',  'Ha1', 'x-OH'],
                 ['Oc',  'Ha2', 'x-OH'],
                 ['Oc',  'Hb1', 'x-OH'],
                 ['Oc',  'Hb2', 'x-OH'],
                 ['Oa',  'Ob',  'x-OO'],
                 ['Oa',  'Oc',  'x-OO'],
                 ['Ob',  'Oc',  'x-OO']);

my $nvar = scalar(@variables);

print "\n<> variables ($nvar) <>\n\n";
for (my $v = 0; $v < $nvar; ++$v) {
    printf "%2d : %3s <===> %3s : %s\n", $v,
        $variables[$v]->[0], $variables[$v]->[1], $variables[$v]->[2];
}

# map variables to integers
my $vtr_to_var = {}; # "vtr" = "variable's triple" (atom, atom, class)
for (my $v = 0; $v < $nvar; ++$v) {
    my $key = join(':', sort @{$variables[$v]});
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

# loop over 1st degree monomials

print "\n<> 1st degree <>\n\n";

my %mono_all;
for (my $k1 = 0; $k1 < $nvar; ++$k1) {
    my @mono = (0) x $nvar;
    $mono[$k1] += 1;
    $mono_all{join(':', @mono)} = undef;
}

print scalar keys %mono_all, " of all possible 1st degree monomials\n";
my $mono1_orb = find_orbits(\%mono_all);
die unless scalar(keys %mono_all) == 0;

print scalar @{$mono1_orb}, " <<== 1st degree 3B orbits\n";

# loop over 2nd degree monomials

print "\n<> 2nd degree <>\n\n";

for (my $k1 = 0; $k1 < $nvar; ++$k1) {
    for (my $k2 = 0; $k2 < $nvar; ++$k2) {
        my @mono = (0) x $nvar;
        $mono[$k1] += 1;
        $mono[$k2] += 1;
        $mono_all{join(':', @mono)} = undef;
    }
}

&filter_monomials(\%mono_all, 'x-intra-OH', 1);
&filter_monomials(\%mono_all, 'x-intra-HH', 1);
print scalar keys %mono_all, " of all possible (filtered) 2nd degree monomials\n";
my $mono2_orb = find_orbits(\%mono_all);
die unless scalar(keys %mono_all) == 0;

print  scalar @{$mono2_orb}, " <<== 2nd degree 3B orbits\n";

# loop over 3rd degree monomials

print "\n<> 3rd degree <>\n\n";

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

&filter_monomials(\%mono_all, 'x-intra-HH', 1);
&filter_monomials(\%mono_all, 'x-intra-OH', 1);
&filter_monomials(\%mono_all, 'x-OO', 3);
#&filter_monomials(\%mono_all, 'x-OH', 1);
&filter_monomials(\%mono_all, 'x-HH', 3);
print scalar keys %mono_all, " of all possible (filtered) 3rd degree monomials\n";
my $mono3_orb = find_orbits(\%mono_all);
die unless scalar(keys %mono_all) == 0;

print scalar @{$mono3_orb}, " <<== 3rd degree 3B orbits\n";

# loop over 4th degree monomials

print "\n<> 4th degree <>\n\n";

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

&filter_monomials(\%mono_all, 'x-intra-HH', 1);
&filter_monomials(\%mono_all, 'x-intra-OH', 1);
&filter_monomials(\%mono_all, 'x-OO', 1);
#&filter_monomials(\%mono_all, 'x-OH', 0);
&filter_monomials(\%mono_all, 'x-HH', 1);

print scalar keys %mono_all, " of all possible (filtered) 4th degree monomials\n";
my $mono4_orb = find_orbits(\%mono_all);
die unless scalar(keys %mono_all) == 0;

print scalar @{$mono4_orb}, " <<== 4th degree 3B orbits\n";

# print it out

print "\n<-> variables <->\n\n";
for (my $v = 0; $v < $nvar; ++$v) {
    print "    x[$v] = \@VAR\@-\|$variables[$v]->[2]\|",
        "($variables[$v]->[0], $variables[$v]->[1]);\n"
}

print "\n<-> monomials <->";

my $i0 = 0;
&print_orbits($i0, $mono1_orb);
$i0 += scalar(@{$mono1_orb});

print "\n\n";

#$i0 = 0;
&print_orbits($i0, $mono2_orb);
$i0 += scalar(@{$mono2_orb});

print "\n\n";

#$i0 = 0;
&print_orbits($i0, $mono3_orb);
$i0 += scalar(@{$mono3_orb});

print "\n\n";

#$i0 = 0;
&print_orbits($i0, $mono4_orb);
$i0 += scalar(@{$mono4_orb});

################################################################################

sub find_orbits {
    my $all = shift;

    my @orbits;

    my $maxtries = scalar(keys %{$all});

    while (1) {
        my @still_there = keys %{$all};
        last if scalar @still_there == 0 or $maxtries-- == 0;
        my $curr = shift @still_there;
        push @orbits, $curr if is_3B($curr);
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
        my ($a1, $a2) = ($atom_idx->{$vtr->[0]}, $atom_idx->{$vtr->[1]});
        my @permuted = ($atom_names->[$p->[$a1]], $atom_names->[$p->[$a2]]);
        push @permuted, $vtr->[2]; # klass stays the same
        my $permuted_vtr = join(':', sort @permuted);
        die " ** something is wrong **\n"
            unless exists($vtr_to_var->{$permuted_vtr});
        $dst[$vtr_to_var->{$permuted_vtr}] = $src[$v];
    }
    return join(':', @dst);
}

################################################################################

sub is_3B {
    my $mono_str = shift;
    my @mono = split(':', $mono_str);

    die unless $nvar == scalar @mono;

    my %mols;

    for (my $v = 0; $v < $nvar; ++$v) {
        next unless $mono[$v] > 0;
        my $var = $variables[$v];
        my $m1 = substr ($var->[0], 1, 1);
        my $m2 = substr ($var->[1], 1, 1);
        next if $m1 eq $m2; # intra-molecular variable => does not count
        $mols{$m1} = undef;
        $mols{$m2} = undef;
    }

    my $nmol = scalar keys %mols;
    return ($nmol == 3 ? 1 : 0);
}

################################################################################

sub filter_monomials {
    my ($monomials, $klass, $max_degree) = @_;

    my @todo = keys(%{$monomials});
    foreach my $m_str (@todo) {
        my $kd = get_klass_degree($m_str, $klass);
        delete ${$monomials}{$m_str} if $kd > $max_degree;
    }
}

################################################################################

sub filter_monomials_at_least {
    my ($monomials, $klass, $min_degree) = @_;

    my @todo = keys(%{$monomials});
    foreach my $m_str (@todo) {
        my $kd = get_klass_degree($m_str, $klass);
        delete ${$monomials}{$m_str} if $kd < $min_degree;
    }
}

################################################################################

sub get_klass_degree {
    my ($m_str, $klass) = @_;

    my @m = split(':', $m_str);
    my $nvar = scalar @m;

    my $klass_degree = 0;
    for (my $v = 0; $v < $nvar; ++$v) {
        next unless $m[$v] > 0;
        my $vtr = $variables[$v];
        $klass_degree += $m[$v] if ($vtr->[2] eq $klass);
    }

    return $klass_degree;
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
    print "    p[$i] = ", join(' + ', keys %terms), ";\n";
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
    print "    p[$i] := ", join('+', keys %terms), ":\n";
}

################################################################################

sub print_orbits {
    my ($i0, $orbs) = @_;

    foreach (@{$orbs}) {
        print_one($i0++, $_);
#        print_one_maple($i0++, $_);
    }
}

################################################################################
