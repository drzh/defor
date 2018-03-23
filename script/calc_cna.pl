#!/usr/bin/perl

# Usage: prog -r FILE_deprat -c FILE_nclust -a FILE_output_depratio

use strict;
use warnings;
use Getopt::Long;

my $filed = '';
my $filec = '';
my $filea = '';
my $para_d = 0.1;
my $para_u = 1;
my $para_b = 0.05;
my $para_l = 1000000;
# my $para_l = 0;
my $para_2l = 100000000;
my $para_meandiff = 0.15;
my $para_meanpairdiff = 0.05;
my $para_freqlow = 0.1;
my $para_freqhigh = 0.9;
my $para_meanlow = 0.4;
my $para_meanhigh = 0.65;
my $para_meanamplow = 0.25;
my $para_meanamphigh = 0.75;
my $para_fold = 3;
my $para_fold_loh = 10;
my $para_len = 10000000;
my $para_n = 10;
my $para_v = 0;
my $para_lcpn2 = 50000000;
my $para_gc = 1;
my $para_prop_ncl = 0.6;
my $para_n_ncl = 3;
my $para_gaplen = 5000000;
my $para_diff_loh_lohamp = 0.1;
my $para_maftrim = 0.45;
my $para_lgrstep = 0.01;
my $para_sdmin = 0.07;
my $para_cpn2lenmin = 50000000;
my $para_gclgrmin = -0.15;
my $para_gclgrmax = 0.15;
my $para_cpn2diffstrict = 0.1;
my $para_cpn2diffloose = 0.15;
my $para_p1diff = 0.05;
my $para_p3diff = 0.05;
my $para_uncdiff = 0.15;

Getopt::Long::Configure("bundling");
GetOptions(
           'r:s' => \$filed,
	   'c:s' => \$filec,
	   'a:s' => \$filea,
	   'd:f' => \$para_d,
	   'u:f' => \$para_u,
	   'b:f' => \$para_b,
	   'l:i' => \$para_l,
	   'y:i' => \$para_2l,
	   'D:f' => \$para_meandiff,
	   'f:f' => \$para_freqlow,
	   'F:f' => \$para_freqhigh,
	   'm:f' => \$para_meanlow,
	   'M:f' => \$para_meanhigh,
	   'n:i' => \$para_n,
	   # 'l:i' => \$para_len,
	   'v:i' => \$para_v
          );

my $lgr2 = 0;
my $lgr1 = log(1 / 2) / log(2);
my $lgr3 = log(3 / 2) / log(2);
my $lgr4 = log(4 / 2) / log(2);

my %statname = (
		'norm' => 'normal',
		'loh.amp' => 'loh-amp',
		'loh' => 'loh',
		# 'loh.norm' => 'loss,normal',
		'tri' => 'amp',
		'qua' => 'amp',
		# 'loh.tri' => 'loss,amp',
		'loss' => 'loss'
	       );

my @rat;
my @ratloo;
my %ratchridxstart;
my %ratchridxend;
my %ncl;
my %cpnlen;
my $flag = 1;
my $sstart = 0;
my $send = 0;
my $estart = 0;
my $eend = 0;
my $chr = 0;
my $flag_inter = 0;
my $cnpre = -1;
my $cn = -1;
my @meanpre = (0, 0, 0, 0);
my @res;
my $flag_pair = -1;

# Read nclust
open FILE, "<$filec"
  || die "Cannot open $filec: $!\n";
while (<FILE>) {
  chomp;
  my @e = split /\t/;
  push @{$ncl{$e[0]}}, [ @e ];
}
close FILE;

exit 0 if (! %ncl);

my %cpn2rat;
my %cpn2len;
my %cpn2end;
my %cpn2ratloo;
my %cpn2lenloo;
my %cpn2endloo;

my $idx = 0;
open FILE, "<$filed"
  || die "Cannot open $filed: $!\n";
while (<FILE>) {
  chomp;
  my @e = split /\t/;
  $e[5] = int($e[4] / $para_b) * $para_b + $para_b / 2;
  my @cpnstat = (0, 0, 0, 0);
  my @frq;
  if (exists $ncl{$e[0]}) {
    foreach my $r (@{$ncl{$e[0]}}) {
      my @c = @{$r};
      if ($e[2] >= $c[1]) {
	if ($e[1] <= $c[2] && $c[3] >= 0) {
	  # $cpnstat[$c[3]]++;
          $cpnstat[ ($c[4] <= $para_cpn2diffstrict) ? 2 : 1 ]++;
	  push @frq, $c[4];
	}
      }
      else {
	last;
      }
    }
  }
  $e[6] = join(',', @cpnstat[0 .. 3]);
  my $cpnstatsum = $cpnstat[0] + $cpnstat[1] + $cpnstat[2] + $cpnstat[3];
  my $cpn = 0;
  my @cpnall;
  foreach my $i (1 .. 3) {
    if ($cpnstat[$i] > $cpnstat[$cpn]) {
      $cpn = $i;
    }
    if ($cpnstat[$i] > 0) {
      push @cpnall, $i;
    }
  }
  if ($cpnstat[$cpn] == 0) {
    $cpn = -1;
  }
  elsif ($cpnstat[$cpn] < $para_n_ncl ||
	 $cpnstat[$cpn] / $cpnstatsum < $para_prop_ncl) {
    $cpn = join(',', @cpnall);
  }
  # push @rat, [ (@e, $cpn) ];
  $rat[$idx] = [ (@e, $cpn, (@frq > 0) ? join(';', @frq) : '.') ];
  if (exists $ratchridxstart{$e[0]}) {
    $ratchridxend{$e[0]} = $idx;
  }
  else {
    $ratchridxstart{$e[0]} = $idx;
  }
  # Select normal region accord strict critria
  my $cpnloo = $cpn;
  if ($cpn eq '2') {
    asigncpn2(\@e, \%cpn2rat, \%cpn2len, \%cpn2end);
  }
  # Select normal region accord loose critria
  if ($cpn eq '2' || $cpn eq '1' && median(@frq) < $para_cpn2diffloose) {
    asigncpn2(\@e, \%cpn2ratloo, \%cpn2lenloo, \%cpn2endloo);
    $cpnloo = 2;
  }
  $ratloo[$idx] = [ (@e, $cpnloo, (@frq > 0) ? join(';', @frq) : '.') ];
  $idx++;
}
close FILE;

# For each chr, calc mean and sd
my $chrminsd = '';
my $minsd = 100;
my $minmedsd = 100;
my %cpn2med;
my %cpn2mean;
my %cpn2sd;
my @medall;
my @sdall;
my @lenall;
my $flag_cpn2 = 0;
my @ratforgc = @rat;
while ($flag_cpn2++ < 2) {
  foreach my $chr (keys %cpn2rat) {
    ($cpn2mean{$chr}, $cpn2sd{$chr}) = meansd(@{$cpn2rat{$chr}});
    $cpn2med{$chr} = median(@{$cpn2rat{$chr}});
    if ($cpn2len{$chr} >= $para_lcpn2) {
      push @medall, $cpn2med{$chr};
      push @sdall, $cpn2sd{$chr};
      push @lenall, $cpn2len{$chr};
      if ($minmedsd > $cpn2med{$chr} + $cpn2sd{$chr}) {
        $chrminsd = $chr;
        $minmedsd = $cpn2med{$chr} + $cpn2sd{$chr};
      }
    }
  }
  if ($chrminsd eq '') {
    %cpn2rat = %cpn2ratloo;
    %cpn2len = %cpn2lenloo;
    %cpn2end = %cpn2endloo;
    @ratforgc = @ratloo;
    %cpn2med = ();
    %cpn2mean = ();
    %cpn2sd = ();
    @medall = ();
    @sdall = ();
    @lenall = ();
  }
  else {
    last;
  }
}
foreach my $chr (keys %cpn2rat) {
  my $rankmedsd = rank($cpn2mean{$chr}, \@medall) + rank($cpn2sd{$chr}, \@sdall) - rank($cpn2len{$chr}, \@lenall);
  if ($cpn2len{$chr} >= $para_lcpn2 && $minmedsd > $rankmedsd) {
    $chrminsd = $chr;
    $minmedsd = $rankmedsd;
  }
}

# test
print STDERR 'Stat before correction:', "\n";
foreach my $chr (sort keys %cpn2mean) {
  print STDERR join("\t", $chr, $cpn2len{$chr}, $cpn2med{$chr}, $cpn2mean{$chr}, $cpn2sd{$chr}, rank($cpn2med{$chr}, \@medall), rank($cpn2sd{$chr}, \@sdall), rank($cpn2len{$chr}, \@lenall)), "\n";
}
print STDERR 'chrmin: ', $chrminsd, "\n"; 
print STDERR '----------', "\n";

# Choose chr (>50MB) with minimal sd for GC correction
my %gcrat;
my %gcscale;

if ($para_gc == 1 && $chrminsd ne '') {
  # foreach my $r (@rat[$ratchridxstart{$chrminsd} .. $ratchridxend{$chrminsd}]) {
  foreach my $r (@ratforgc[$ratchridxstart{$chrminsd} .. $ratchridxend{$chrminsd}]) {
    my @e = @{$r};
    if ($e[7] eq '2') {
      push @{$gcrat{$e[5]}}, $e[3];
    }
  }
  gcscale(\%gcrat, \%gcscale);
}
applygcscale(\@rat, \%gcscale);

# Calc mean and sd
my $mean = 0;
my $med = 0;
my $lgrsd = $para_sdmin;
my $lgrscale = 0;
($lgrscale, $lgrsd) = cpn2stat(\@rat, \%ratchridxstart, \%ratchridxend);

# Define low and high limit for copy number
my ($low1, $high1, $low2, $high2, $low3, $high3, $low4, $high4, $lowloh, $highloh, $lowlohamp, $highlohamp, $lowlohtri, $highlohtri, $highlohnorm);
define_cutoff();

my @block;
my @uncer;
my $start = 0;
my $end = 0;
$chr = '';
my $stat = '';
my $idxstart = 0;
my $idxend = 0;
my $cpn = 0;
my @lgrs;
my @idx;

# Preliminary calling
callcn(\@rat);
# my ($sumstat, $lenlof, $lentri) = calcstat($lgrscale, \@rat);
# print STDERR join("\t", $stat, $lenlof, $lentri), "\n";

# # Refine GC correction
# %gcrat = ();
# %gcscale = ();
# foreach my $i (0 .. $#rat) {
#   my ($cchr, $cstart, $cend, $clgru, $gc, $gcg, $clu, $ccpn, $cfrq, $gcgc, $clgr, $cstat) = @{$rat[$i]};
#   if ($ccpn eq 2 && $clgr - $lgrscale > $para_gclgrmin && $clgr < $para_gclgrmax) {
#     push @{$gcrat{$gcg}}, $clgru;
#   }
# }
# gcscale(\%gcrat, \%gcscale);
# applygcscale(\@rat, \%gcscale);

# # Recalculate mean and sd
# ($lgrscale, $lgrsd) = cpn2stat(\@rat, \%ratchridxstart, \%ratchridxend);

# # Redefine cutoff
# define_cutoff();

# # Recall
# callcn(\@rat);

# # Test
# foreach my $r (@rat) {
#   print STDERR join("\t", @{$r}), "\n";
# }

# Continue calling
foreach my $i (0 .. $#rat) {
  my ($cchr, $cstart, $cend, $clgru, $gc, $gcg, $clu, $ccpn, $cfrq, $gcgc, $clgr, $cstat) = @{$rat[$i]};
  if ($cstat ne '') {
    if ($stat eq '') {
      $chr = $cchr;
      $start = $cstart;
      $end = $cend;
      $stat = $cstat;
      $idxstart = $i;
      $idxend = $i;
      $cpn = $ccpn;
      push @lgrs, $clgr;
      push @idx, $i;
    } else {
      if (($cstat eq $stat
           || ($cstat eq 'loh' || $cstat eq 'loh.amp' || $cstat eq 'loss')
	   && ($stat eq 'loh' || $stat eq 'loh.amp' || $cstat eq 'loss')
	   && $clgr - median(@lgrs) < $para_diff_loh_lohamp
	   && median(@lgrs) - $clgr < $para_diff_loh_lohamp
	  )
	  && $cchr eq $chr 
	  && $cstart < $end + $para_gaplen) {
	$end = $cend;
	$idxend = $i;
	push @lgrs, $clgr;
	push @idx, $i;
	if ($cstat eq 'loh') {
	  $stat = 'loh';
	  $cpn = 1
	}
      } else {
	push @block, [ ($chr, $start, $end, $stat, $cpn, $idxstart, $idxend, [ @idx ]) ];
	@lgrs = ();
	@idx = ();
	$chr = $cchr;
	$start = $cstart;
	$end = $cend;
	$stat = $cstat;
	$idxstart = $i;
	$idxend = $i;
	$cpn = $ccpn;
	push @lgrs, $clgr;
	push @idx, $i;
      }
    }
  } else {
    push @uncer, $i;
  }
}
if ($stat ne '') {
  push @block, [ ($chr, $start, $end, $stat, $cpn, $idxstart, $idxend, [ @idx ]) ];
}

# foreach my $chr (sort keys %cpn2len) {
#   print join("\t", $chr, $cpn2len{$chr}, $cpn2med{$chr}, $cpn2mean{$chr}, $cpn2sd{$chr}), "\n";
# }
# print join("\t", 'chrminsd', $chrminsd), "\n";
# foreach my $gc (sort {$a <=> $b} keys %gcscale) {
#   print $gc, "\t", $gcscale{$gc}, "\n";
# }
# print join("\t", 'chrminsd', $chrminsd), "\n";
# foreach my $r (@rat) {
#   print join("\t", @{$r}), "\n";
# }

# Process uncertained segments
if ($para_u == 1) {
  my @uncer2;
  foreach my $iu (@uncer) {
    my @u = @{$rat[$iu]};
    my $flag = 1;
    foreach my $rb (@block) {
      my @b = @{$rb};
      if ($u[0] eq $b[0] && $u[1] < $b[2] && $u[2] > $b[2]) {
      	my @r = @{$rat[$b[6]]};
      	my $lgr = $u[10];
        my $blocklgr = ratlgrmedian(@rat[@{$b[7]}]);
      	# if ($b[4] eq '2' && $lgr > $low2 && $lgr < $high2 ||
      	#     $b[4] eq '1' && $lgr < $high1 ||
      	#     $b[4] eq '3' && $lgr > $low3 ||
      	#     $b[4] eq '4' && $lgr > $low4 ||
	#     $b[4] eq '0' && $lgr < $high1
      	#    ) {
      	if (abs($lgr - $blocklgr) < $para_uncdiff) {
          $rb -> [2] = $u[2];
      	  $rb -> [6] = $iu;
      	  push @{$rb -> [7]}, $iu;
      	  $flag = 0;
      	  last;
      	}
      }
    }
    if ($flag == 1) {
      unshift @uncer2, $iu;
    }
  }
  
  foreach my $iu (@uncer2) {
    my @u = @{$rat[$iu]};
    foreach my $rb (reverse @block) {
      my @b = @{$rb};
      if ($u[0] eq $b[0] && $u[1] < $b[1] && $u[2] > $b[1]) {
  	my @r = @{$rat[$b[5]]};
  	# if ($b[4] eq '2' && $u[10] > $low2 && $u[10] <= $high2 ||
  	#     $b[4] eq '1' && $u[10] < $high1 ||
  	#     ($b[4]eq '3' || $b[4] eq '4') && $u[10] > $low3 ||
	#     $b[4]
  	#    ) {
      	my $lgr = $u[10];
        my $blocklgr = ratlgrmedian(@rat[@{$b[7]}]);
      	# if ($b[4] eq '2' && $lgr > $low2 && $lgr < $high2 ||
      	#     $b[4] eq '1' && $lgr < $high1 ||
      	#     $b[4] eq '3' && $lgr > $low3 ||
      	#     $b[4] eq '4' && $lgr > $low4 ||
	#     $b[4] eq '0' && $lgr < $high1
      	#    ) {
      	if (abs($lgr - $blocklgr) < $para_uncdiff) {
  	  $rb -> [1] = $u[1];
  	  $rb -> [5] = $iu;
  	  push @{$rb -> [7]}, $iu;
  	  last;
  	}
      }
    }
  }
}

# Process overlap
foreach my $i (0 .. ($#block - 1)) {
  if ($block[$i][0] eq $block[$i + 1][0] && $block[$i][2] > $block[$i + 1][1]) {
    $block[$i][2] = int(($block[$i][2] + $block[$i + 1][1]) / 2);
    $block[$i + 1][1] = $block[$i][2] + 1;
  }
}

# Merge
@block = merge(@block);

# Drop short segments
my @idxstart = (0);
foreach my $i (1 .. $#block) {
  next if ($block[$i][2] - $block[$i][1] + 1 >= $para_l);
  if ($i == $#block || $block[$i][0] ne $block[$i + 1][0]) {
    # Processing last segment for each chr
    if ($block[$i][0] ne $block[$i - 1][0] || 
	$block[$i][3] ne $block[$i - 1][3]
       ) {
      $block[$i][3] = ($block[$i][3] eq $block[$i - 1][3]) ? $block[$i][3] : 'NA';
    }
  }
  elsif ($block[$i][0] ne $block[$i - 1][0]) {
    # Processing first segment for each chr
    push @idxstart, $i;
  }
  else {
    # Processing internal segment
    if ($block[$i][3] ne $block[$i - 1][3] ||
	$block[$i][3] ne $block[$i + 1][3]
       ) {
      $block[$i][3] = 'NA';
    }
  }
}
# Processing first segment for each chr
foreach my $i (@idxstart) {
  if ($block[$i][2] - $block[$i][1] + 1 < $para_l &&
      ($block[$i][0] ne $block[$i + 1][0] ||
       $block[$i][3] ne $block[$i + 1][3]
      )
     ) {
    $block[$i][3] = ($block[$i][3] eq $block[$i + 1][3]) ? $block[$i][3] : 'NA';
  }
}

# Merge
my $nbpre = 0;
my $nb = scalar(@block);
while ($nb != $nbpre) {
  $nbpre = $nb;
  @block = merge(@block);
  $nb = scalar(@block);
}

# # Merge loh
# $nbpre = 0;
# $nb = scalar(@block);
# while ($nb != $nbpre) {
#   $nbpre = $nb;
#   @block = merge_loh_lohamp(@block);
#   $nb = scalar(@block);
# }

# Output
foreach my $r (@block) {
  my @b = @{$r};
  if ($b[3] ne 'NA') {
    # Calc log ratio
    my @brat;
    foreach my $i (@{$b[7]}) {
      push @brat, sprintf("%.2f", $rat[$i][10]);
    }
    if (@brat > 0) {
      my $medrat = median(@brat) - $lgrscale;
      print join("\t", $b[0], $b[1], $b[2], $statname{$b[3]}, $b[4], sprintf("%.2f",$medrat)), "\n";
      # print join("\t", $b[0], $b[1], $b[2], $statname{$b[3]}, $b[4], sprintf("%.2f",$medrat), join(',', @brat)), "\n";
    }
  }
}

sub meansd {
  my $sum = 0;
  my $n = 0;
  my $mean = 0;
  my $sd = 0;
  foreach my $e (@_) {
    $sum += $e;
    $n++;
  }
  if ($n == 1) {
    $mean = $sum;
  }
  elsif ($n > 1) {
    $mean = $sum / $n;
    $sum = 0;
    foreach my $e (@_) {
      $sum += ($e - $mean) ** 2;
    }
    $sd = ($sum / ($n - 1)) ** (1/2);
  }
  return ($mean, $sd);
}

sub median {
  my @e = sort {$a <=> $b} @_;
  return (@e % 2 == 1) ? $e[$#e / 2] : ($e[($#e - 1) / 2] + $e[($#e + 1) / 2]) / 2;
}

sub quantile {
  my ($ref, $p) = @_;
  my @e = sort {$a <=> $b} @{$ref};
  return $e[int($p * @e) - 1];
}

# rank($v, \@array)
sub rank {
  my $v = $_[0];
  my @e = sort {$a <=> $b} @{$_[1]};
  foreach my $i (0 .. $#e) {
    if ($v == $e[$i]) {
      return $i;
    }
  }
  return -1;
}

sub merge {
  my @block = @_;
  my @blockfinal;
  my $k = -1;
  foreach my $i (0 .. $#block) {
    next if ($block[$i][3] eq 'NA');
    if ($k == -1 || 
	$blockfinal[$k][0] ne $block[$i][0] || 
	$blockfinal[$k][2] + $para_gaplen <= $block[$i][1] ||
	$blockfinal[$k][3] ne $block[$i][3]
       ) {
      $k++;
      $blockfinal[$k] = [ @{$block[$i]} ];
    }
    else {
      $blockfinal[$k][2] = $block[$i][2];
      push @{$blockfinal[$k][7]}, @{$block[$i][7]};
    }
  }
  return @blockfinal;
}

sub merge_loh_lohamp {
  my @block = @_;
  my @blockfinal;
  my $k = -1;
  foreach my $i (0 .. $#block) {
    next if ($block[$i][3] eq 'NA');
    # Recalculate lgr median
    # my $blockfinallgr = ratlgrmedian(@rat[@{$blockfinal[$k][7]}]);
    # my $blocklgr = ratlgrmedian(@rat[@{$block[$k][7]}]);
    my ($blockfinallgr, $blocklgr) = (0, 0);
    if ($k != -1) {
      $blockfinallgr = ratlgrmedian(@rat[@{$blockfinal[$k][7]}]);
      $blocklgr = ratlgrmedian(@rat[@{$block[$k][7]}]);
    }
    if ($k != -1 &&
        $blockfinal[$k][0] eq $block[$i][0]
        && $blockfinal[$k][2] + $para_gaplen > $block[$i][1]
        && ($blockfinal[$k][3] eq 'loh' || $blockfinal[$k][3] eq 'loh.amp' || $blockfinal[$k][3] eq 'loss')
        && ($block[$i][3] eq 'loh' || $block[$i][3] eq 'loh.amp' || $block[$i][3] eq 'loss')
        && $blockfinallgr - $blocklgr < $para_diff_loh_lohamp
        && $blocklgr - $blockfinallgr < $para_diff_loh_lohamp
        # && $blockfinal[$k][5] - $block[$i][5] < $para_diff_loh_lohamp
        # && $block[$i][5] - $blockfinal[$k][5] < $para_diff_loh_lohamp
       ) {
      $blockfinal[$k][2] = $block[$i][2];
      if ($blockfinal[$k][3] ne $block[$i][3]) {
        $blockfinal[$k][3] = 'loh';
      }
      push @{$blockfinal[$k][7]}, @{$block[$i][7]};
    }
    else {
      $k++;
      $blockfinal[$k] = [ @{$block[$i]} ];
    }
  }
  return @blockfinal;
}

sub maf {
  my ($cn, $ref) = @_;
  my @frq = @{$ref};
  if (@frq != 4) {
    return '.';
  }
  if ($cn == 2) {
    if (($frq[1] + $frq[2]) / 2 > $para_meanlow && ($frq[1] + $frq[2]) / 2 < $para_meanhigh) {
      my $maf = $frq[1] + 0.5 - ($frq[1] + $frq[2]) / 2;
      $maf = ($maf < 0.5) ? $maf : (1 - $maf);
      return ($maf < $para_maftrim) ? $maf : 0.5;
    }
    elsif (abs($frq[1] - 0.5) < abs($frq[2] - 0.5)) {
      return ($frq[1] < 0.5) ? $frq[1] : (1 - $frq[1]);
    }
    else {
      return ($frq[2] < 0.5) ? $frq[2] : (1 - $frq[2]);
    }
  }
  if ($cn == 1 || $cn == 3) {
    if (($frq[1] + $frq[2]) / 2 > $para_meanlow && ($frq[1] + $frq[2]) / 2 < $para_meanhigh) {
      my $maf = $frq[1] + 0.5 - ($frq[1] + $frq[2]) / 2;
      $maf = ($maf < 0.5) ? $maf : (1 - $maf);
      return ($maf < $para_maftrim) ? $maf : 0.5;
    }
    else {
      return 0;
    }
  }
  return '.';
}

sub callcn {
  my $ref = shift;
  my @rat = @{$ref};
  foreach my $i (0 .. $#rat) {
    my $cstat = '';
    my ($cchr, $cstart, $cend, $clgru, $gc, $gcg, $clu, $ccpn, $cfrq, $gcgc, $clgr) = @{$rat[$i]};
    if ($cfrq ne '.') {
      $cfrq = median(split /;/, $cfrq);
      if ($cfrq < $para_cpn2diffloose) {
        if ($clgr > $low2 && $clgr < $high2) {
          $cstat = 'norm';
          $ccpn = 2;
        } elsif ($clgr > $low4) {
          $cstat = 'qua';
          $ccpn = 4;
        } elsif ($clgr < $low2) {
          $cstat = 'loss';
          $ccpn = 0;
        }
      }
      elsif ($cfrq > $para_cpn2diffstrict) {
        if ($clgr < $low3) {
          if ($clgr < $highloh) {
            $cstat = 'loh';
            $ccpn = 1;
          } else {
            $cstat = 'loh.amp';
            $ccpn = 2;
          }
        } elsif ($clgr > $low3) {
          $cstat = 'tri';
          $ccpn = 3;
        }
      }
      else {
        $cstat = 'norm';
        $ccpn = 2;
        my $lgr = $clgr - $lgrscale;
        if ($lgr < 0) {
          my $plgr = 2 - 2 * (2 ** $lgr);
          my $pmaf = 2 - 1 / (0.5 + $cfrq / 2);
          if (abs($plgr - $pmaf) < $para_p1diff) {
            $cstat = 'loh';
            $ccpn = 1;
          }
        }
        elsif ($lgr > 0) {
          my $plgr = 2 * (2 ** $lgr - 1);
          my $pmaf = 1 / (0.5 - $cfrq / 2) - 2;
          if (abs($plgr - $pmaf) < $para_p3diff) {
            $cstat = 'tri';
            $ccpn = 3;
          }
        }
      }
    }
    @{$rat[$i]} = ($cchr, $cstart, $cend, $clgru, $gc, $gcg, $clu, $ccpn, $cfrq, $gcgc, $clgr, $cstat);
  }
}

sub calcstat {
  my ($lgrscale, $ref) = @_;
  my @rat = @{$ref};
  my $sum = 0;
  my $lenlof = 0;
  my $lentri = 0;
  foreach my $i (0 .. $#rat) {
    my ($cchr, $cstart, $cend, $clgru, $gc, $gcg, $clu, $ccpn, $cfrq, $gcgc, $clgr, $cstat) = @{$rat[$i]};
    my @frq = split /,/, $cfrq;
    print STDERR $cstat, "\t", $cfrq, "\n";
    next if (@frq != 4);
    if ($cstat eq 'loh') {
      my $pdep = 2 - 2 * (2 ** ($clgr - $lgrscale));
      my $pmaf = 2 - 1 / (1 - maf(@frq));
      $sum += ($pdep - $pmaf) ** 2;
      $lenlof += $cend - $cstart + 1;
    }
    elsif ($cstat eq 'tri') {
      my $pdep = 2 * (2 ** ($clgr - $lgrscale) - 1);
      my $pmaf = 1 / maf(@frq) - 2;
      $sum += ($pdep - $pmaf) ** 2;
      $lentri += $cend - $cstart + 1;
    }
  }
  return ($sum, $lenlof, $lentri);
}

sub gcscale {
  my ($gcrat, $gcscale) = @_;
  # my %gcrat = %{$refgcrat};
  # my %gcscale = %{$refgcscale};
  foreach my $gc (sort keys %{$gcrat}) {
    $gcscale -> {$gc} = median(@{$gcrat -> {$gc}});
  }
  my @gcs = sort {$a <=> $b} keys %{$gcscale};
  my $gc = 0;
  while ($gc <= 1) {
    if (! exists $gcscale -> {$gc}) {
      if ($gc < $gcs[0]) {
      	$gcscale -> {$gc} = ($gcscale -> {$gcs[0]} - $gcscale -> {$gcs[1]}) / ($gcs[0] - $gcs[1]) * ($gc - $gcs[0]) + $gcscale -> {$gcs[0]};
      }
      elsif ($gc > $gcs[$#gcs]) {
      	$gcscale -> {$gc} = ($gcscale -> {$gcs[$#gcs - 1]} - $gcscale -> {$gcs[$#gcs]}) / ($gcs[$#gcs - 1] - $gcs[$#gcs]) * ($gc - $gcs[$#gcs]) + $gcscale -> {$gcs[$#gcs]};
      }
      else {
      	my $i = 0;
      	while ($gc > $gcs[++$i]){};
      	$gcscale -> {$gc} = ($gcscale -> {$gcs[$i - 1]} - $gcscale -> {$gcs[$i]}) / ($gcs[$i - 1] - $gcs[$i]) * ($gc - $gcs[$i]) + $gcscale -> {$gcs[$i]};
      }
    }
    $gc += 0.001;
    $gc = sprintf("%.3f", $gc);
  }
}

sub define_cutoff {
  $low1 = $lgr1 - 3 * $lgrsd + $lgrscale;
  $high1 = $lgr2 - 2 * $lgrsd + $lgrscale;
  $low2 = $lgr2 - 3 * $lgrsd + $lgrscale;
  $high2 = $lgr2 + 3 * $lgrsd + $lgrscale;
  $low3 = $lgr2 + $lgrsd + $lgrscale;
  $high3 = $lgr3 + 3 * $lgrsd + $lgrscale;
  $low4 = $lgr2 + 4 * $lgrsd + $lgrscale;
  $high4 = $lgr4 + 3 * $lgrsd + $lgrscale;
  $lowloh = $lgr1 - 3 * $lgrsd + $lgrscale;
  $highloh = $lgr2 - 3 * $lgrsd + $lgrscale;
  $lowlohamp = $lgr2 - 2 * $lgrsd + $lgrscale;
  $highlohamp = $lgr2 + 2 * $lgrsd + $lgrscale;
  $lowlohtri = $lgr2 + 3 * $lgrsd + $lgrscale;
  $highlohtri = $lgr3 + 3 * $lgrsd + $lgrscale;
  $highlohnorm = $lgr2 + $lgrscale;
}

sub asigncpn2 {
  my ($e, $cpn2rat, $cpn2len, $cpn2end) = @_;
  push @{$cpn2rat -> {$e -> [0]}}, $e -> [3];
  if (exists $cpn2end -> {$e -> [0]}) {
    $cpn2len -> {$e -> [0]} += ($e -> [1] <= $cpn2end -> {$e -> [0]}) ? ($e -> [2] - $cpn2end -> {$e -> [0]}) : ($e -> [2] - $e -> [1] + 1);
  }
  else {
    $cpn2len -> {$e -> [0]} = $e -> [2] - $e -> [1] + 1;
  }
  $cpn2end -> {$e -> [0]} = $e -> [2];
}

sub applygcscale {
  my ($rat, $gcscale) = @_;
  if ($filea ne '') {
    open FILEA, ">$filea"
      || die "Cannot open $filea: $!\n";
  }
  foreach my $r (@{$rat}) {
    if (exists $gcscale -> {$r -> [4]}) {
      $r -> [9] = $gcscale -> {$r -> [4]};
      $r -> [10] = $r -> [3] - $gcscale -> {$r -> [4]};
    }
    else {
      $r -> [9] = 0;
      $r -> [10] = $r -> [3] - 0;
      print STDERR join("\t", 'No GC category found: ', @{$r}), "\n";
    }
    if ($filea ne '') {
      my @frq = split /,/, $r -> [8];
      print FILEA join("\t", @{$r}, maf($r -> [7], \@frq)), "\n";
    }
  }
  if ($filea ne '') {
    close FILEA;
  }
}

sub cpn2stat {
  my @rat = @{$_[0]};
  my %ratchridxstart = %{$_[1]};
  my %ratchridxend = %{$_[2]};
  my $lgrscale = 0;
  my $lgrsd = 0;
  my @cpn2stat;
  my $cpn2sdmin = 100;
  foreach my $chr (keys %ratchridxstart) {
    my $len = 0;
    my $end = 0;
    my @lgrs;
    foreach my $r (@rat[$ratchridxstart{$chr} .. $ratchridxend{$chr}]) {
      my @e = @{$r};
      if ($e[7] eq '2') {
        push @lgrs, $e[10];
        if ($e[1] > $end) {
          $len += $e[2] - $e[1] + 1;
          $end = $e[2];
        }
        else {
          $len += $e[2] - $end + 1;
          $end = $e[2];
        }
      }
    }
    next if (@lgrs == 0);
    my ($mean, $sd) = meansd(@lgrs);
    my $med = median(@lgrs);
    push @cpn2stat, [ $chr, $len, $med, $mean, $sd ];
    if ($cpn2sdmin > $sd) {
      $cpn2sdmin = $sd;
    }
    print STDERR join("\t", $chr, $len, sprintf("%.3f", $med), sprintf("%.3f", $mean), $sd), "\n";
  }
  print STDERR '--------------', "\n";
  if (@cpn2stat > 0) {
    my @cpn2rat;
    my @cpn2sd;
    foreach my $ref (@cpn2stat) {
      my ($chr, $len, $med, $mean, $sd) = @{$ref};
      if ($sd <= 2 * $cpn2sdmin && $len >= $para_cpn2lenmin) {
        push @cpn2rat, $med;
        push @cpn2sd, $sd;
        print STDERR join("\t", $chr, $len, sprintf("%.3f", $med), sprintf("%.3f", $mean), $sd), "\n";
      }
    }
    if (@cpn2rat > 0) {
      $lgrscale = median(@cpn2rat);
      $lgrsd = median(@cpn2sd);
      if ($lgrsd < $para_sdmin) {
        $lgrsd = $para_sdmin;
      }
    }
  }
  print STDERR '-----------', "\n", 'lgrscale: ', $lgrscale, "\n", 'sd: ', $lgrsd, "\n", '------------------', "\n";
  if ($lgrsd < $para_sdmin) {
    $lgrsd = $para_sdmin;
  }
  return ($lgrscale, $lgrsd)
}

sub ratlgrmedian {
  my @rat = @_;
  my @lgrs;
  foreach my $r (@rat) {
    next if ($r -> [10] eq '.' || $r -> [10] eq 'NA');
    push @lgrs, $r -> [10];
  }
  return (@lgrs > 0) ? median(@lgrs): 'NA';
}

