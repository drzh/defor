#!/usr/bin/perl

# Usage: prog -i FILE_nclust(default: 'stdin') -d FLOAT_mean_diff -n INT_min_count -f FLOAT_freqlow -F FLOAT_freqhigh -m FLOAT_meanlow -M FLOAT_meanhigh -l INT_min_length -v INT_verbose | <STDOUT>

use strict;
use warnings;

use Getopt::Long;

my $filei = '/dev/stdin';
my $para_meandiff = 0.1;
my $para_meanpairdiff = 0.05;
my $para_freqlow = 0.1;
my $para_freqhigh = 0.9;
my $para_meanlow = 0.45;
my $para_meanhigh = 0.6;
my $para_meanamplow = 0.25;
my $para_meanamphigh = 0.75;
my $para_meanampdiff = 0.4;
my $para_fold = 3;
my $para_fold_loh = 5;
my $para_len = 10000000;
my $para_n = 10;
my $para_v = 0;
my $para_mean12diff = 0.05;
my $para_cncutoff = 0.15;
my $para_mean03diff = 0.8;
my $para_fold_adj = 5;
my $para_fold_loh_adj = 1;

Getopt::Long::Configure("bundling");
GetOptions(
           'i:s' => \$filei,
	   'd:f' => \$para_meandiff,
	   'f:f' => \$para_freqlow,
	   'F:f' => \$para_freqhigh,
	   'm:f' => \$para_meanlow,
	   'M:f' => \$para_meanhigh,
	   'n:i' => \$para_n,
	   'l:i' => \$para_len,
	   'v:i' => \$para_v
          );

my $flag = 1;
my $sstart = 0;
my $send = 0;
my $estart = 0;
my $eend = 0;
my $chr = '';
my $flag_inter = 0;
my $cnpre = -1;
my $cn = -1;
my $nmidpre = -1;
my @meanpre = (0, 0, 0, 0);
my @res;
my $flag_pair = -1;
my @mean0;
my @mean1;
my @mean2;
my @mean3;
my $mean12diff;
my $mean12diffpre;
my @mean12diffs;

open FILE, "<$filei"
  || die "Can not open $filei: $!\n";
while (<FILE>) {
  chomp;
  my @e = split /\t/;
  if ($flag_pair == -1) {
    my $n = scalar(@e);
    if ($n == 6) {
      $flag_pair = 0;
      $para_meanpairdiff = 0;
    }
    elsif ($n == 9) {
      $flag_pair = 1;
    }
    else {
      print STDERR "Wrong input file!\n";
    }
  }
  my @mean = split /,/, $e[3];
  my @sd = split /,/, $e[4];
  my @cnt = split /,/, $e[5];
  my @meannorm = @mean;
  my @sdnorm = @sd;
  my @cntnorm = @cnt;
  if ($flag_pair == 1) {
    @mean = split /,/, $e[6];
    @sd = split /,/, $e[7];
    @cnt = split /,/, $e[8];
  }
  my $nlow = 0;
  my $nhigh = 0;
  my $nmid = 0;
  my $meanmid = 0;
  my $flag_mean = 0;
  foreach my $i (0 .. 3) {
    if ($mean[$i] > 0) {
      if ($mean[$i] < $para_freqlow) {
	$nlow++;
      }
      elsif ($mean[$i] > $para_freqhigh) {
        $nhigh++;
      }
      else {
	$nmid++;
	$meanmid = $mean[$i];
	if ($mean[$i] == $meanpre[$i]) {
	  $flag_mean++;
	}
      }
    }
  }
  if ($nmid > 2) {
    $nmid = 2;
  }
  my $flag_next = 0;
  my $nmidnorm = 0;
  my $meanmidnorm = 0;
  if ($flag_pair == 1) {
    foreach my $i (0 .. 3) {
      if ($sdnorm[$i] > 0) {
        if ($meannorm[$i] >= $para_freqlow && $meannorm[$i] <= $para_freqhigh) {
          $nmidnorm++;
          $meanmidnorm = $meannorm[$i];
        }
      }
    }
    if ($nmidnorm == 0) {
      $flag_next = 1;
    }
    elsif ($nmidnorm == 2) {
      if ($cntnorm[1] < $para_fold * $cntnorm[2]) {
        $meanmidnorm = $meannorm[2];
        $nmidnorm = 1;
      }
      elsif ($cntnorm[2] < $para_fold * $cntnorm[1]) {
        $meanmidnorm = $meannorm[1];
        $nmidnorm = 1;
      }
      elsif (($meannorm[1] + $meannorm[2]) / 2 < $para_meanlow ||
          ($meannorm[1] + $meannorm[2]) / 2 > $para_meanhigh
          # $meannorm[2] - $meannorm[1] > $para_meandiff
         ) {
        $flag_next = 1;
      }
    }
  }

  # next if ($flag_next == 1);
  if ($nmidnorm == 1) {
    $cntnorm[1] = 1;
    $cntnorm[2] = 1;
  }
  
# Estimate number of clusters
  if ($nmid == 2) {
    $mean12diff = abs($mean[2] - $mean[1]);
    # if ($cnt[1] < $para_fold * $cnt[2] && $cnt[2] < $para_fold * $cnt[1]) {
    if (($cntnorm[1] < $cntnorm[2] &&
         ($cnt[1] + $para_fold_adj) / ($cnt[1] + $cnt[2] + $para_fold_adj) < $para_fold * ($cntnorm[1] + $para_fold_adj) / ($cntnorm[1] + $cntnorm[2] + $para_fold_adj) && 
         ($cntnorm[1] + $para_fold_adj) / ($cntnorm[1] + $cntnorm[2] + $para_fold_adj) < $para_fold * ($cnt[1] + $para_fold_adj) / ($cnt[1] + $cnt[2] + $para_fold_adj)
        )
        ||
        ($cntnorm[1] >= $cntnorm[2] &&
         ($cnt[2] + $para_fold_adj) / ($cnt[1] + $cnt[2] + $para_fold_adj) < $para_fold * ($cntnorm[2] + $para_fold_adj) / ($cntnorm[1] + $cntnorm[2] + $para_fold_adj) && 
         ($cntnorm[2] + $para_fold_adj) / ($cntnorm[1] + $cntnorm[2] + $para_fold_adj) < $para_fold * ($cnt[2] + $para_fold_adj) / ($cnt[1] + $cnt[2] + $para_fold_adj)
        )
       ) {
      if (($cnt[1] + $cnt[2] + $para_fold_loh_adj) / ($cnt[0] + $cnt[1] + $cnt[2] + $cnt[3] + $para_fold_loh_adj) * $para_fold_loh > ($cntnorm[1] + $cntnorm[2] + $para_fold_loh_adj) / ($cntnorm[0] + $cntnorm[1] + $cntnorm[2] + $cntnorm[3] + $para_fold_loh_adj)) {
        if (($mean[1] + $mean[2]) / 2 < $para_meanlow ||
            ($mean[1] + $mean[2]) / 2 > $para_meanhigh) {
          $nmid = -1;
        }
      }
      else {
        $nmid = 0;
      }
    }
    else {
      $nmid = 1;
      my $idx = $cnt[1] > $cnt[2] ? 1 : 2;
      $meanmid = $mean[$idx];
    }
  }
  if ($nmid == 1) {
    $mean12diff = $para_mean12diff;
    if ($meanmid < $para_meanlow || $meanmid > $para_meanhigh || ($cnt[1] + $cnt[2] + $para_fold_loh_adj) / ($cnt[0] + $cnt[1] + $cnt[2] + $cnt[3] + $para_fold_loh_adj) * $para_fold_loh < ($cntnorm[1] + $cntnorm[2] + $para_fold_loh_adj) / ($cntnorm[0] + $cntnorm[1] + $cntnorm[2] + $cntnorm[3] + $para_fold_loh_adj)) {
      $nmid = 0;
    }
  }
  if ($nmid == 0) {
    $mean12diff = $para_mean03diff;
    if ($nlow + $nhigh == 0) {
      $nmid = -1;
    }
  }

  if ($para_v == 1) {
    print STDERR $_, "\t", $nmid, "\n";
  }

  # # Test
  # print join("\t", @e, $nmidnorm, $nmid, $mean12diff), "\n";

  next if ($nmid > 2 || $nmid < 0);

  # Call cluser segment
  my $flag_rep = 1;
  while ($flag_rep == 1) {
    $flag_rep = 0;
    if ($e[0] eq $chr) {
      if ($e[1] < $eend &&
          abs($mean12diff - median(@mean12diffs)) < $para_meandiff
	 ) {
	$estart = $e[1];
	$eend = $e[2];
        push @mean12diffs, $mean12diff;
      }
      else {
	$flag = segment($chr, $sstart, $send, $estart, $eend, median(@mean12diffs));
	if ($flag == 0) {
	  if (@res) {
	    my $r = pop @res;
	    ($chr, $sstart, $send, $estart, $eend) = @{$r};
	    $flag_rep = 1;
	  }
	  else {
	    $flag = 1;
	  }
	}
      }
    }
    else {
      segment($chr, $sstart, $send, $estart, $eend, median(@mean12diffs));
      $flag = 1;
    }
  }

  if ($flag == 1) {
    $chr = $e[0];
    $sstart = $e[1];
    $send = $e[2];
    $estart = $sstart;
    $eend = $send;
    $flag = 0;
    @mean12diffs = ();
    push @mean12diffs, $mean12diff;
  }
  $nmidpre = $nmid;
  @meanpre = @mean;
}
close FILE;

segment($chr, $sstart, $send, $estart, $eend, median(@mean12diffs));

if (@res > 0) {
  # Drop short segment
  $res[0][7] = $res[0][1];
  my @resfinal;
  foreach my $i (0 .. ($#res - 1)) {
    if ($res[$i][0] eq $res[$i + 1][0] && $res[$i][4] > $res[$i + 1][1]) {
      $res[$i][8] = int(($res[$i][4] + $res[$i + 1][1]) / 2);
      $res[$i + 1][7] = $res[$i][8] + 1;
    }
    else {
      $res[$i][8] = $res[$i][4];
      $res[$i + 1][7] = $res[$i + 1][1];
    }
    if ($res[$i][8] - $res[$i][7] >= $para_len) {
      push @resfinal, [ @{$res[$i]} ];
    }
  }
  $res[$#res][8] = $res[$#res][4];
  if ($res[$#res][8] - $res[$#res][7] >= $para_len) {
    push @resfinal, [ @{$res[$#res]} ];
  }
  @res = @resfinal;
  
  # Address overlap and output
  $res[0][7] = $res[0][1];
  foreach my $i (0 .. ($#res - 1)) {
    if ($res[$i][0] eq $res[$i + 1][0] && $res[$i][4] > $res[$i + 1][1]) {
      $res[$i][8] = int(($res[$i][4] + $res[$i + 1][1]) / 2);
      $res[$i + 1][7] = $res[$i][8] + 1;
    }
    else {
      $res[$i][8] = $res[$i][4];
      $res[$i + 1][7] = $res[$i + 1][1];
    }
    output($res[$i]);
    # print join("\t", $res[$i][0], $res[$i][7], $res[$i][8], sprintf("%.3f", $res[$i][5]), $res[$i][6]), "\n";
  }
  $res[$#res][8] = $res[$#res][4];
  output($res[$#res]);
  # print join("\t", $res[$#res][0], $res[$#res][7], $res[$#res][8], sprintf("%.3f", $res[$#res][5]), $res[$#res][6]), "\n";
}

sub segment {
  my ($chr, $sstart, $send, $estart, $eend, $med) = @_;
  if ($chr ne '') {
    push @res, [ @_ ];
    return 1;
  }
  else {
    return 0;
  }
}

sub median {
  my @e;
  foreach my $i (0 .. $#_) {
    if (defined $_[$i]) {
      push @e, $_[$i];
    }
  }
  # my @e = sort {$a <=> $b} @_;
  if (@e) {
    @e = sort {$a <=> $b} @e;
    return (@e % 2 == 1) ? $e[$#e / 2] : ($e[($#e - 1) / 2] + $e[($#e + 1) / 2]) / 2;
  }
  else {
    return 0;
  }
}

sub output {
  my $ref = $_[0];
  my $cn = '.';
  if ($ref -> [5] <= $para_cncutoff) {
    $cn = 2;
  }
  elsif ($ref -> [5] > $para_cncutoff) {
    $cn = 1;
  }
  print join("\t", $ref -> [0], $ref -> [7], $ref -> [8], $cn, sprintf("%.3f", $ref -> [5])), "\n";
}
