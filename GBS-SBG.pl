#!/usr/bin/perl -w
#
# Script for assembly-based calling of Group B Streptococcus serotypes for WGS data
#
# Source data should always be up-to-date at https://github.com/swainechen/GBS-SBG
#
use warnings;
use strict;
use Cwd;
use File::Copy;
use File::Spec;
use File::Fetch;
use File::Temp;
use File::Basename;
use IPC::Open3;
use Getopt::Long;
Getopt::Long::Configure("pass_through");

my $BLASTN_BIN = "";
my $MIN_BLASTN_VERSION = "2.7.0";
my $BLASTN_COMMAND;
my $GITHUB_REPO = "https://github.com/swainechen/GBS-SBG";
my $GITHUB_REF_FASTA = "https://raw.githubusercontent.com/swainechen/GBS-SBG/main/GBS-SBG.fasta";
my $REF_FASTA_BASE = "GBS-SBG.fasta";
my $REF_FASTA_FULL = "";
my $USE_STDIN = 0;
my $ORIGDIR = File::Spec->rel2abs(getcwd);
my $TEMPDIR;
my $INFILE;
my $BLASTOUT;
my $NAME = "";
# this hash ref holds the cutoff values for reporting uncertainty
# these are different from the cutoffs used for making a serotype call, which
# are generally 90% identity over 90% coverage
my $UNCERTAINTY = { 'MAX_NEXT' => 0.95,
                    'MIN_COV' => 0.99,
                    'MIN_PID' => 99,	# this is 0-100 as per BLASTN output
                    'MAX_HITS' => 1,
                    'MAX_CONTIGS' => 1,
                    'MIN_NEXT_PRINT' => 0.9	# this is for when to not even
                                                # bother printing the next best
                  };
my ($d, $i, $t);
my (@f, @g);
my @index;
my $ff;
my (@start, @end, @id);
my (@start_c, @end_c, @id_c);
my $contig_count;
my $contig_sum;
my $serotype;

# keyed on serotype, then should have SCORE, LENGTH, COVERAGE, FOOTPRINT, ID, HITS, NEXTBEST, COVERAGES_CONTIG
my $sdata;

my @uncertainty;
my $contig;
my ($pid, $chld_in, $chld_out, $chld_err);
my ($max_cov, $max_id);
my $pass;
my $show_help = 0;
my $DEBUG = 0;
my $VERBOSE = 0;	# show all hits
my $BEST_ONLY = 0;	# better for automation if this is set

# Some internal parameters
my $BLAST_MIN_PID = 90;	# minimum percent ID to consider an individual hit
my $BLAST_MIN_LEN = 100;	# minimum alignment length to consider an individual hit
my $TOTAL_MIN_COVERAGE = 0.9;	# minimum fraction of the serotype to give a positive call

GetOptions (
  'help!' => \$show_help,
  'name=s' => \$NAME,
  'blast=s' => \$BLASTN_BIN,
  'ref=s' => \$REF_FASTA_FULL,
  'best!' => \$BEST_ONLY,
  'verbose!' => \$VERBOSE,
  'debug!' => \$DEBUG
);

# figure out the input filename - mostly if it's on STDIN
if (defined $ARGV[0] && $ARGV[0] ne "" && $ARGV[0] ne "-") {
  if (!-f $ARGV[0]) {
    $show_help = 1;
    print STDERR "Can't find input file $ARGV[0]\n\n";
  } else {
    $INFILE = File::Spec->rel2abs($ARGV[0]);
    if ($NAME eq "") {
      $NAME = basename($INFILE);
    }
  }
} else {
  # we'll try STDIN - set this here so we can check everything else first
  $USE_STDIN = 1;
  if ($NAME eq "") {
    $NAME = "STDIN";
  }
}

if ($show_help) {
  &print_help;
  exit;
}
$VERBOSE = 1 if $DEBUG;

# make sure we have everything we need
if ($DEBUG) {
  $TEMPDIR = File::Temp::tempdir( CLEANUP => 0 );
  print STDERR "[GBS-SBG INFO] Using temp directory $TEMPDIR\n";
} else {
  $TEMPDIR = File::Temp::tempdir( CLEANUP => 1 );
}
$BLASTN_BIN = get_blastn($BLASTN_BIN);
$DEBUG && print STDERR "[GBS-SBG INFO] Found blastn $BLASTN_BIN\n";
if (-f $REF_FASTA_BASE && $REF_FASTA_FULL eq "") {
  $REF_FASTA_FULL = File::Spec->rel2abs($REF_FASTA_BASE);
} else {
  if ($REF_FASTA_FULL eq "") {
    # try a few guesses
    if (eval "use FindBin; 1") {
      use FindBin;
      if (-f File::Spec->catfile($FindBin::Bin, $REF_FASTA_BASE)) {
        $REF_FASTA_FULL = File::Spec->catfile($FindBin::Bin, $REF_FASTA_BASE);
      } elsif (-f File::Spec->catfile($FindBin::Bin, "..", "lib", $REF_FASTA_BASE)) {
        $REF_FASTA_FULL = File::Spec->catfile($FindBin::Bin, "..", "lib", $REF_FASTA_BASE);
      }
    }
  } elsif (-d $REF_FASTA_FULL) {
    if (-f File::Spec->catfile($REF_FASTA_FULL, $REF_FASTA_BASE)) {
      $REF_FASTA_FULL = File::Spec->catfile($FindBin::Bin, $REF_FASTA_BASE);
    }
  } elsif (-f $REF_FASTA_FULL) {
    $REF_FASTA_FULL = File::Spec->rel2abs($REF_FASTA_FULL);
  }
}
if (-f $REF_FASTA_FULL) {
  if (!-f "$REF_FASTA_FULL.nin") {
    copy($REF_FASTA_FULL, File::Spec->catfile($TEMPDIR, $REF_FASTA_BASE));
    $REF_FASTA_FULL = File::Spec->catfile($TEMPDIR, $REF_FASTA_BASE);
    chdir($TEMPDIR);
    $pid = open3($chld_in, $chld_out, $chld_err, "makeblastdb -in $REF_FASTA_FULL -dbtype nucl");
    waitpid($pid, 0);
    die("Some problem with running makeblastdb...exit $?\n") if $?;
  }
}
if ($REF_FASTA_FULL eq "" || !-f $REF_FASTA_FULL) {
  $ff = File::Fetch->new(uri => "$GITHUB_REF_FASTA");
  $REF_FASTA_FULL = $ff->fetch( to => $TEMPDIR );
  chdir($TEMPDIR);
  $pid = open3($chld_in, $chld_out, $chld_err, "makeblastdb -in $REF_FASTA_FULL -dbtype nucl");
  waitpid($pid, 0);
  die("Some problem with running makeblastdb...exit $?\n") if $?;
}
if (!-f $REF_FASTA_FULL) {
  die("Couldn't find $REF_FASTA_BASE and couldn't download it; set the location with the -ref command line parameter\n");
}
$DEBUG && print STDERR "[GBS-SBG INFO] Using reference file $REF_FASTA_FULL\n";

# run the blast
chdir($TEMPDIR);
if ($USE_STDIN) {
  @f = <>;
  open F, ">" . File::Spec->rel2abs(File::Spec->catfile($TEMPDIR, "STDIN_INPUT"));
  print F @f;
  close F;
  $INFILE = File::Spec->rel2abs(File::Spec->catfile($TEMPDIR, "STDIN_INPUT"));
}
$BLASTOUT = File::Spec->rel2abs(File::Spec->catfile($TEMPDIR, "BLAST_OUTPUT"));
$BLASTN_COMMAND = "blastn -db $REF_FASTA_FULL -query $INFILE -out $BLASTOUT -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send sstrand evalue bitscore'";
$DEBUG && print STDERR "[GBS-SBG INFO] Running blast: $BLASTN_COMMAND\n";
$pid = open3($chld_in, $chld_out, $chld_err, $BLASTN_COMMAND);
waitpid($pid, 0);
die ("Error $? with blast command...\n") if $?;

# parse the results
open F, $BLASTOUT;
@f = <F>;
close F;
$d = ();
$DEBUG && print "== BLASTN Results ==\n", @f, "== END BLASTN Results ==\n";

# only take hits that are >90% identical
# data structure: $d->{Serotype}->{Query}->{i}->{SSTART}
# data structure: $d->{Serotype}->{Query}->{i}->{SEND}
# data structure: $d->{Serotype}->{Query}->{i}->{PID}
# data structure: $d->{Serotype}->{Query}->{i}->{LEN}
# data structure: $d->{Serotype}->{Query}->{i}->{BITSCORE}
# i is just a counter, i.e. line number in the blast output, can treat like a unique ID
# columns:
# 0 - qseqid
# 1 - sseqid
# 2 - pident (0-100)
# 3 - length
# 4 - mismatch
# 5 - gapopen
# 6 - qlen
# 7 - qstart
# 8 - qend
# 9 - slen
# 10 - sstart
# 11 - send
# 12 - sstrand
# 13 - evalue
# 14 - bitscore
foreach $i (0..$#f) {
  chomp $f[$i];
  @g = split /\t/, $f[$i];
  next if ($g[2] <= $BLAST_MIN_PID);
  next if ($g[3] <= $BLAST_MIN_LEN);
  $serotype = parse_serotype($g[1]);
  if (!defined $sdata->{$serotype}->{SCORE}) {
    $sdata->{$serotype}->{SCORE} = 0;
    $sdata->{$serotype}->{LENGTH} = $g[9];
  }
  $sdata->{$serotype}->{SCORE} += $g[14];	# bitscore
  $d->{$serotype}->{$g[0]}->{$i}->{PID} = $g[2];
  $d->{$serotype}->{$g[0]}->{$i}->{LEN} = $g[3];
  if ($g[12] eq "minus" || $g[12] eq "-" || $g[11] < $g[10]) {
    $d->{$serotype}->{$g[0]}->{$i}->{SSTART} = $g[11];
    $d->{$serotype}->{$g[0]}->{$i}->{SEND} = $g[10];
  } else {
    $d->{$serotype}->{$g[0]}->{$i}->{SSTART} = $g[10];
    $d->{$serotype}->{$g[0]}->{$i}->{SEND} = $g[11];
  }
  $d->{$serotype}->{$g[0]}->{$i}->{BITSCORE} = $g[14];
}
# main call is based on highest total bitscore, but have to check for >90% coverage
# collapse to highest ID per serotype to do some checking
foreach $serotype (keys %$sdata) {
  $sdata->{$serotype}->{COVERAGE} = 0;
  $sdata->{$serotype}->{CONTIGCOV} = ();
  @start = ();
  @end = ();
  @id = ();
  foreach $contig (keys %{$d->{$serotype}}) {
    $t = $d->{$serotype}->{$contig};
    # we have to get per-contig data in this loop, then accumulate it for overall
    @start_c = ();
    @end_c = ();
    @id_c = ();
    foreach $i (sort {$t->{$a}->{SSTART} <=> $t->{$b}->{SSTART}} keys %{$t}) {
      push @start_c, $t->{$i}->{SSTART};
      push @end_c, $t->{$i}->{SEND};
      push @id_c, $t->{$i}->{PID};
      push @start, $t->{$i}->{SSTART};
      push @end, $t->{$i}->{SEND};
      push @id, $t->{$i}->{PID};
    }
    # merge from the back
    foreach $i (reverse 1..$#start_c) {
      # @start_c should already be sorted by start coordinates
      if ($start_c[$i] <= $end_c[$i-1]) {
        # aggregate the identity here
        if ($id_c[$i-1] > $id_c[$i]) {
          $id_c[$i-1] = ( ($end_c[$i-1] - $start_c[$i-1] + 1) * $id_c[$i-1]  + 
                          ($end_c[$i] - $end_c[$i-1]) * $id_c[$i]          ) /
                        ( $end_c[$i] - $start_c[$i-1] + 1 );
        } else {
          $id_c[$i-1] = ( ($start_c[$i] - $start_c[$i-1]) * $id_c[$i-1]  +
                          ($end_c[$i] - $start_c[$i] + 1) * $id_c[$i]  ) /
                        ( $end_c[$i] - $start_c[$i-1] + 1 );
        }
        $end_c[$i-1] = $end_c[$i];
        splice @start_c, $i, 1;
        splice @end_c, $i, 1;
        splice @id_c, $i, 1;
      }
    }
    # we should have nonoverlapping intervals now
    foreach $i (0..$#start_c) {
      $sdata->{$serotype}->{CONTIGCOV}->{$contig} += ($end_c[$i] - $start_c[$i] + 1) / $sdata->{$serotype}->{LENGTH};
    }
  }
  # need to re-sort the coordinates because we aggregated across contigs
  @index = sort {$start[$a] <=> $start[$b]} 0..$#start;
  @start = @start[@index];
  @end = @end[@index];
  @id = @id[@index];
  $sdata->{$serotype}->{HITS} = scalar @start;
  # merge from the back
  foreach $i (reverse 1..$#start) {
    # @start should already be sorted by start coordinates
    if ($start[$i] <= $end[$i-1]) {
      if ($end[$i-1] < $end[$i]) {
        # aggregate the identity here
        if ($id[$i-1] > $id[$i]) {
          $id[$i-1] = ( ($end[$i-1] - $start[$i-1] + 1) * $id[$i-1]  +
                        ($end[$i] - $end[$i-1]) * $id[$i]          ) /
                      ( $end[$i] - $start[$i-1] + 1 );
        } else {
          $id[$i-1] = ( ($start[$i] - $start[$i-1]) * $id[$i-1]  +
                        ($end[$i] - $start[$i] + 1) * $id[$i]  ) /
                      ( $end[$i] - $start[$i-1] + 1 );
        }
        $end[$i-1] = $end[$i];
        splice @start, $i, 1;
        splice @end, $i, 1;
        splice @id, $i, 1;
      } else {
      # interval $i is completely contained within interval $i-1
      # if $id[$i] > $id[$i-1] we could increase the overall identity...but
      # this doesn't seem to make sense, since at this point interval $i-1 is
      # a single blast hit, so instead just get rid of interval $i
        splice @start, $i, 1;
        splice @end, $i, 1;
        splice @id, $i, 1;
      }
    }
  }
  # we should have nonoverlapping intervals now
  foreach $i (0..$#start) {
    $sdata->{$serotype}->{COVERAGE} += ($end[$i] - $start[$i] + 1);
    $sdata->{$serotype}->{ID} += ($end[$i] - $start[$i] + 1) * $id[$i];
  }
  $sdata->{$serotype}->{ID} /= $sdata->{$serotype}->{COVERAGE};
  $sdata->{$serotype}->{COVERAGE} /= $sdata->{$serotype}->{LENGTH};
}
$DEBUG && print "[GBS-SBG INFO] Main output follows below:\n";
if ($VERBOSE) {
  print join ("\t", "# Name", "Serotype", "Total BitScore", "Total Coverage", "Percent ID", "Number of BLASTN Hits", "Number Contigs", "NextBest BitScore"), "\n";
} else {
  print join ("\t", "# Name", "Serotype", "Uncertainty"), "\n";
}
# calculate the bitscore differential to the next best hit
$i = 0;
foreach $serotype (sort {$sdata->{$a}->{SCORE} <=> $sdata->{$b}->{SCORE}} keys %$sdata) {
  $sdata->{$serotype}->{NEXTBEST} = $i/$sdata->{$serotype}->{SCORE};
  $i = $sdata->{$serotype}->{SCORE};
}
# first ensure we don't have a nontypeable
$pass = 0;
$max_cov = 0;
$max_id = 0;
foreach $serotype (sort {$sdata->{$b}->{SCORE} <=> $sdata->{$a}->{SCORE}} keys %$sdata) {
  if ($sdata->{$serotype}->{COVERAGE} >= $TOTAL_MIN_COVERAGE &&
      $sdata->{$serotype}->{ID} >= $BLAST_MIN_PID) {
    $pass = 1;
  }
  $max_cov = $sdata->{$serotype}->{COVERAGE} if $max_cov < $sdata->{$serotype}->{COVERAGE};
  $max_id = $sdata->{$serotype}->{ID} if $max_id < $sdata->{$serotype}->{ID};
}
if (!$pass) {
  if ($VERBOSE) {
    print join ("\t", $NAME, "NT", "0", "0", "0", "0", "0", "0"), "\n";
  } else {
    print join ("\t", $NAME, "NT", join (";", "MaxCov:$max_cov", "MaxID:$max_id")), "\n";
  }
}

# we will print in decreasing order of total bitscore
# $i here is used to count output lines if not debug (max 2)
$i = 0;
foreach $serotype (sort {$sdata->{$b}->{SCORE} <=> $sdata->{$a}->{SCORE}} keys %$sdata) {
  # count # of contigs required
  $contig_count = 0;
  $contig_sum = 0;
  foreach $contig (sort {$sdata->{$serotype}->{CONTIGCOV}->{$b} <=> $sdata->{$serotype}->{CONTIGCOV}->{$a}} keys %{$sdata->{$serotype}->{CONTIGCOV}}) {
    $contig_sum += $sdata->{$serotype}->{CONTIGCOV}->{$contig};
    $contig_count++;
    last if $contig_sum >= $sdata->{$serotype}->{COVERAGE};
  }
  if ($VERBOSE) {
    print join ("\t", $NAME, $serotype, $sdata->{$serotype}->{SCORE}, $sdata->{$serotype}->{COVERAGE}, $sdata->{$serotype}->{ID}, $sdata->{$serotype}->{HITS}, $contig_count, $sdata->{$serotype}->{NEXTBEST}), "\n";
  } else {
    last if !$pass && $BEST_ONLY;
    $i = 1 if !$pass && $i == 0;
    print join ("\t", $NAME, $serotype);
    @uncertainty = ();
    if ($sdata->{$serotype}->{NEXTBEST} > $UNCERTAINTY->{MAX_NEXT}) {
      push @uncertainty, "NextBitScore:$sdata->{$serotype}->{NEXTBEST}";
    }
    if ($sdata->{$serotype}->{COVERAGE} < $UNCERTAINTY->{MIN_COV}) {
      push @uncertainty, "Coverage:$sdata->{$serotype}->{COVERAGE}";
    }
    if ($sdata->{$serotype}->{ID} < $UNCERTAINTY->{MIN_PID}) {
      push @uncertainty, "Pidentity:$sdata->{$serotype}->{ID}";
    }
    if ($sdata->{$serotype}->{HITS} > $UNCERTAINTY->{MAX_HITS}) {
      push @uncertainty, "BLASTNHits:$sdata->{$serotype}->{HITS}";
    }
    if ($contig_count > $UNCERTAINTY->{MAX_CONTIGS}) {
      push @uncertainty, "Contigs:$contig_count";
    }
    if (scalar @uncertainty) {
      print "\t", join (";", @uncertainty), "\n";
    } else {
      print "\n";
    }
    $i++;
    last if $BEST_ONLY;
    last if ($i >= 2);
    next if ($sdata->{$serotype}->{NEXTBEST} > $UNCERTAINTY->{MIN_NEXT_PRINT});
    last if ($sdata->{$serotype}->{COVERAGE} >= $TOTAL_MIN_COVERAGE && $sdata->{$serotype}->{ID} >= $BLAST_MIN_PID && $contig_count == 1);
  }
}

#
# subroutines below here
#
sub parse_serotype {
  # assuming SRST2 format here
  my ($s) = @_;
  my @f = split /__/, $s;
  if (defined $f[2]) {
    return $f[2];
  }
  return $s;
}

sub get_blastn {
  my ($init) = @_;
  my $full = "";
  my $output = "";
  my $pid;
  my ($chld_in, $chld_out, $chld_err);
  my @h;
  if ($init eq "") {
    $full = "blastn";
    $pid = open3($chld_in, $chld_out, $chld_err, "$full -h");
    $output = join ("", <$chld_out>);
    waitpid($pid, 0);
    if ($?) {
      die("Can't find blastn executable, can set this with the -blast command line parameter\n");
    } else {
      if (check_blastn_version($output)) {
        return($full);
      }
    }
  } else {
    $full = File::Spec->rel2abs($init);
    if (-d $full) {
      $full = File::Spec->catfile($full, "blastn");
    }
    if (-f $full) {
      $pid = open3($chld_in, $chld_out, $chld_err, "$full -h");
      $output = join ("", <$chld_out>);
      waitpid($pid, 0);
      if ($?) {
        die ("Can't find blastn executable at $init, can set this with the -blast command line parameter");
      }
    }
  }
  return($full);

  sub check_blastn_version {
    my ($output) = @_;
    my @h = split /\n/, $output;
    my $i;
    my $version_string;
    my ($major, $minor, $patch);
    my ($minmajor, $minminor, $minpatch) = split /\./, $MIN_BLASTN_VERSION;
    for ($i=0; $i < scalar @h; $i++) {
      last if $h[$i] =~ /^DESCRIPTION/;
    }
    if ($i < $#h) {
      if ($h[$i+1] =~ /BLAST (\d+)\.(\d+)\.(\d+)\+?/) {
        $version_string = $h[$i];
        $version_string =~ s/^.*BLAST //;
        $major = $1;
        $minor = $2;
        $patch = $3;
        if ($major > $minmajor) {
          return (1);
        } elsif ($major == $minmajor) {
          if ($minor > $minminor) {
            return (1);
          } elsif ($minor == $minminor) {
            if ($patch >= $minpatch) {
              return (1);
            }
          }
        }
        die ("Need minimum blastn version $MIN_BLASTN_VERSION - could only find version $version_string - can set a custom blastn path with the -blast command line parameter\n");
      }
    }
    die ("Could not find blastn or couldn't determine the version - can set a custom blastn path with the -blastn command line parameter\n");
  }
}

sub print_help {
  print <<__HELP__;
Usage: $0 <assembly_fasta_file> [ -name <string> ] [ -best ] [ -blastn <path_to_blastn> ] [ -ref <GBS-SBG references> ] [ -debug ] [ -verbose ]

<assembly_fasta_file> should be a regular multi-fasta file with assembled contigs or a complete genome.
You should specify the -name parameter, all output will be prefixed by that string. Defaults to the input filename.
The -best option only prints out one call (with possible uncertainty information). Default behavior is to also print the next best call if any of the following are true for the "best" serotype call:
 - Overall coverage < $TOTAL_MIN_COVERAGE
 - Overall percent ID < $BLAST_MIN_PID
 - Number of contigs > 1
 - Number of BLASTN hits > 1
 - BitScore for next best serotype is >$UNCERTAINTY->{MIN_NEXT_PRINT}*BitScore for the best serotype

Requires BLAST+ version $MIN_BLASTN_VERSION or above. Will look on your path, or you can specify the blastn binary with -blast.
Requires GBS-SBG references (typically $REF_FASTA_BASE), will try looking in a few places, otherwise will try to pull directly from $GITHUB_REPO.

The -debug switch prints debugging info to STDERR.
The -verbose switch prints more information and outputs data for all hits, even if they don't pass the cutoffs.

More complete documentation available at $GITHUB_REPO.
__HELP__
}
