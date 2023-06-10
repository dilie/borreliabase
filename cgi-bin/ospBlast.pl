#!/usr/bin/perl
#use warnings;
use strict;
use JSON;
use CGI;
use Bio::SearchIO;
use Bio::Seq;
use Bio::SeqUtils;

my $cgi = new CGI;
print $cgi->header(-type => "application/json", -charset => "utf-8");
my $seq = $cgi->param("seq");
my $ev = 0.0000000001;

my $seqobj = Bio::Seq->new(-id=>'tmp', -seq=>$seq);
my $type = $seqobj->alphabet();

if ($type eq 'rna'){ print to_json({rna=>1}); exit }
$seq = &trans($seqobj) if $type eq 'dna';

use lib $ENV{'DOCUMENT_ROOT'}.'/cgi-bin';
use BBConfig;

my %config = BBConfig::get_config();
my $dir = $config{'tmp_path'};
&printErr("No tmp_path defined in config") unless $dir;

my $dirDown = $config{'download_dir'};
&printErr("No download_dir defined in config") unless $dirDown;

my $blast_path = $config{'blastp_path'};
&printErr("no blast_path defined in config") unless $blast_path;

open IN, '>', $dir.'/blast.in';
print IN '>', 'tmp', "\n", $seq, "\n";
close IN;

my $db = 'ospC.pep';

my $blastout = `$blast_path -query $dir/blast.in -db $dirDown/$db -evalue $ev`;
&printErr("blast error") unless $blastout;
open my $fh, '<', \$blastout;

my $searchio = Bio::SearchIO->new(-format=>'blast', -fh=>$fh);

my @hsps;
my $result = $searchio->next_result();
while (my $hit = $result->next_hit){
  my $id = $hit->name();
  $id =~ s/lcl\|//;

  while (my $hsp = $hit->next_hsp) {
    my %obj = (
        id => $id * 1,
        ev => $hsp->evalue(),
        Positives => $hsp->frac_conserved(),
        Identities => $hsp->frac_identical(),
        Gaps => $hsp->gaps(),
        );
    push @hsps, \%obj
  }
}

my %h;
if (scalar @hsps){
  @hsps = sort {$a->{ev} <=> $b->{ev}} @hsps;
  %h = %{shift @hsps};
} else {
  %h = (no=>1)
}

print to_json(\%h);

sub printErr {
  print to_json({error=>shift}); exit
}

sub trans {
  my $obj = shift;
  my $s = $obj->translate()->seq();
  return $s if $s !~ /\*.+/;
  my $util = Bio::SeqUtils->new();
  my @pros = $util->translate_6frames($obj);
  shift @pros;
  foreach (@pros){
    $s = $_->seq();
    return $s if $s !~ /\*.+/
  }
  print to_json({error=>'No valid reading frame'});
  exit
}
