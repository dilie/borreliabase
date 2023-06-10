#!/usr/bin/perl
#use warnings;
use strict;
use JSON;
use CGI;
use Bio::SearchIO;
use Bio::Seq;
use DBI;

my $cgi = new CGI;
print $cgi->header(-type => "application/json", -charset => "utf-8");
my $ev = $cgi->param("ev");
my $gid = $cgi->param("gid");
my $seq = $cgi->param("seq");

my $seq_obj = Bio::Seq->new(-id=>'tmp', -seq => $seq);
my $type = $seq_obj->alphabet();

if ($type eq 'rna'){ print to_json({rna=>1}); exit }
#use lib $ENV{'DOCUMENT_ROOT'}.'/cgi-bin';
use lib 'cgi-bin';
use BBConfig;
my %config = BBConfig::get_config();
my $dir = $config{'tmp_path'} or die "No tmp_path defined in config";
my $dirDown = $config{'download_dir'} or die "No download_dir defined in config";
my $blast_path = $config{$type eq 'protein'? 'blastp_path' : 'blastn_path'} or die "no blast_path defined in config";

open IN, '>', $dir.'/blast.in';
print IN '>', 'tmp', "\n", $seq, "\n";
close IN;

my $db = $gid.'.'. ($type eq 'protein'? 'pep' : ($gid eq '100'? 'genome' : 'nuc'));

my $blastout;
if (length($seq)<=27 && $type eq 'dna'){
  $blastout = `$blast_path -task 'blastn-short' -query $dir/blast.in -db $dirDown/$db -evalue $ev` 
} else {
  $blastout = `$blast_path -query $dir/blast.in -db $dirDown/$db -evalue $ev`
}
open my $fh, '<', \$blastout;

my $searchio = Bio::SearchIO->new(-format=>'blast', -fh=>$fh);

my @hsps;
my $result = $searchio->next_result();
while(my $hit = $result->next_hit) {
  my $locus = $hit->name();
  $locus =~ s/lcl\|//;
  my $lenS = $hit->length();
    while (my $hsp = $hit->next_hsp) {
      my $hit_start = $hsp->start('hit');
      my $hit_stop = $hsp->end('hit');
      my @orfs;
      @orfs = @{&locateHit($locus, $hit_start, $hit_stop)} if $gid eq '100' && $type eq 'dna';

      my %obj = (
          locus => $locus,
          lenS => $lenS,
          ev => $hsp->evalue(),
          Positives => $hsp->frac_conserved(),
          Identities => $hsp->frac_identical(),
#          Gaps => $hsp->gaps(),
          seqQ => [$hsp->start('query'), $hsp->end('query'), $hsp->query_string],
          seqS => [$hit_start, $hit_stop, $hsp->hit_string],
          seqH => $hsp->homology_string,
#          strandQ => $hsp->strand('query'),
         );
      $obj{rev} = 1 if $hsp->strand('hit')==-1;
      $obj{orfs} = \@orfs if scalar @orfs;
      push @hsps, \%obj
     }
}

print to_json(scalar @hsps? \@hsps : {no=>1});
exit;

sub locateHit {
  my ($acc, @pos) = @_;
  @pos = sort {$a<=>$b} @pos;
  $acc =~ s/_.+$//;

  my $db_name = $config{'db_name'} or die "no db_name defined in config gc\n";
  my $db_host = $config{'db_host'} or die "no db_host defined in config\n";
  my $db_user = $config{'db_user'} or die "no db_user defined in config\n";
  my $db_pass = $config{'db_pass'} or die "no db_pass defined in config\n";
  my $dbh = DBI->connect("dbi:Pg:dbname=$db_name;host=$db_host", $db_user, $db_pass) or die;
  my $sth1=$dbh->prepare("select locus, start+shift_start, stop+shift_stop, strand, b.rep_id from orf4 JOIN contig4 b USING (con_acc) where con_acc=? and start+shift_start<=? order by start desc limit 1");
  my $sth2=$dbh->prepare("select locus, start+shift_start, stop+shift_stop, strand, b.rep_id from orf4 JOIN contig4 b USING (con_acc) where con_acc=? and stop+shift_stop>=? order by start limit 1");
  my $sth3=$dbh->prepare("select locus, start+shift_start, stop+shift_stop, strand from orf4 where con_acc=? and start+shift_start>=? and start+shift_start<=? order by start");

  $sth1->execute($acc, $pos[0]);
  my @data1 = $sth1->fetchrow();
  $sth1->finish();

  $sth2->execute($acc, $pos[1]);
  my @data2 = $sth2->fetchrow();
  $sth2->finish();

  my $cont;
  if (!$data1[0]){
    $cont = pop @data2;
    return [$cont, 0, \@data2]
  }
  $cont = pop @data1;
  if (!$data2[0]){
    return [$cont, \@data1, 0]
  }
  pop @data2;
  if ($data1[0] eq $data2[0]){
    if ($data1[-1]){ $data1[2] -= 3 } else { $data1[1] += 3 }
    return [$cont, \@data1]
  }

  my @res = ($cont);
  $sth3->execute($acc, $data1[1], $data2[1]);
  while (my @data = $sth3->fetchrow()) {
    if ($data[-1]){ $data[2] -= 3 } else { $data[1] += 3 }
    push @res, \@data
  }
  $sth3->finish();
  return \@res
}
