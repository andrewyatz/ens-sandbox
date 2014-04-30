#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=pod

=head1 CHAIN FILES

This code writes chain files from Ensembl assembly records. Chain files are defined at

  http://genome.ucsc.edu/goldenPath/help/chain.html

They can be used with liftover (http://genome.ucsc.edu/cgi-bin/hgLiftOver) and crossmap (http://crossmap.sourceforge.net). The file format says the following:

  chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id 
  size dt dq
  size

Majority of fields are self explanatory expect for the declaration of target & query sequences. In a mapping of NCBI36 (hg18) to GRCh37 (hg19) the following roles are designated:

  - target == NCBI36
  - query  == GRCh37

You can see this in the the example chain files for human mappings http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz and http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz.

dt & dq are meant to be the differences between the end of this current block and the start of the next in the target and query. The final line is just a size to indicate no further offset.

tSize and qSize are the total sizes of the target and reference sequences.

=cut

use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::IO qw/work_with_file/;
use feature qw/say/;
use File::Path qw/mkpath/;
use File::Spec;

my $chain_id = 1;
my %global_mapings;

sub get_options {
  my ($db_name, $db_host, $db_user, $db_pass, $db_port, $help, $species, $group, $release, @external_db);
  # my ($asm, $cmp) = ('GRCh37', 'NCBI36');
  $db_port = 3306;
  $species = 'human';
  $group = 'core';

  GetOptions(
    "db_name|dbname|database=s"       => \$db_name,
    "db_host|dbhost|host=s"           => \$db_host,
    "db_user|dbuser|user|username=s"  => \$db_user,
    "db_pass|dbpass|pass|password=s"  => \$db_pass,
    "db_port|dbport|port=s"           => \$db_port,
    "species=s"                       => \$species,
    "version|release=i"               => \$release,
    "external_db=s@"                   => \@external_db,
    "h!"                              => \$help,
    "help!"                           => \$help,
  );

  die "No -external_db given" unless @external_db;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -HOST => $db_host, -PORT => $db_port, 
    -USER => $db_user, -PASS => $db_pass,
    -DB_VERSION => $release,
  );
  my $core_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
  return ($core_dba, @external_db);
}

run();

sub run {
  my ($core_dba, @external_dbs) = get_options();
  my $prod_name = $core_dba->get_MetaContainer->get_production_name();

  say "Working with ${prod_name}";
  say "Fetching slices";
  my $slices = slices($core_dba);
  say '';

  foreach my $external_db (@external_dbs) {
    say "Processing ${external_db}";

    if(! should_run($core_dba, $external_db)) {
      say "\tNo synonyms found. Nothing to convert";
      return;
    }

    write_mappings($prod_name, $slices, $external_db, 1);
    write_mappings($prod_name, $slices, $external_db, 0);
    say "\tDone";
    say q{};
  }
  write_tab_lookup($core_dba, $prod_name, @external_dbs);
}

sub write_tab_lookup {
  my ($core_dba, $prod_name, @external_dbs) = @_;

  my $assembly = $core_dba->get_GenomeContainer->get_version();
  my $file = "${assembly}_region_lookup.tsv";
  my $dir = File::Spec->catdir(File::Spec->curdir(), $prod_name);
  mkpath($dir);
  my $path = File::Spec->catfile($dir, $file);

  say "Writing tab lookup to $path";
  work_with_file($path, 'w', sub {
    my ($fh) = @_;
    say $fh '#'.join("\t", ('Ensembl', @external_dbs));
    foreach my $slice (sort keys %global_mapings) {
      say $fh join("\t", $slice, map { defined($_) ? $_ : '.' } map {$global_mapings{$slice}{$_}} @external_dbs);
    }
  });
  return;
}

sub should_run {
  my ($core_dba, $external_db) = @_;
  my $sql = <<'SQL';
select count(*) 
from seq_region 
join seq_region_synonym using(seq_region_id) 
join external_db using (external_db_id) 
join coord_system using (coord_system_id)
where db_name = ?
and species_id =?
SQL
  my $params = [$external_db, $core_dba->species_id()];
  return $core_dba->dbc->sql_helper->execute_single_result(-SQL => $sql, -PARAMS => $params);
}

sub write_mappings {
  my ($prod_name, $slices, $external_db, $to_ensembl) = @_;
  
  my $direction = ($to_ensembl) ? "from ${external_db} to Ensembl" : "from Ensembl to ${external_db}";
  say "\tBuilding chain mappings $direction";
  my $chains = build_chains($slices, $external_db, $to_ensembl);

  # Setup path
  my $file = ($to_ensembl) ? "${external_db}ToEnsembl.chain" : "EnsemblTo${external_db}.chain";
  my $dir = File::Spec->catdir(File::Spec->curdir(), $prod_name);
  mkpath($dir);
  my $path = File::Spec->catfile($dir, $file);

  say "\tWriting mappings to $path";
  work_with_file($path, 'w', sub {
    my ($fh) = @_;
    print_chains($fh, $chains);
  });
  return;
}

sub slices {
  my ($core_dba) = @_;
  my %hash; 
  my $slices = $core_dba->get_SliceAdaptor->fetch_all('toplevel', undef, 'nonref');
  while(my $slice = shift @{$slices}) {
    $hash{$slice->seq_region_name()} = $slice;
  }
  return [values %hash];
}

sub build_chains {
  my ($slices, $external_db, $to_ensembl) = @_;
  my @chains;
  foreach my $slice (@{$slices}) {
    my $chain = build_chain($slice, $external_db, $to_ensembl);
    next unless defined $chain;
    push(@chains, $chain);
  }
  return \@chains;
}

# trg is the source data & q is the target i.e. converting from UCSC to Ensembl means trg == chr1 & q == 1
sub build_chain {
  my ($slice, $external_db, $to_ensembl) = @_;
  my ($synonym) = @{$slice->get_all_synonyms($external_db)};
  my $seq_region_name = $slice->seq_region_name();
  return unless $synonym;
  $global_mapings{$seq_region_name}{$external_db} = $synonym->name();
  my ($t_name, $q_name) = ($to_ensembl) ? ($synonym->name(), $seq_region_name) : ($seq_region_name, $synonym->name());
  my $size = $slice->seq_region_length();
  return {
    # header => ['chain', $chain_score, $t_name, $t_size, $t_strand, $t_start, $t_end, $q_name, $q_size, $q_strand, $q_start, $q_end, $chain_id],
    header => ['chain', 1, $t_name, $size, '+', 0, $size, $q_name, $size, '+', 0, $size, $chain_id++],
    # gaps => [[match t_gap_to_next_block q_gap_to_next_block]]
    gaps => [[$size]]
  };
}

sub print_chains {
  my ($fh, $chains) = @_;
  while(my $chain_def = shift @{$chains}) {
    my $header = $chain_def->{header};
    my $gaps = $chain_def->{gaps};
    print $fh join(q{ }, @{$header});
    print $fh "\n";
    foreach my $gap (@{$gaps}) {
      print $fh join(q{ }, @{$gap});
      print $fh "\n";
    }
    print $fh "\n";
  }
}
