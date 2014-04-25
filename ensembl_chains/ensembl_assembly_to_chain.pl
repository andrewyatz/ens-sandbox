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
use Scalar::Util qw/looks_like_number/;

sub get_options {
  my ($db_name, $db_host, $db_user, $db_pass, $db_port, $help, $species, $group, $release);
  my ($asm, $cmp) = ('GRCh37', 'NCBI36');
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
    "asm=s"                           => \$asm,
    "cmps=s"                          => \$cmp,
    "h!"                              => \$help,
    "help!"                           => \$help,
  );

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -HOST => $db_host, -PORT => $db_port, 
    -USER => $db_user, -PASS => $db_pass,
    -DB_VERSION => $release,
  );
  my $core_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
  return ($core_dba, $asm, $cmp);
}

run();

sub run { 
  my ($core_dba, $asm_cs, $cmp_cs) = get_options();
  my $asm_to_cmp_mappings = get_assembly_mappings($core_dba, $asm_cs, $cmp_cs);
  write_mappings($asm_cs, $cmp_cs, $asm_to_cmp_mappings);
  my $cmp_to_asm_mappings = get_reverse_assembly_mappings($core_dba, $asm_cs, $cmp_cs);
  write_mappings($cmp_cs, $asm_cs, $cmp_to_asm_mappings);
}

sub write_mappings {
  my ($source_cs, $target_cs, $mappings) = @_;
  my $file = "${source_cs}To${target_cs}.chain";
  my $chains = build_chain_mappings($mappings);
  work_with_file($file, 'w', sub {
    my ($fh) = @_;
    print_chains($fh, $chains);
  });
  return;
}

sub build_chain_mappings {
  my ($assembly_mappings) = @_;
  my @chain_mappings;
  
  my ($t_name, $t_size, $t_strand, $t_start, $t_end);
  my ($q_name, $q_size, $q_strand, $q_start, $q_end);
  my $chain_id = 1;
  my @chain_gaps;
  
  my $length = scalar(@{$assembly_mappings});
  for (my $i = 0; $i < $length; $i++) {

    my $current = $assembly_mappings->[$i];
    my $next = ($i+1 != $length) ? $assembly_mappings->[$i+1] : undef;

    my $ori = $current->{ori};
    my ($asm_diff, $cmp_diff);
    if($next) {
      $asm_diff = ($next->{asm_start} - $current->{asm_end})-1;
      # Rev strands means the next cmp region has a lower start than the 
      # current end (because it's running in reverse). Rember length in 1-based
      # coords always is (high-low)-1
      $cmp_diff = ($ori == 1) ? ($next->{cmp_start} - $current->{cmp_end})-1 : ($current->{cmp_start} - $next->{cmp_end})-1;
    }

    if(! $t_name) {
      # Reset variables to current
      @chain_gaps = ();
      $chain_id++;
      ($t_name, $t_size, $t_strand, $t_start) = ($current->{asm_name}, $current->{asm_length}, 1, $current->{asm_start});
      ($q_name, $q_size, $q_strand) = ($current->{cmp_name}, $current->{cmp_length}, $current->{ori});
      $q_start = ($ori == 1) ? $current->{cmp_start} : $current->{cmp_end};
    }

    # Can mean we are into a new chromsome, strand has swapped or we have run out of mappings
    if( ! defined $next || 
        $t_name ne $next->{asm_name} || 
        $ori != $next->{ori} ||
        $cmp_diff < 0) {
      # $DB::single = 1;;
      # Add the last gap on which is just the length of this alignment
      push(@chain_gaps, [$current->{length}]);
      # Set the ends of the chain since this is the last block
      $t_end = $current->{asm_end};
      $q_end = ($ori == 1) ? $current->{cmp_end} : $current->{cmp_start};

      if ($i != 0) {
        #If strand was negative we need to represent all data as reverse complemented regions
        if($q_strand == -1) {
          # $t_start = ($t_size - $t_start)+1;
          # $t_end = ($t_size - $t_end)+1;
          $q_start = ($q_size - $q_start)+1;
          $q_end = ($q_size - $q_end)+1;
        }
        # Convert to UCSC formats (0-based half-open intervals and +/- strands)
        $t_start--;
        $q_start--;
        $t_strand = ($t_strand == 1) ? '+' : '-';
        $q_strand = ($q_strand == 1) ? '+' : '-';

        #Store the chain
        my $chain_score = 1;
        push(@chain_mappings, {
          header => ['chain', $chain_score, $t_name, $t_size, $t_strand, $t_start, $t_end, $q_name, $q_size, $q_strand, $q_start, $q_end, $chain_id],
          gaps => [@chain_gaps]
        });
      }

      if(! defined $next) {
        last;
      }

      # Clear variables
      ($t_name, $t_size, $t_strand, $t_start) = ();
      ($q_name, $q_size, $q_strand, $q_start) = ();
    }

    push(@chain_gaps, [$current->{length}, $asm_diff, $cmp_diff]);
  }

  return \@chain_mappings;
}

sub get_assembly_mappings {
  my ($dba, $asm_version, $cmp_version) = @_;
  my $sql = get_sql(
    ['sr1.name as asm_name', 
    'sr1.length as asm_length', 
    'sr2.name as cmp_name', 
    'sr2.length as cmp_length', 
    'asm.asm_start', 'asm.asm_end', 'asm.cmp_start', 'asm.cmp_end'],
    'order by cs2.version, sr2.name, asm.cmp_start');
  return $dba->dbc->sql_helper->execute( -SQL => $sql, -PARAMS => [$asm_version, $cmp_version], -USE_HASHREFS => 1);
}

sub get_reverse_assembly_mappings {
  my ($dba, $asm_version, $cmp_version) = @_;
  # Reverse the columns to get the reverse mapping
  my $sql = get_sql(
  ['sr1.name as cmp_name', 
  'sr1.length as cmp_length', 
  'sr2.name as asm_name', 
  'sr2.length as asm_length', 
  'asm.cmp_start as asm_start', 'asm.cmp_end as asm_end', 'asm.asm_start as cmp_start', 'asm.asm_end as cmp_end'],
  'order by cs2.version, sr2.name, asm.cmp_start');
  return $dba->dbc->sql_helper->execute( -SQL => $sql, -PARAMS => [$asm_version, $cmp_version], -USE_HASHREFS => 1);
}

sub get_sql {
  my ($columns, $order_by) = @_;
  my $select = join(q{, }, @{$columns}, 'asm.ori', '(asm.asm_end - asm.asm_start)+1 as length');
  return <<SQL;
select ${select}
from coord_system cs1
join seq_region sr1 on (cs1.coord_system_id = sr1.coord_system_id)
join assembly asm on (sr1.seq_region_id = asm.asm_seq_region_id)
join seq_region sr2 on (sr2.seq_region_id = asm.cmp_seq_region_id)
join coord_system cs2 on (sr2.coord_system_id = cs2.coord_system_id)
where cs1.version =?
and cs2.version =?
$order_by
SQL
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