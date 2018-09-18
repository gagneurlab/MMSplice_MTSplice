=head1 LICENSE

Copyright [1999-2015] Wellcome EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

 MMSplice

=head1 SYNOPSIS

 mv MMSplice.pm ~/.vep/Plugins
 pip install mmsplice
 ./vep -i variants.vcf --plugin MMSplice,[port_of_mmsplice_server=5000],[intronl_len=100],[intronr_len=80],[exon_cut_l=0],[exon_cut_r=0],[acceptor_intron_cut=6],[donor_intron_cut=3],[acceptor_intron_len=20],[acceptor_exon_len=3],[donor_exon_len=3],[donor_intron_len=6],[acceptor_intronM],[acceptorModelFile],[exonModelFile],[donorModelFile],[donor_intronModelFile]


=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 runs MMSplice.

Add more doc.

=cut

package MMSplice;

use strict;
use warnings;

use List::Util qw(max);
use LWP::UserAgent;

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

our $CACHE_SIZE = 50;

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  $self->init_params();
  $self->init_api();

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
    return {
      ref_acceptor_intron => "acceptor intron ref",
      ref_acceptor => "acceptor ref",
      ref_exon => "exon ref",
      ref_donor => "donor ref",
      ref_donor_intron => "donor intron ref",
      alt_acceptor_intron => "acceptor intron alt",
      alt_acceptor => "acceptor alt",
      alt_exon => "exon alt",
      alt_donor => "alt donor",
      alt_donor_intron => "alt donor intron"
    };
}

sub init_params {
  my $self = shift;
  my $params = $self->params;

  $self->{api_port} = shift @$params || 5000;

  $self->{overhang_l} = shift @$params || 100;
  $self->{overhang_r} = shift @$params || 80;
  $self->{exon_cut_l} = shift @$params || 0;
  $self->{exon_cut_r} = shift @$params || 0;
  $self->{acceptor_intron_cut} = shift @$params || 6;
  $self->{donor_intron_cut} = shift @$params || 6;
  $self->{acceptor_intron_len} = shift @$params || 50;
  $self->{acceptor_exon_len} = shift @$params || 3;
  $self->{donor_exon_len} = shift @$params || 5;
  $self->{donor_intron_len} = shift @$params || 13;

  $self->{acceptor_intronM} = shift @$params || "";
  $self->{acceptorM} = shift @$params || "";
  $self->{exonM} = shift @$params || "";
  $self->{donorM} = shift @$params || "";
  $self->{donor_intronM} = shift @$params || "";
}

sub init_api {
  my $self = shift;

  $self->{ua} = LWP::UserAgent->new;
  $self->{server_endpoint} = "http://localhost:$self->{api_port}";

  $self->start_api();
  $self->create_model();
}

sub start_api {
  my $self = shift;

  my $pid = fork();
  die if not defined $pid;

  if ($pid == 0) {
      exec("mmsplice run_api --port=$self->{api_port}");
  }
  else {
      $self->{api_pid} = $pid;
  }

  sleep(10);
}

sub create_model {
  my ($self, $seq) = @_;

  my $req = HTTP::Request->new(POST => $self->{server_endpoint} . "/create-model");
  $req->header('content-type' => 'application/json');

  my $content = qq{{
  "acceptor_intronM": "$self->{acceptor_intronM}",
  "acceptorM": "$self->{acceptorM}",
  "exonM": "$self->{exonM}",
  "donorM": "$self->{donorM}",
  "donor_intronM": "$self->{donor_intronM}",

  "exon_cut_l": $self->{exon_cut_l},
  "exon_cut_r": $self->{exon_cut_r},
  "acceptor_intron_cut": $self->{acceptor_intron_cut},
  "donor_intron_cut": $self->{donor_intron_cut},
  "acceptor_intron_len": $self->{acceptor_intron_len},
  "acceptor_exon_len": $self->{acceptor_exon_len},
  "donor_exon_len": $self->{donor_exon_len},
  "donor_intron_len": $self->{donor_intron_len}
}};

  $req->content($content);

  my $resp = $self->{ua}->request($req);
  if ($resp->is_success) {
      return;
  }
  else {
      print "HTTP GET error code: ", $resp->code, "\n";
      print "HTTP GET error message: ", $resp->message, "\n";
  }
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  my $tv = $tva->transcript_variation;
  my $tr = $tva->transcript;
  my $tr_strand = $tr->strand;

  foreach my $exon(@{$self->overlap_exons($tr, $vf)}) {

    my ($splicing_start, $splicing_end) = @{$self->overhanged_interval($exon, $tr_strand)};
    my $ref_seq = $self->fetch_seq($tva, $splicing_start, $splicing_end);
    my $alt_seq = $self->fetch_variant_seq($tva, $exon, $splicing_start, $splicing_end);

    my @ref_scores = $self->req_psi_score($ref_seq);
    my @alt_scores = $self->req_psi_score($alt_seq);

    return {
      ref_acceptor_intron => $ref_scores[0],
      ref_acceptor => $ref_scores[1],
      ref_exon => $ref_scores[2],
      ref_donor => $ref_scores[3],
      ref_donor_intron => $ref_scores[4],
      alt_acceptor_intron => $alt_scores[0],
      alt_acceptor => $alt_scores[1],
      alt_exon => $alt_scores[2],
      alt_donor => $alt_scores[3],
      alt_donor_intron => $alt_scores[4]
    }
  }

  return {};
}

sub variant_ref {
  my ($self, $tva) = @_;
  my $vf = $tva->variation_feature;

  if ($vf->ref_allele_string eq '-') {
    return '';
  }
  return $vf->ref_allele_string;
}

sub variant_alt {
  my ($self, $tva) = @_;

  if ($tva->feature_seq eq '-') {
    return '';
  }
  return $tva->feature_seq;
}

sub overhanged_interval {
  my ($self, $exon, $tr_strand) = @_;

  my ($splicing_start, $splicing_end);

  if($tr_strand > 0) {
    $splicing_start = $exon->start - $self->{overhang_l};
    $splicing_end = $exon->end + $self->{overhang_r};
  }
  else {
    $splicing_start = $exon->start - $self->{overhang_r};
    $splicing_end = $exon->end + $self->{overhang_l};
  }

  return [$splicing_start, $splicing_end];
}

sub overlap_exons {
    my ($self, $tr, $vf) = @_;

    my @exons;
    my ($vf_start, $vf_end) = ($vf->start, $vf->end);

    foreach my $exon(@{$tr->get_all_Exons()}) {

      my ($low, $high) = ($exon->{start}, $exon->{end});
      my ($start, $end) = ($vf_start - $self->{overhang_l}, $vf_end + $self->{overhang_r});

      unless ($high < $start or $low > $end) {
          push @exons, $exon;
      }
    }

    return \@exons;
}

sub len_diff {
  my ($self, $tva) = @_;
  return (length $self->variant_alt($tva)) - (length $self->variant_ref($tva));
}

sub variant_side {
  my ($self, $vf, $tr_strand, $exon) = @_;

  if ($tr_strand > 0) {
    if ($vf->{start} < $exon->start) {
      return "5'";
    }
    elsif ($vf->{start} > $exon->end) {
      return "3'";
    }
  }
  else {
    if ($vf->{start} < $exon->start) {
      return "3'";
    }
    elsif ($vf->{start} > $exon->end) {
      return "5'";
    }
  }

  return "n'";
}

sub fetch_variant_seq {
  my ($self, $tva, $exon, $splicing_start, $splicing_end) = @_;
  my $vf = $tva->variation_feature;
  my $tr_strand = $tva->transcript->strand;

  my ($ibefore_start, $ibefore_end) = ($splicing_start, $vf->{start} - 1);
  my ($iafter_start, $iafter_end) = ($vf->{end} + 1, $splicing_end);

  my $v_side = $self->variant_side($vf, $tr_strand, $exon);

  if (($tr_strand > 0 && $v_side eq "3'") || ($tr_strand < 0 && $v_side eq "5'")) {
    $iafter_end -= $self->len_diff($tva);
  }
  elsif (($tr_strand > 0 && $v_side eq "5'") || ($tr_strand < 0 && $v_side eq "3'")) {
    $ibefore_start += $self->len_diff($tva);
  }

  my $before_seq = $self->fetch_seq($tva, $ibefore_start, $ibefore_end);
  my $after_seq = $self->fetch_seq($tva, $iafter_start, $iafter_end);

  if ($tr_strand < 0){
    ($before_seq, $after_seq) = ($after_seq, $before_seq);
  }

  return $before_seq . $self->variant_alt($tva) . $after_seq;
}

sub req_psi_score {
  my ($self, $seq) = @_;

  my $req = HTTP::Request->new(POST => $self->{server_endpoint} . "/psi-score");
  $req->header('content-type' => 'application/json');

  my $content = qq{{
  "intronl_len": $self->{overhang_l},
  "intronr_len": $self->{overhang_r},
  "seq": "$seq"
}};

  $req->content($content);

  my $resp = $self->{ua}->request($req);
  if ($resp->is_success) {

      my $decoded_content = $resp->decoded_content;

      my @scores = split(',', $decoded_content);

      return @scores;
  }
  else {
      print "HTTP POST error code: ", $resp->code, "\n";
      print "HTTP POST error message: ", $resp->message, "\n";
  }
}

sub fetch_seq {
  my ($self, $tva, $start, $end) = @_;
  my $vf = $tva->variation_feature;
  my $tr_strand = $tva->transcript->strand;

  my $seq = $vf->{slice}->sub_Slice(
    $start,
    $end,
    $tr_strand
  )->seq;

  return $seq;
}

sub DESTROY {
    my $self = shift;

    kill 2, $self->{api_pid};
    print('api closed');
}

1;
