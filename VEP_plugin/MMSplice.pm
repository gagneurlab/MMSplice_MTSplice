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
 ./vep -i variants.vcf --plugin MMSplice,[intronl_len=100],[intronr_len=100],[exon_cut_l=0],[exon_cut_r=0],[acceptor_intron_cut=6],[donor_intron_cut=3],[acceptor_intron_len=20],[acceptor_exon_len=3],[donor_exon_len=3],[donor_intron_len=6],[acceptor_intronM],[acceptorModelFile],[exonModelFile],[donorModelFile],[donor_intronModelFile]


=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 runs MMSplice (modular modeling of splicing) which performs a set of prediction
 on splicing.

 The plugin requires MMSplice python package as an external dependency since it wraps mmsplice package as vep plugin.
 Thus, MMSplice package should be installed with `pip install mmsplice`.
 Then, it automatically runs python server in background and analysis variant with python server.

 The plugin predicts delta_logit_psi and pathogenicity values of variants in addition to the score of each component, for both reference and variant sequences, such as acceptor_intron, acceptor, exon, donor, and donor_intron.

 The plugin don't filters any variant. Some of the variants may not have prediction because they are not matched. In this case, emtpy values are returned.
=cut

package MMSplice;

use strict;
use warnings;
use diagnostics;
use IPC::Open3;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);

    $self->init_params();
    $self->init_python();

    return $self;
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        mmsplice_ref_acceptor_intron => "acceptor intron score of reference sequence",
        mmsplice_ref_acceptor => "acceptor score of reference sequence",
        mmsplice_ref_exon => "exon score of reference sequence",
        mmsplice_ref_donor => "donor score of reference sequence",
        mmsplice_ref_donor_intron => "donor intron score of reference sequence",
        mmsplice_alt_acceptor_intron => "acceptor intron score of variant sequence ",
        mmsplice_alt_acceptor => "acceptor score of variant sequence",
        mmsplice_alt_exon => "exon score of variant sequence",
        mmsplice_alt_donor => "donor score of variant sequence",
        mmsplice_alt_donor_intron => "donor intron score of variant sequence",
        mmsplice_delta_logit_psi => "delta logit psi score of variant",
        mmsplice_pathogenicity => "pathogenicity effect of variant"
    };
}

sub init_params {
    my $self = shift;
    my $params = $self->params;

    $self->{overhang_l} = shift @$params || 100;
    $self->{overhang_r} = shift @$params || 100;
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

sub call_python {
    my ($self, $content, $timeout) = @_;
    $timeout = $timeout || 1;
    my $python_stdout = $self->{python_stdout};
    my $python_stdin = $self->{python_stdin};
    my $python_selout = $self->{python_selout};
    my $response_keyword = "MMSPLICE-RESPONSE:";

    print $python_stdin "$content\n";

    my $result = '';
    while($python_selout->can_read($timeout)) {

        chomp($result = <$python_stdout>);

        if ($result eq "")
        {
            next;
        }

        if(substr($result, 0, length($response_keyword)) eq $response_keyword) {
            $result = substr($result, length($response_keyword), length($result));
            last;
        }

        print "$result\n";
    }
    return $result;
}

sub init_python {
    my $self = shift;

    $self->{api_pid} = open3(my $python_stdin, my $python_stdout,  my $python_stderr, "mmsplice run");

    my $python_selout = new IO::Select();
    $python_selout->add($python_stdout);
    my $python_selerr = new IO::Select();
    $python_selerr->add($python_stderr);

    $self->{python_selout} = $python_selout;
    $self->{python_selerr} = $python_selerr;
    $self->{python_stdin} = $python_stdin;
    $self->{python_stdout} = $python_stdout;
    $self->{python_stderr} = $python_stderr;

    my $content = qq[{"acceptor_intronM": "$self->{acceptor_intronM}", "acceptorM": "$self->{acceptorM}", "exonM": "$self->{exonM}", "donorM": "$self->{donorM}", "donor_intronM": "$self->{donor_intronM}", "exon_cut_l": $self->{exon_cut_l}, "exon_cut_r": $self->{exon_cut_r}, "acceptor_intron_cut": $self->{acceptor_intron_cut}, "donor_intron_cut": $self->{donor_intron_cut}, "acceptor_intron_len": $self->{acceptor_intron_len}, "acceptor_exon_len": $self->{acceptor_exon_len}, "donor_exon_len": $self->{donor_exon_len}, "donor_intron_len": $self->{donor_intron_len}}];

    my $status = $self->call_python($content, 60);
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

        my @scores = $self->get_psi_score($ref_seq, $alt_seq);

        return {
            mmsplice_ref_acceptor_intron => $scores[0],
            mmsplice_ref_acceptor => $scores[1],
            mmsplice_ref_exon => $scores[2],
            mmsplice_ref_donor => $scores[3],
            mmsplice_ref_donor_intron => $scores[4],
            mmsplice_alt_acceptor_intron => $scores[5],
            mmsplice_alt_acceptor => $scores[6],
            mmsplice_alt_exon => $scores[7],
            mmsplice_alt_donor => $scores[8],
            mmsplice_alt_donor_intron => $scores[9],
            mmsplice_delta_logit_psi => $scores[10],
            mmsplice_pathogenicity => $scores[11]
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

sub get_psi_score {
    my ($self, $ref_seq, $alt_seq) = @_;

    my $content = qq[{"intronl_len": $self->{overhang_l}, "intronr_len": $self->{overhang_r}, "ref_seq": "$ref_seq", "alt_seq": "$alt_seq"}];
    my $response = $self->call_python($content);

    my @scores = split(',', $response);
    return @scores;
}

sub fetch_seq {
    my ($self, $tva, $start, $end) = @_;
    my $vf = $tva->variation_feature;
    my $tr_strand = $tva->transcript->strand;

    my $seq = $vf->{slice}->sub_Slice($start, $end, $tr_strand)->seq;
    return $seq;
}

sub DESTROY {
    my $self = shift;
    kill 2, $self->{api_pid} if (defined $self->{api_pid});
}

1;
