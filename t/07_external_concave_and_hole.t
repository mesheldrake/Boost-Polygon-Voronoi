#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 21;
use Boost::Polygon::Voronoi qw(calc_Q tangent_angle get_curve_focus);

eval {require "t/svg_test_helpers.pl"}; 
if ($@) {eval{require "svg_test_helpers.pl"}}
if ($@) {fail("couldn't find svg_test_helpers.pl");}

{

    my $b = new Boost::Polygon::Voronoi::builder;

    load_test_svg($b);

    my $ma = get_medial_axis($b);
    #my $ma = get_voronoi_diagram($b);

    my @edges = @{$ma->edges()};
    my @cells = @{$ma->cells()};

    make_result_svg(\@edges,\@cells,undef,1.0,'medial axis');
    #make_result_svg(\@edges,\@cells,undef,1.0,'Voronoi diagram');

    my @external = grep $_->is_external, @edges;
    my @internal = grep $_->is_internal, @edges;

    my $external_count = scalar(@external);
    my $internal_count = scalar(@internal);
    

    # actually these aren't final
    #warn "( $external_count  and $internal_count )";
    # thought I was doing the 02 test here
    # figure equiv for that, and fix these when you get this looking right
    ok($external_count == 212, "$external_count of 284 external edges");
    ok($internal_count == 262, "$internal_count of 190 internal edges");


    my @vertices = @{$ma->vertices()};

    my @internal_from_verts = grep $_->incident_edge()->is_internal(), @vertices;
    my @internal_infinite_from_verts = grep $_->incident_edge()->is_infinite(), @internal_from_verts;
    ok(scalar(@internal_infinite_from_verts) == 0, "internal edges cannot be infinite (from verts)");


}

__END__
