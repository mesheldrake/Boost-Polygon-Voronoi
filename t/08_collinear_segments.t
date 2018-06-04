#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 24;
use Boost::Polygon::Voronoi qw(calc_Q tangent_angle);

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

    make_result_svg(\@edges,\@cells,undef,5.0,'medial axis');
    #make_result_svg(\@edges,\@cells,undef,5.0,'Voronoi diagram');
    
    my $finite_count = scalar(grep $_->is_finite(), @edges);
    ok($finite_count == 26, "26 finite edges? ($finite_count)");

    $finite_count = scalar(grep $_->is_finite() && $_->is_secondary(), @edges);
    ok($finite_count == 0, "0 finite secondary edges? ($finite_count)");

    $finite_count = scalar(grep $_->is_finite() && $_->is_primary(), @edges);
    ok($finite_count == 26, "26 finite primary edges? ($finite_count)");

    $finite_count = scalar(grep $_->is_finite() && $_->is_primary() && $_->is_curved(), @edges);
    ok($finite_count == 0, "0 finite primary curved edges? ($finite_count)");

    my @internal = grep $_->is_internal, @edges;

    # main symptom of collinear input segments was that you
    # would get what looked like a corner, and would be labeled as
    # "external", but it would be on the main center line of the medial
    # axis, in a sequence of edges labeled "internal". So we loop
    # over what should be all internal edge sequences, looking for
    # any marked external.
    # Fix for this was to add a vertex and relink these corner edges
    # (in the C++ code) so that collinear input segments get "flat corners" 
    # at all those collinear segment meeting points, which should work 
    # fine with any other code that looks for or deals with corners.

    my $limit_fail='';
    my $external_fail='';
    my %seen_edges;

    foreach my $edge (@internal) {
        my $external_count=0;
        my $start_e = $edge->next()->prev();
        next if $seen_edges{''.($$start_e)};
        my $e = $edge;
        my $limit=100; # arbitrary.
        do {
            $seen_edges{''.($$e)}=1;
            $e = $e->next();
            if ($e->is_external()) {$external_count++;}
        } while (($$e) != ($$start_e) && $limit-- > 0);

        if (!($limit > 0)) {
            $limit_fail='(reached edge looping safety limit)';
            warn "start edge that led to too much looping:\n",
            '  ['.($$start_e).']',' n['.(${$start_e->next()}).']',' p['.(${$start_e->prev()}).']',"\n",
            'np['.(${$start_e->next()->prev()}).']',"\n",
            'pn['.(${$start_e->prev()->next()}).']',"\n",
            $start_e->is_infinite() ? 'INF':'FIN',' ',
            $start_e->is_primary() ? 'PRI':'2ND',' ',
            $start_e->is_internal() ? 'INT':'EXT',' ',
            (!($start_e->vertex0() && $start_e->vertex1()) && $start_e->vertex0()) ? 'out':'in',"\n",
            $start_e->vertex0() ? 'v0['.$start_e->vertex0()->x().', '.$start_e->vertex0()->y().']':'v0[]',"\n",
            $start_e->vertex1() ? 'v1['.$start_e->vertex1()->x().', '.$start_e->vertex1()->y().']':'v1[]',"\n"
            ;
        }
        if ($external_count) {
            $external_fail='(found external edge in an internal cell loop)';
        }
    }
    ok(!$limit_fail, "Internal edge cell loops are closed loops. $limit_fail");
    ok(!$external_fail, "No external edges in internal cell edge loops. $external_fail");


}

__END__
