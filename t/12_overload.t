#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 109;
use Boost::Polygon::Voronoi qw(calc_Q tangent_angle);

eval {
    require "t/svg_test_helpers.pl";
    }; 
if ($@) {eval{
    require "svg_test_helpers.pl";
};}
ok(! $@ , "loaded test helpers");


# test whether overloaded equality and bool operators
# work for edges, vertices, cells

# using scopy of test shape for test 02_rectangle.t

foreach my $vd_or_ma ('Voronoi diagram','medial axis') {
    
    my $tpre = $vd_or_ma=~/Voronoi/ ? 'VD:':'MA:';

    my $b = new Boost::Polygon::Voronoi::builder;

    load_test_svg($b);

    my $vd = $vd_or_ma eq 'Voronoi diagram'
           ? get_voronoi_diagram($b) 
           : get_medial_axis($b);


    my @edges = @{$vd->edges()};
    my @cells = @{$vd->cells()};
    my @vertices = @{$vd->vertices()};
    
    my $e = $edges[0];
    my $c = $cells[0];
    my $v = (grep {$_->is_finite() && $_->is_primary} @edges)[0]->vertex0();

    # first see that things check out without triggering overloaded operators
    ok(defined $e, "$tpre an edge exists and is defined");
    
    my $edge_class_expected = "Boost::Polygon::".($vd_or_ma=~/Voronoi/?'Voronoi':'MedialAxis')."::Edge";
    my $edge_class = ref($e);
    my $edge_class_match = $e->isa($edge_class_expected);
    ok($edge_class_match,"$tpre edge isa $edge_class_expected"
        . ($edge_class_match ? '':" (".$edge_class.")"));
    
    ok(defined $e->twin(),"$tpre an edge has a twin");
    
    ok(${$e->twin()->twin()} == ${$e},
        "$tpre a deref'd edge is the same as it's twin's twin, deref'd");
    
    ok(!(${$e->twin()} == ${$e}),
        "$tpre a deref'd edge is NOT the same as it's twin, deref'd");
    
    # now exercise the overloaded operators
    
    ok((1 && $e),"$tpre an edge is boolean-true");
    
    ok((1 && $e->twin()),"$tpre an edge's twin is boolean-true");
    
    # these two work whether or not they really work, don't they?
    ok(!($e->twin() == $e),
        "$tpre edge == edge->twin is false");    
    ok($e->twin() != $e,
        "$tpre edge != edge->twin is true");

    # this is the key one to test
    ok($e->twin()->twin() == $e,
        "$tpre edge == edge->twin->twin is true");

    # now make sure the overload is inherited when we upgrade to
    # auto-scaling class
    
    bless $e, $edge_class_expected."Scaled";
    ok($e->isa($edge_class_expected."Scaled"), 
        "$tpre subclassed edge");
    ok($e->twin()->twin() == $e,
        "$tpre edge == edge->twin->twin is still true");

    
    
    
    # now similar for cells
    
    # first see that things check out without triggering overloaded operators
    ok(defined $c, "$tpre a cell exists and is defined");
    my $c_same = $c->incident_edge()->cell();
    my $c_diff = $c->incident_edge()->twin->cell();
    ok(defined $c_same,
        "$tpre a cell retrieved through references is defined");
    ok(defined $c_diff,
        "$tpre a different cell retrieved through references is defined");
    
    my $cell_class_expected = "Boost::Polygon::".($vd_or_ma=~/Voronoi/?'Voronoi':'MedialAxis')."::Cell";
    my $cell_class = ref($c);
    my $cell_class_match = $c->isa($cell_class_expected);
    ok($cell_class_match,"$tpre cell isa $cell_class_expected"
        . ($cell_class_match ? '':" (".$cell_class.")"));
    
    ok(${$c_same} == ${$c},
        "$tpre two references to the same cell are == when deref'd");
    
    ok(!(${$c_diff} == ${$c}),
        "$tpre two references to different cells fail to be == when deref'd");
    
    # now exercise the overloaded operators
    
    ok((1 && $c),"$tpre a cell is boolean-true");
    ok((1 && $c_same),"$tpre a different ref to the same cell is boolean-true");
    ok((1 && $c_diff),"$tpre a different cell is boolean-true");
    
    # these two work whether or not they really work, don't they?
    ok(!($c_diff == $c),
        "$tpre cell == different cell is false");    
    ok($c_diff != $c,
        "$tpre cell != different cell is true");

    # this is the key one to test
    ok($c_same == $c,
        "$tpre cell == same cell from different ref is true");

    # now make sure the overload is inherited when we upgrade to
    # auto-scaling class
    
    bless $c, $cell_class_expected."Scaled";
    ok($c->isa($cell_class_expected."Scaled"), 
        "$tpre subclassed cell");
    ok($c_same == $c,
        "$tpre subclassed cell == same cell from different ref with original class");
    
    
    # now similar for vertices
    
    # first see that things check out without triggering overloaded operators
    ok(defined $v, "$tpre a vertex exists and is defined");
    my $v_same = $v->incident_edge()->vertex0();
    my $v_diff = $v->incident_edge()->vertex1();
    ok(defined $v_same,
        "$tpre a vertex retrieved through references is defined");
    ok(defined $v_diff,
        "$tpre a different vertex retrieved through references is defined");
    
    my $vertex_class_expected = "Boost::Polygon::".($vd_or_ma=~/Voronoi/?'Voronoi':'MedialAxis')."::Vertex";
    my $vertex_class = ref($v);
    my $vertex_class_match = $v->isa($vertex_class_expected);
    ok($vertex_class_match,"$tpre vertex isa $vertex_class_expected"
        . ($vertex_class_match ? '':" (".$vertex_class.")"));
    
    ok(${$v_same} == ${$v},
        "$tpre two references to the same vertex are == when deref'd");
    
    ok(!(${$v_diff} == ${$v}),
        "$tpre two references to different verticess fail to be == when deref'd");
    
    # now exercise the overloaded operators
    
    ok((1 && $v),"$tpre a vertex is boolean-true");
    ok((1 && $v_same),"$tpre a different ref to the same vertex is boolean-true");
    ok((1 && $v_diff),"$tpre a different vertex is boolean-true");
    
    # these two work whether or not they really work, don't they?
    ok(!($v_diff == $v),
        "$tpre vertex == different vertex is false");    
    ok($v_diff != $v,
        "$tpre vertex != different vertex is true");

    # this is the key one to test
    ok($v_same == $v,
        "$tpre vertex == same vertex from different ref is true");

    # now make sure the overload is inherited when we upgrade to
    # auto-scaling class
    
    bless $v, $vertex_class_expected."Scaled";
    ok($v->isa($vertex_class_expected."Scaled"), 
        "$tpre subclassed vertex");
    ok($v_same == $v,
        "$tpre subclassed vertex == same vertex from different ref with original class");
    

    $b->clear();

}

__END__
