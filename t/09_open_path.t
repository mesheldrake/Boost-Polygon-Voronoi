#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 71;
use Boost::Polygon::Voronoi qw(calc_Q tangent_angle);

eval {
    require "t/svg_test_helpers.pl";
    require "t/general_topology_checks.pl";
    }; 
if ($@) {eval{
    require "svg_test_helpers.pl";
    require "general_topology_checks.pl";
};}
ok(! $@ , "loaded general tests");


foreach my $vd_or_ma ('Voronoi diagram','medial axis') {
    
    my $isma = $vd_or_ma=~/medial/;

    my $b = new Boost::Polygon::Voronoi::builder;

    load_test_svg($b);

    my $vd = $vd_or_ma eq 'Voronoi diagram'
           ? get_voronoi_diagram($b) 
           : get_medial_axis($b);

    general_topology_checks(
        $vd,
        {
        edges              => 138 + ($isma ? 20 : 0),
        finite             =>  96 + ($isma ? 20 : 0),
        infinite           =>  42,
        finite_secondary   =>  28 + ($isma ? 20 : 0),
        infinite_secondary =>  32,
        finite_primary     =>  68,
        infinite_primary   =>  10,
        curved             =>   8,
        vertices           =>  35 + ($isma ? 10 : 0),
        cells              =>  35,
        degenerate_cells   =>   0,
        _cell_loop_limit   => 1000
        });

    make_result_svg($vd->edges(),$vd->cells(),$vd->vertices(),1.5,$vd_or_ma);

    $b->clear();

}

__END__
