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
    
    my $b = new Boost::Polygon::Voronoi::builder;

    load_test_svg($b);

    my $vd = $vd_or_ma eq 'Voronoi diagram'
           ? get_voronoi_diagram($b) 
           : get_medial_axis($b);
           
    my $isma = $vd_or_ma eq 'Voronoi diagram' ? 0:1;

    general_topology_checks(
        $vd,
        {
        edges              => 74 + ($isma ? 2 : 0),
        finite             => 48 + ($isma ? 2 : 0),
        infinite           => 26,
        finite_secondary   => 12,
        infinite_secondary => 24,
        finite_primary     => 36 + ($isma ? 2 : 0),
        infinite_primary   => 2,
        curved             => 12,
        vertices           => 20 + ($isma ? 1 : 0),
        cells              => 18,
        degenerate_cells   => 0,
        _cell_loop_limit   => 1000
        });

    make_result_svg($vd->edges(),$vd->cells(),$vd->vertices(),1.5,$vd_or_ma);
    
    $b->clear();

}

__END__
