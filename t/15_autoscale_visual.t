#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 31;
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
        
    # using Voronoi builder 
    my $b = new Boost::Polygon::Voronoi::builder;

    # using wrapped api
    my $bpv = Boost::Polygon::Voronoi->new({auto_scale => 0});
    load_test_svg($bpv->{builder});

    if ($vd_or_ma eq 'Voronoi diagram' ) {
        $bpv->voronoi_diagram();
    }
    else {
        $bpv->medial_axis();
    }

    # make sure visual test output works with wrapped api's cached results
    make_result_svg($bpv->edges(),$bpv->cells(),$bpv->vertices(),1.5,$vd_or_ma);

    $bpv->clear();

}

__END__
