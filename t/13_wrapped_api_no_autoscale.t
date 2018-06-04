#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 113;
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
    
    # compare the lower level interface results
    # to results from our wrapped api
    
    # using Voronoi builder 
    my $b = new Boost::Polygon::Voronoi::builder;
    load_test_svg($b);

    # using wrapped api
    my $bpv = Boost::Polygon::Voronoi->new({auto_scale => 0});
    # make this match what's loaded above from svg
    my $poly=[[0,0],[400,0],[400,500],[0,500],[0,0]];
    # match scale too
    for (@$poly) {$_->[0]*=10000;$_->[1]*=10000;}
    $bpv->addPolygon($poly);


    # builder
    my $vd = $vd_or_ma eq 'Voronoi diagram'
           ? get_voronoi_diagram($b) 
           : get_medial_axis($b);

    general_topology_checks(
        $vd,
        {
        edges              => 26,
        finite             => 10,
        infinite           => 16,
        finite_secondary   => 0,
        infinite_secondary => 16,
        finite_primary     => 10,
        infinite_primary   => 0,
        curved             => 0,
        vertices           => 6,
        cells              => 8,
        degenerate_cells   => 0,
        _cell_loop_limit   => 1000
        });

    # wrapped
    if ($vd_or_ma eq 'Voronoi diagram' ) {
        $bpv->voronoi_diagram();
    }
    else {
        $bpv->medial_axis();
    }
    general_topology_checks(
        $bpv,
        {
        edges              => 26,
        finite             => 10,
        infinite           => 16,
        finite_secondary   => 0,
        infinite_secondary => 16,
        finite_primary     => 10,
        infinite_primary   => 0,
        curved             => 0,
        vertices           => 6,
        cells              => 8,
        degenerate_cells   => 0,
        _cell_loop_limit   => 1000
        });
        

    my @builder_coords = map {[$_->x(),$_->y()]} @{$vd->vertices()};
    my @wrapped_coords = map {[$_->x(),$_->y()]} @{$bpv->vertices()};

    is_deeply(\@wrapped_coords, \@builder_coords,
        "same result vertex coordinates from lower & higher level interfaces");

    ok(scalar(grep {ref($_) eq "Boost::Polygon::"
                     .($vd_or_ma eq 'Voronoi diagram' ? 'Voronoi':'MedialAxis')
                     ."::Vertex"
       } @{$bpv->vertices()})
       ==
       scalar(@{$bpv->vertices()})
       , "all vertices blessed into correct package"
       
    );
    ok(scalar(grep {ref($_) eq "Boost::Polygon::"
                     .($vd_or_ma eq 'Voronoi diagram' ? 'Voronoi':'MedialAxis')
                     ."::Edge"
       } @{$bpv->edges()})
       ==
       scalar(@{$bpv->edges()})
       , "all edges blessed into correct package"
       
    );
    ok(scalar(grep {ref($_) eq "Boost::Polygon::"
                     .($vd_or_ma eq 'Voronoi diagram' ? 'Voronoi':'MedialAxis')
                     ."::Cell"
       } @{$bpv->cells()})
       ==
       scalar(@{$bpv->cells()})
       , "all cells blessed into correct package"
       
    );
    
    # make sure visual test output works with these
    # new results
    make_result_svg($bpv->edges(),$bpv->cells(),$bpv->vertices(),1.5,$vd_or_ma);

    $b->clear();
    $bpv->clear();

}

__END__
