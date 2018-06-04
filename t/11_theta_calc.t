#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 22;
use Boost::Polygon::Voronoi qw(calc_Q tangent_angle);

eval {require "t/svg_test_helpers.pl"}; 
if ($@) {eval{require "svg_test_helpers.pl"}}
if ($@) {fail("couldn't find svg_test_helpers.pl");}

{

    # TODO : We have mostly the same enge tangent angle finding code in Perl,
    #        in the module's .pm, for use with the stock Voronoi diagram,
    #        and in C++ for the medial axis. For the medial axis, we figure
    #        those angles and store them with the medial axis edges.
    #        Tests here should make sure all these angles are calculable
    #        for all the cases we target in the code in both perl and c++. 


    my $b = new Boost::Polygon::Voronoi::builder;

    load_test_svg($b);

    #my $ma = get_medial_axis($b);
    my $ma = get_voronoi_diagram($b);

    my @edges = @{$ma->edges()};
    my @cells = @{$ma->cells()};

    make_result_svg(\@edges,\@cells,undef,1.0,'Voronoi diagram');
    
    my $finite_count = scalar(grep $_->is_finite(), @edges);
    ok($finite_count == 18, "18 finite edges? ($finite_count)");

    $finite_count = scalar(grep $_->is_finite() && $_->is_secondary(), @edges);
    ok($finite_count == 8, "8 finite secondary edges? ($finite_count)");

    $finite_count = scalar(grep $_->is_finite() && $_->is_primary(), @edges);
    ok($finite_count == 10, "10 finite primary edges? ($finite_count)");

    $finite_count = scalar(grep $_->is_finite() && $_->is_primary() && $_->is_curved(), @edges);
    ok($finite_count == 4, "4 finite primary curved edges? ($finite_count)");

    

}

__END__
