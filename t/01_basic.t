#!/usr/bin/perl

use strict;
use warnings;



# already cut out most
# now just need to change to B:P:V stuff and build out again

use Test::More tests => 2;
use Boost::Polygon::Voronoi;

{
    my $square = [  # ccw
        [10, 10],
        [20, 10],
        [20, 20],
        [10, 20],
    ];
    my $hole_in_square = [  # cw
        [14, 14],
        [14, 16],
        [16, 16],
        [16, 14],
    ];
    my $polygon = [$square, $hole_in_square];
    my $linestring = [ [5, 15], [30, 15] ];
    my $linestring2 = [ [40, 15], [50, 15] ];  # external
    my $multilinestring = [ [ [5, 15], [30, 15] ], [ [40, 15], [50, 15] ] ];
    
    {
        my $b = new Boost::Polygon::Voronoi::builder;
        $b->insert_point(0,50);
        $b->insert_point(50,50);
        my $ma = $b->medial_axis();
        my @trigger = @{$ma->edges()};
        is 1,1,'fake test 1';

        #is polygon_area([$square, $hole_in_square]), 10*10 - 2*2, 'polygon area';
    }
    {
        is 1,1,'fake test 2';
        #my $intersection =
        #    polygon_multi_linestring_intersection($polygon, [$linestring]);
        #is_deeply $intersection, [
        #    [ [10, 15], [14, 15] ],
        #    [ [16, 15], [20, 15] ],
        #], 'line is clipped to square with hole';
    }
    
}


__END__
