#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 51;
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

# does autoscale shift and scale input geometry as expected
# and unscale and unshift
# and does it wrap/re-bless all topological objects into classes that
# do the scaling and re-bless all their returned topological objects
# all while leaving the underlying geometry in it's scaled-up state

foreach my $vd_or_ma ('Voronoi diagram','medial axis') {
        
    my $bpv = Boost::Polygon::Voronoi->new({auto_scale => 1});
    my $poly=[[0,0],[400,0],[400,500],[0,500],[0,0]]; # simple rectangle
    $bpv->addPolygon($poly);

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
        
    # min, max x is 0, 400, so xoff is minx+xspan/2 = 200
    # min, max y is 0, 500, so yoff is miny+yspan/2 = 250
    # the 32 bit signed integer limit we're using is 2147483647, +/-
    # use the larger yspan to figure the scale
    my $xoff_should_be = -200;
    my $yoff_should_be = -250;
    my $scale_should_be = int(2147483647 / (500/2));

    is_deeply($bpv->{scale},[$xoff_should_be,$yoff_should_be,$scale_should_be],
        "scale values correct");

    my $vert_ref = ref($bpv->vertices()->[0]);
    $vert_ref =~ /(Boost::Polygon::[^:]+::[^:]+::)(.*)$/;
    my $base_scale_class = $1;
    my $scale_symbol = $2;
    #diag("vert ref: ".ref($bpv->vertices()->[0]));

    # this duplicates what's done in scale setup code
    my $scale_symbol_string_should_be = join('::', map { 
        tr/\.\-\+0123456789eE/DNPabcdefghijEE/r # note e and E map to E
    } @{$bpv->{scale}});

    is($scale_symbol, $scale_symbol_string_should_be,
        "encoded scale values in package name correspond to scale values");

    # this decodes scale values from the end of the class name
    my ($rec_xoff, $rec_yoff, $rec_scale) = map { 
        tr/DNPabcdefghijE/\.\-\+0123456789e/r  # note E just maps to e
    } split('::', $scale_symbol);
    
    is_deeply([$rec_xoff, $rec_yoff, $rec_scale], $bpv->{scale},
        "scale values match values decoded from package name");

    my $package_base = "Boost::Polygon::"
                     .($vd_or_ma eq 'Voronoi diagram' ? 'Voronoi':'MedialAxis');
    my $vert_package = $package_base
                     .'::VertexScaled::'
                     .$scale_symbol_string_should_be;
    my $edge_package = $package_base
                     .'::EdgeScaled::'
                     .$scale_symbol_string_should_be;
    my $cell_package = $package_base
                     .'::CellScaled::'
                     .$scale_symbol_string_should_be;

    ok(scalar(grep {ref($_) eq $vert_package} @{$bpv->vertices()})
       == scalar(@{$bpv->vertices()})
       , "all vertices blessed into correct package"
       
    );
    ok(scalar(grep {ref($_) eq $edge_package} @{$bpv->edges()})
       == scalar(@{$bpv->edges()})
       , "all edges blessed into correct package"
       
    );
    ok(scalar(grep {ref($_) eq $cell_package} @{$bpv->cells()})
       == scalar(@{$bpv->cells()})
       , "all cells blessed into correct package"
    );


    my %orig_points = (map {$_->[0].'_'.$_->[1] => 0} @$poly);
    my @orig_points_keys = keys %orig_points;
    my %orig_points_scaled = (
        # note that we int() these scaled up coords
        # because we're going to compare them to 
        # coords that had that happen on their way into C++
        # even though they come out as floats
        # (output float vertices that are copies of input integer points)
        map {      int(($_->[0] + $xoff_should_be) * $scale_should_be)
              .'_'.int(($_->[1] + $yoff_should_be) * $scale_should_be)
             => 0
            } @$poly);
    my @orig_points_scaled_keys = keys %orig_points_scaled;

    # The sprintf() business here is because the scaling with floats
    # has lost precision in the round trip.
    # Need to look into that.
    # Do we loose that precision when scaled up floats get truncated?
    # Are there scale values that will limit or avoid that?
    # Can we at least get the bounds of error caused by that, and note
    # that in the docs.
    # Kinda stinks that you output vertices might not match up exactly with
    # your input points. But maybe that's the price for using the auto
    # scale feature. And maybe it doesn't matter for most applications.
    # If it does matter, users can turn autoscaling off and take care to 
    # set up their own coords, as scaled up integers maybe, and those
    # should be okay round trip.
    # 
    $orig_points{sprintf("%.0f",$_->x()).'_'.sprintf("%.0f",$_->y())}++ for @{$bpv->vertices()};
    $orig_points_scaled{$_->x().'_'.$_->y()}++ for @{$bpv->{vd}->vertices()};

    is(scalar((grep {$_ == 0} map $orig_points{$_}, @orig_points_keys)),
       0,
       "all input points found in autoscaled output vertices");

    is(scalar((grep {$_ == 0} map $orig_points_scaled{$_}, @orig_points_scaled_keys)),
       0,
       "all input points found in underlying output vertices, in scaled up form");
    
    $bpv->clear();

}

__END__
