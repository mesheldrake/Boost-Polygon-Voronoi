package Boost::Polygon::Voronoi;
{
{
  $Boost::Polygon::Voronoi::VERSION = '0.01';
  $Boost::Polygon::Voronoi::AUTO_VALIDATE_SEGMENTS = 1;
  $Boost::Polygon::Voronoi::AUTO_SCALE = 1;
}
# ABSTRACT: Bindings for the Boost Polygon Voronoi library
use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);

use XSLoader;
XSLoader::load('Boost::Polygon::Voronoi', $Boost::Polygon::Voronoi::VERSION);

our @EXPORT_OK = qw(addPolygon addPolygons addPolyline addPolylines
                    addSegment addSegments addPoint addPoints
                    validateSegments calc_Q tangent_angle get_curve_focus);

my $int_limit = 2147483647; # signed int32 is -2147483648 to 2147483647

my $pi = 4 * atan2(1,1);
sub tan {return sin($_[0])/cos($_[0]);}

sub new {
    my $class = shift;
    my $opts = @_ && ref($_[0])=~/HASH/i ? shift : {};
    my $self = {};
    bless $self, $class;
    
    $self->{builder}   = new Boost::Polygon::Voronoi::builder();
    
    $self->clear();

    $self->{auto_scale} = exists $opts->{auto_scale}
                           ? ($opts->{auto_scale} ? 1:0)
                           : $Boost::Polygon::Voronoi::AUTO_SCALE;
    $self->{auto_validate} = exists $opts->{auto_validate}
                           ? ($opts->{auto_validate} ? 1:0)
                           : $Boost::Polygon::Voronoi::AUTO_VALIDATE_SEGMENTS;
    
    return $self;
}

sub clear {
    my $self = shift;
    
    #input staging arrays
    $self->{polylines} = []; # both closed polygons and open polylines
    $self->{points}    = [];
    $self->{segments}  = [];

    # output cache, so you only have trigger the bulk C to Perl conversion once
    # for $self->{vd}->cells(), ->vertices(), ->edges();
    $self->{cells}     = [];
    $self->{edges}     = [];
    $self->{vertices}  = [];

    # Voronoi diagram result from $self->{builder} will go here
    $self->{vd} = undef;

    # for re-blessing things when autoscaling is in effect
    $self->{cell_class} = undef;
    $self->{edge_class} = undef;
    $self->{vertex_class} = undef;
    $self->{foot_class} = undef; # for medial axis only

    # clears out objects "living in C++"
    $self->{builder}->clear();
}

### adding geometry

sub addPolygon { # [ [x,y], ... ]
    my ($self, $poly) = @_;
    if (scalar(@{$poly}) > 0) {
        if    (@{$poly} == 1 ) {$self->addPoint([$poly->[0]->[0],$poly->[0]->[1]]);}
        elsif (@{$poly} == 2 ) {$self->addSegment([$poly->[0]->[0],$poly->[0]->[1],$poly->[1]->[0],$poly->[1]->[1]]);}
        else { push @{$self->{polylines}}, ['polygon',$poly]; }
    }
}

sub addPolygons { # [ [ [x,y], ... ] , ... ]
    my ($self, $polygons) = @_;
    warn "addPolygons requires a reference to an array of polygons" if (!$polygons || !ref($polygons));
    $self->addPolygon($_) for @$polygons;
}

sub addPolyline { # [ [x,y], ... ]
    my ($self, $poly) = @_;
    if (scalar(@{$poly}) > 0) {
        if    (@{$poly} == 1 ) {$self->addPoint([$poly->[0]->[0],$poly->[0]->[1]]);}
        elsif (@{$poly} == 2 ) {$self->addSegment([$poly->[0]->[0],$poly->[0]->[1],$poly->[1]->[0],$poly->[1]->[1]]);}
        push @{$self->{polylines}}, ['polyline',$poly];
    }
}

sub addPolylines { # [ [ [x,y], ... ] , ... ]
    my ($self, $polylines) = @_;
    warn "addPolylines requires a reference to an array of polylines" if (!$polylines || !ref($polylines));
    $self->addPolyline($_) for @$polylines;
}

sub addSegment { # [x1, y1, x2, y2]
    my ($self, $segment) = @_;
    if (scalar(@$segment) == 2 && ref($segment->[0]) && ref($segment->[1])) {
        push @{$self->{segments}}, [$segment->[0]->[0],$segment->[0]->[1],
                                    $segment->[1]->[0],$segment->[1]->[1]];
    }
    elsif (scalar(@$segment) == 4) {
        push @{$self->{segments}}, $segment;
    }
    else {
        warn "addSegment requires array ref of four coordintates, "
            ,"e.g. [x1,y1,x2,y2], or of two points, e.g. [[x1,y1],[x2,y2]]";
    }
}

sub addSegments { # [ [x1, y1, x2, y2], ... ]
    my ($self, $segments) = @_;
    warn "addSegments requires a reference to an array of segments" if (!$segments || !ref($segments));
    push @{$self->{segments}}, @$segments;
}

sub addPoint { # [x,y] point unaffiliated with PLSG segments
    my ($self, $point) = @_;
    # points might be more than 2 dimensional, so just take the first two coords
    push @{$self->{points}}, [$point->[0],$point->[1]];
}

sub addPoints { # [ [x,y], ... ] points unaffiliated with PLSG segments
    my ($self, $points) = @_;
    # points might be more than 2 dimensional, so just take the first two coords
    push @{$self->{points}}, map [$_->[0],$_->[1]], @$points;
    }

### preparing geometry for Voronoi builder

sub prep_poly_segs {
    my $self = shift;
    foreach my $pl (@{$self->{polylines}}) {
        my $is_gon = $pl->[0] eq 'polygon' ? 1 : 0;
        my $poly = $pl->[1];
        if ($is_gon && (   $poly->[-1]->[0] ne $poly->[0]->[0] 
                        || $poly->[-1]->[1] ne $poly->[0]->[1])) {
        push @$poly, [$poly->[0]->[0],$poly->[0]->[1]];
        }
        my @segset;
        for (my $i = 0; $i < $#$poly; $i++) {
            push @segset, [$poly->[ $i ]->[0], $poly->[ $i ]->[1],
                           $poly->[$i+1]->[0], $poly->[$i+1]->[1]];
        }
        if (scalar(@segset) > 0) { push @{$self->{segments}}, @segset; }
    }
}

sub insert_segments {
    my $self = shift;
    my $segs = $self->{segments};
    my $insert_last_index = -1;
    my $bad_seg_report = '';
    for (my $i=0;$i<@{$segs};$i++) {
        eval {
            my ($x1,$y1,$x2,$y2) = @{$segs->[$i]};
            if (defined $self->{scale}) {
                for ($x1,$x2) { $_ += $self->{scale}->[0]; 
                                $_ *= $self->{scale}->[2]; }
                for ($y1,$y2) { $_ += $self->{scale}->[1]; 
                                $_ *= $self->{scale}->[2]; }
            }
            $insert_last_index = $self->{builder}->insert_segment($x1,$y1,$x2,$y2);
        };
        if ($@) {
            $bad_seg_report = 'trouble inserting this segment: ['.
                              join(',',@{$segs->[$i]}).
                              ']'.$@;
            last;
        }
    }
    if ($bad_seg_report) { warn "insert_segment error: $bad_seg_report"; }
    if (defined wantarray) {return $insert_last_index;} 
}

sub insert_points {
    my $self = shift;
    my $points = $self->{points};
    my $insert_last_index = -1;
    my $bad_point_report = '';
    for (my $i=0;$i<@{$points};$i++) {
        eval {
            my ($x,$y) = @{$points->[$i]};
            if (defined $self->{scale}) {
                $x += $self->{scale}->[0]; 
                $x *= $self->{scale}->[2];
                $y += $self->{scale}->[1]; 
                $y *= $self->{scale}->[2];
            }
            $insert_last_index = $self->{builder}->insert_point($x,$y);
        };
        if ($@) {
            $bad_point_report = 'trouble inserting this point: ['.
                              join(',', @{$points->[$i]}).
                              ']'.$@;
            last;
        }
    }
    if ($bad_point_report) { warn "insert_point error: $bad_point_report"; }
    if (defined wantarray) {return $insert_last_index;} 
}

sub configure_scale {
    my ($self, $isma) = @_;
    my $polycnt   = scalar(@{$self->{polylines}});
    my $segcnt    = scalar(@{$self->{segments}});
    my $pointcnt  = scalar(@{$self->{points}});
    return [0,0,1] if (   $polycnt  == 0
                       && $segcnt   == 0
                       && $pointcnt == 0
                       );
    my $a_point = $polycnt 
                ? $self->{polylines}->[0]->[1]->[0]
                : ($segcnt
                  ? [$self->{segments}->[0]->[0],$self->{segments}->[0]->[1]]
                  : $self->{points}->[0]);
    my ($minx, $maxx, $miny, $maxy) = ($a_point->[0],$a_point->[0],
                                       $a_point->[1],$a_point->[1]);
    foreach my $p ((map @{$_->[1]}, @{$self->{polylines}}),
                   map {([$_->[0],$_->[1]],[$_->[2],$_->[3]])} @{$self->{segments}}, 
                   @{$self->{points}}
                 ) {
        $minx = $p->[0] if $minx > $p->[0];
        $maxx = $p->[0] if $maxx < $p->[0];
        $miny = $p->[1] if $miny > $p->[1];
        $maxy = $p->[1] if $maxy < $p->[1];
    }
    
    my $xspan = $maxx - $minx;
    my $xoff = $minx + $xspan/2;
    my $yspan = $maxy - $miny;
    my $yoff = $miny + $yspan/2;
    my $x_scale_up = int($int_limit / ($xspan/2));
    my $y_scale_up = int($int_limit / ($yspan/2));
    my $scale_up = ($x_scale_up < $y_scale_up) ? $x_scale_up : $y_scale_up;
    
    $self->{scale} = [-$xoff, -$yoff, $scale_up];
    
    # Why this dynamic package/class/symbol table mucking nonsense?
    # So we can still work with the underlying C++ objects
    # yet add auto-scaling behavior, simply by blessing the
    # refs we get to those objects.
    # This means we don't have to copy over or recreate
    # all the topology of the edges, vertices and cells, which are
    # all cross referenced to eachother with pointers in the C++ objects.
    # For vertices, this means we don't have to
    # store our three scaling values with each vertex object.
    # The intent is to preserve and use the data "living in C++" as
    # much as possible, to take advantage of the reduced memory footprint.
    # Also, we only want to do this extra blessing when needed for autoscaling.
    # We still have the option of turning autoscaling off and getting
    # result objects "closer" to their underlying C++ objects, which should be
    # "faster".

    # create unique package names based on the specific scale config
    my $diagram_class = 'Boost::Polygon::' . ($isma ? 'MedialAxis' : 'Voronoi');
    my $scale_symbol_string = join('::', map { 
        tr/\.\-\+0123456789eE/DNPabcdefghijEE/r
        } @{$self->{scale}});

    my $vertex_class = $diagram_class.'::VertexScaled::' . $scale_symbol_string;
    my $edge_class   = $diagram_class.'::EdgeScaled::'   . $scale_symbol_string;
    my $cell_class   = $diagram_class.'::CellScaled::'   . $scale_symbol_string;
    my $foot_class   = $diagram_class.'::FootScaled::'   . $scale_symbol_string if $isma;

    # remember these class names for when it comes time to import and 
    # re-bless results from C++ objects
    $self->{vertex_class} = $vertex_class;
    $self->{edge_class} = $edge_class;
    $self->{cell_class} = $cell_class;
    $self->{foot_class} = $foot_class                                  if $isma;

    no strict 'refs';

    if (!%{"${vertex_class}::"}) { # avoid redefining if package exists
            
        # create the packages by modifying the symbol table

        # these scaled subclasses need to override
        # all the methods that return coordinates
        # to apply the scale transforms
        *{"${vertex_class}::x"} = sub {$_[0]->super_x() / $self->{scale}->[2] - $self->{scale}->[0]};
        *{"${vertex_class}::y"} = sub {$_[0]->super_y() / $self->{scale}->[2] - $self->{scale}->[1]};
        if ($isma) {
            *{"${vertex_class}::r"} = sub {$_[0]->super_r() / $self->{scale}->[2]};
            *{"${foot_class}::x"} = sub {$_[0]->super_x() / $self->{scale}->[2] - $self->{scale}->[0]};
            *{"${foot_class}::y"} = sub {$_[0]->super_y() / $self->{scale}->[2] - $self->{scale}->[1]};
        }
        
        # inherit from the base classes that provide the super_*() methods
        # just used above and that ensure methods that return objects from 
        # underlying topology will bless those objects into the appropriate 
        # auto-scaling class
        push @{"${vertex_class}::ISA"}, $diagram_class.'::VertexScaled';
        push @{"${edge_class}::ISA"}  , $diagram_class.'::EdgeScaled';
        push @{"${cell_class}::ISA"}  , $diagram_class.'::CellScaled';
        push @{"${foot_class}::ISA"}  , $diagram_class.'::FootScaled'      if $isma;

    }
}

sub validateSegments {
    my $segments = shift;
    my %dups;
    return 0 if scalar(@$segments) == 0;
    for (my $i=0;$i<@$segments;$i++) {
        my $k = join('^',@{$segments->[$i]});
        if ( ! exists $dups{$k} ) {
            $dups{$k} = 1;
        }
        else {
            return 0;
        }
    }
    return 1;
}

sub cleanSegments {
    my $self = shift;
    my %dups;
    my @newsegs;
    for (my $i=0;$i<@{$self->{segments}};$i++) {
        my $k = join('^',@{$self->{segments}->[$i]});
        if ( ! exists $dups{$k} ) {
            $dups{$k} = 1;
            push @newsegs, $self->{segments}->[$i];
        }
        else {
            #return 0;
        }
    }
    if (scalar(@{$self->{segments}}) != scalar(@newsegs)) {
        #warn "removed ",
        #     (scalar(@{$self->{segments}}) - scalar(@newsegs)),
        #     " duplicate segments";
        $self->{segments} = \@newsegs;
    }
}


### execute 

sub medial_axis {
    my $self = shift;
    if ($self->{auto_scale}) {
        $self->configure_scale('ma');
    } else { delete $self->{scale}; }
    $self->prep_poly_segs();
    if (!validateSegments($self->{segments})) {
        #warn "cleaning out duplicate segments\n";
        $self->cleanSegments();
    }
    $self->insert_segments();
    $self->insert_points();
    $self->{vd} = $self->{builder}->medial_axis();
}
sub voronoi_diagram {
    my $self = shift;
    if ($self->{auto_scale}) {
        $self->configure_scale();
    } else { delete $self->{scale}; }
    $self->prep_poly_segs();
    $self->insert_segments();
    $self->insert_points();
    $self->{vd} = $self->{builder}->voronoi_diagram();
}

### access results

sub cells {
    my $self = shift;
    if (scalar(@{$self->{cells}}) == 0) {
        @{$self->{cells}} = @{$self->{vd}->cells()}; # cache
        if ($self->{cell_class}) { # upgrade if scaling in effect
            bless($_, $self->{cell_class}) for @{$self->{cells}};
        }
    }
    return $self->{cells};
}
sub edges {
    my $self = shift;
    if (scalar(@{$self->{edges}}) == 0) {
        @{$self->{edges}} = @{$self->{vd}->edges()}; # cache
        if ($self->{edge_class}) { # upgrade if scaling in effect
            bless($_, $self->{edge_class}) for @{$self->{edges}};
        }
    }
    return $self->{edges};
}
sub vertices {
    my $self = shift;
    if (scalar(@{$self->{vertices}}) == 0) {
        @{$self->{vertices}} = @{$self->{vd}->vertices()}; # cache
        if ($self->{vertex_class}) { # upgrade if scaling in effect
            bless($_, $self->{vertex_class}) for @{$self->{vertices}};
        }
    }
    return $self->{vertices};
}

# This is part of what tangent_angle() does.
# Keep this in sync with that, or have that call this.
# Currently this is used in one place in
# the visual helper code for tests.
# It is a fundamental component of recovering parabola
# data from the stock Boost::Polygon::Voronoi Voronoi diagram.
sub get_curve_focus {
    my $edge = shift;
    my $focus;
    if ($edge->is_curved()) {
        my $twin_has_point = $edge->cell()->contains_point() ? 0 : 1;
        my $edge_with_point = $twin_has_point
                     ? $edge->twin()
                     : $edge;

        # A curved primary segment should always have secondaries
        # before and after it, coming from the parabola focus and
        # going to it. Around vertices at the ends of infinite segments,
        # one of these secondaries may be missing, so use the other.

        my $secondary_with_point = $edge_with_point;

        my $focus;
        
        if (   $edge_with_point->prev_vd()->is_secondary() && $edge_with_point->prev_vd()->is_finite()
            && $edge_with_point->prev_vd()->prev_vd()->is_secondary() && $edge_with_point->prev_vd()->prev_vd()->is_finite()
            &&    ${$edge_with_point->prev_vd()->vertex0()}
               == ${$edge_with_point->prev_vd()->prev_vd()->vertex1()}
           ) {
            $focus = $edge_with_point->prev_vd()->vertex0();
        }
        elsif ($edge_with_point->next_vd()->is_secondary() && $edge_with_point->next_vd()->is_finite()
            && $edge_with_point->next_vd()->next_vd()->is_secondary() && $edge_with_point->next_vd()->next_vd()->is_finite()
            &&    ${$edge_with_point->next_vd()->vertex1()}
               == ${$edge_with_point->next_vd()->next_vd()->vertex0()}
           ) {
           $focus = $edge_with_point->next_vd()->vertex1();
        }
        else {
            my $l=300;
            do {
                $secondary_with_point = $secondary_with_point->prev_vd();
            } while ($secondary_with_point->is_primary() && $$secondary_with_point != $$edge_with_point && $l-- > 0);
            die "safety while limit reached" if ! ($l > 0);
            $focus = $secondary_with_point->vertex0();

            if (!$focus) {
                $l=300;
                $secondary_with_point = $edge_with_point;
                do {
                    $secondary_with_point = $secondary_with_point->next_vd();
                } while ($secondary_with_point->is_primary() && $$secondary_with_point != $$edge_with_point && $l-- > 0);
                die "safety while limit reached" if ! ($l > 0);
                $focus = $secondary_with_point->vertex1();
            }
        }
    }
    return $focus;
}
# Given an Voronoi edge, determine the tangent angle at the start of that edge
# heading out from edge->vertex0(). The main task here is either to calculate
# the tangent angle for a curved edge (parabola math), or infer the angle
# for an infinite edge from finite edges nearby.

# catch infinite recursion
$Boost::Polygon::Voronoi::last_tangent_angle_case = -1;

sub tangent_angle {
    my $edge = shift;
    
    my $theta;

    if ($edge->is_finite()) {

        if ($edge->is_curved()) {

            my $twin_has_point = $edge->cell()->contains_point() ? 0 : 1;
            my $edge_with_point = $twin_has_point
                         ? $edge->twin()
                         : $edge;

            # A curved primary segment should always have secondaries
            # before and after it, coming from the parabola focus and
            # going to it. Around vertices at the ends of infinite segments,
            # one of these secondaries may be missing, so use the other.

            my $secondary_with_point = $edge_with_point;

            my $focus;
            
            if (   $edge_with_point->prev_vd()->is_secondary() && $edge_with_point->prev_vd()->is_finite()
                && $edge_with_point->prev_vd()->prev_vd()->is_secondary() && $edge_with_point->prev_vd()->prev_vd()->is_finite()
                &&    ${$edge_with_point->prev_vd()->vertex0()}
                   == ${$edge_with_point->prev_vd()->prev_vd()->vertex1()}
               ) {
                $focus = $edge_with_point->prev_vd()->vertex0();
            }
            elsif ($edge_with_point->next_vd()->is_secondary() && $edge_with_point->next_vd()->is_finite()
                && $edge_with_point->next_vd()->next_vd()->is_secondary() && $edge_with_point->next_vd()->next_vd()->is_finite()
                &&    ${$edge_with_point->next_vd()->vertex1()}
                   == ${$edge_with_point->next_vd()->next_vd()->vertex0()}
               ) {
               $focus = $edge_with_point->next_vd()->vertex1();
            }
            else {
                my $l=300;
                do {
                    $secondary_with_point = $secondary_with_point->prev_vd();
                } while ($secondary_with_point->is_primary() && $$secondary_with_point != $$edge_with_point && $l-- > 0);
                die "safety while limit reached" if ! ($l > 0);
                $focus = $secondary_with_point->vertex0();

                if (!$focus) {
                    $l=300;
                    $secondary_with_point = $edge_with_point;
                    do {
                        $secondary_with_point = $secondary_with_point->next_vd();
                    } while ($secondary_with_point->is_primary() && $$secondary_with_point != $$edge_with_point && $l-- > 0);
                    die "safety while limit reached" if ! ($l > 0);
                    $focus = $secondary_with_point->vertex1();
                }
            }

            my $p0 = $edge_with_point->vertex0();
            my $p2 = $edge_with_point->vertex1();

            #warn 'foc: ',$focus->x(),",",$focus->y(),"\n";
            #warn 'p0 : ',$p0->x(),",",$p0->y(),"\n";
            #warn 'p2 : ',$p2->x(),",",$p2->y(),"\n";
            
            # directrix angle = alpha - gamma
            # alpha = angle of line p0 to p2
            # gamma = angle opposite delta_r in a right triangle with b as base, r as height, p0 to p2 line as hypotenuse
            # base b is derived from known delta_r and hypotenuse
            my $delta_x = $p2->x() - $p0->x();
            my $delta_y = $p2->y() - $p0->y();
            my $alpha = atan2($delta_y,$delta_x);
            
            my $delta_r = (  sqrt(($p2->x() - $focus->x())**2 + ($p2->y() - $focus->y())**2)
                           - sqrt(($p0->x() - $focus->x())**2 + ($p0->y() - $focus->y())**2));
            my $h = sqrt($delta_x**2 + $delta_y**2);
            my $b = sqrt($h**2 - $delta_r**2);
            my $gamma = atan2($delta_r, $b);



            my $directrix_angle = $alpha - $gamma;

            while ($directrix_angle >  $pi) {$directrix_angle -= ($pi * 2);}
            while ($directrix_angle < -$pi) {$directrix_angle += ($pi * 2);}
            
            # Now ready to figure the tangent angle at requested edge's vertex0
            # although, if we had to use the twin to figure the above,
            # we're really figuring the angle at the twin's vertex1, which is
            # 180 degrees opposite of what we want.
            # So we'll turn it around when we're done, if that's that case.
            
            my $angle_focus_p = atan2($edge->vertex0()->y() - $focus->y(), $edge->vertex0()->x() - $focus->x());

            my $angle_focus_p_un_dir_rot = $angle_focus_p - $directrix_angle;
            while ($angle_focus_p_un_dir_rot >  $pi) {$angle_focus_p_un_dir_rot -= ($pi * 2);}
            while ($angle_focus_p_un_dir_rot < -$pi) {$angle_focus_p_un_dir_rot += ($pi * 2);}

            # Fixes 180 degree wrongness in 2nd quadrant.
            # Consider that this signed angle $angle_focus_p_un_dir_rot has value -90 degrees
            # when the point is at the bottom of the parabola. We then have a 180 degree sweep
            # both positive and negative from there, for the parabola's two legs (of course
            # never reaching those limits). So on the right side that's -90 to 90.
            # On the left that's -90 to -270. But atan2() isn't going to give us the negative
            # angles we want in the 2nd quadrant. So convert all positive 2nd quadrant angles to
            # to negatives by subtracting 2*PI.
            if ($angle_focus_p_un_dir_rot > $pi/2) {$angle_focus_p_un_dir_rot-=2*$pi;}

            # sketched this a couple ways for + and - angles and seems right
            my $theta_un_dir_rot = $pi/4 + $angle_focus_p_un_dir_rot/2;
            
            $theta = $theta_un_dir_rot + $directrix_angle;

            # angle reduce
            while ($theta >  $pi) {$theta -= ($pi * 2);}
            while ($theta < -$pi) {$theta += ($pi * 2);}

            # The 180 degree turn-around, if needed.
            if ($twin_has_point) {
                $theta += $pi;
            }
        }
        else {
            my $p0 = $edge->vertex0();
            my $p2 = $edge->twin()->vertex0();
            $theta = atan2($p2->y() - $p0->y(),
                           $p2->x() - $p0->x());
        }
    }
    else { # infinite

        my $inf_theta_debug = 0;

        if ($edge->is_primary()) {
            # infinite primary linear
            # always has a primary behind that we can get the angle from
            # 3 cases
            
            my $outgoing = $edge->vertex0() ? $edge : $edge->twin();
            
            # four primaries meet
            if (   $outgoing->rot_next_vd()->is_primary()
                && $outgoing->rot_prev_vd()->is_primary()
                && $outgoing->rot_next_vd()->rot_next_vd()->is_primary()
                &&    ${$outgoing->rot_next_vd()->rot_next_vd()}
                   != ${$outgoing}
                &&    ${$outgoing->rot_next_vd()->rot_next_vd()}
                   != ${$outgoing->rot_next_vd()}
                &&    ${$outgoing->rot_next_vd()->rot_next_vd()}
                   != ${$outgoing->rot_prev_vd()}
               ) {
                print "CASE P0" if $inf_theta_debug;
                if ($outgoing->rot_next_vd()->rot_next_vd()->is_finite()) {
                    print ".1 " if $inf_theta_debug;
                    $theta = tangent_angle($outgoing->rot_next_vd()->rot_next_vd()) + $pi;
                }
                elsif (   $outgoing->rot_next_vd()->is_finite()
                       && $outgoing->rot_prev_vd()->is_finite()
                      ) {
                    print ".2 " if $inf_theta_debug;
                    my $a1 = tangent_angle($outgoing->rot_next_vd());
                    my $a2 = tangent_angle($outgoing->rot_prev_vd());
                    # abs() here because we're figuring the counterclockwise
                    # sweep angle, which we add half of to prev() angle to 
                    # get the a1-a2 bisector in the right direction.
                    $theta = $a2 + (abs($a1 - $a2) / 2);
                }
                else {die "unhandled infinite primary tangent angle finding case";}
            }

          #// three primaries meet at a point
            elsif ($outgoing->rot_next_vd()->is_primary()
                && $outgoing->rot_prev_vd()->is_primary()
                &&    ${$outgoing->rot_next_vd()->rot_next_vd()->rot_next_vd()}
                   == ${$outgoing}
                 ) {
              if ($inf_theta_debug) {print "CASE P0.5";}
              if (   $outgoing->rot_next_vd()->is_finite()
                  && $outgoing->rot_prev_vd()->is_finite()
                 ) {
                if ($inf_theta_debug) {print ".1 ";}
                my $a1 = tangent_angle($outgoing->rot_prev_vd());
                my $a2 = tangent_angle($outgoing->rot_next_vd());
                $theta = $a1 + ((2*$pi - abs($a2 - $a1)) / 2);
            }
            else {
                die("unhandled infinite primary tangent angle finding case\n");
                }
            }
            elsif ($outgoing->rot_next_vd()->is_primary()) {
                print "CASE P1" if $inf_theta_debug;
                if ($outgoing->rot_prev_vd()->rot_prev_vd()->is_primary()) {
                    if (   ${$outgoing->rot_prev_vd()->rot_prev_vd()}
                        == ${$outgoing->rot_next_vd()}) {
                        print ".1 " if $inf_theta_debug;
                        $theta = tangent_angle($outgoing->rot_next_vd()) + $pi;
                    }
                    else {
                        print ".2 " if $inf_theta_debug;
                        my $a1 = tangent_angle($outgoing->rot_next_vd());
                        my $a2 = tangent_angle($outgoing->rot_prev_vd()->rot_prev_vd());
                        $theta = $a2 + (abs($a1 - $a2) / 2);
                    }
                }
                else {die "unhandled infinite primary tangent angle finding case";}
            }
            elsif ($outgoing->rot_prev_vd()->is_primary()) {
                print "CASE P2" if $inf_theta_debug;
                if ($outgoing->rot_next_vd()->rot_next_vd()->is_primary()) {
                    if (   ${$outgoing->rot_next_vd()->rot_next_vd()}
                        == ${$outgoing->rot_prev_vd()}) {
                        print ".1 " if $inf_theta_debug;
                        $theta = tangent_angle($outgoing->rot_prev_vd());
                        }
                    else {
                        print ".2 " if $inf_theta_debug;
                        my $a1 = tangent_angle($outgoing->rot_prev_vd());
                        my $a2 = tangent_angle($outgoing->rot_next_vd()->rot_next_vd());
                        $theta = ($a1 + $a2) / 2;
                    }

                }
                else {die "unhandled infinite primary tangent angle finding case";}
            }
            elsif (   $outgoing->rot_next_vd()->rot_next_vd()->is_primary()
                   && $outgoing->rot_prev_vd()->rot_prev_vd()->is_primary()) {
                print "CASE P3" if $inf_theta_debug;
                if (   ${$outgoing->rot_next_vd()->rot_next_vd()}
                    == ${$outgoing->rot_prev_vd()->rot_prev_vd()}) {
                    print ".1 " if $inf_theta_debug;
                    $theta = tangent_angle($outgoing->rot_prev_vd()->rot_prev_vd()) + $pi;
                }
                else {
                    print ".2 " if $inf_theta_debug;
                    my $a1 = tangent_angle($outgoing->rot_next_vd()->rot_next_vd());
                    my $a2 = tangent_angle($outgoing->rot_prev_vd()->rot_prev_vd());
                    $theta = $a2 + (abs($a1 - $a2) / 2);
                }
            }
            else {
                die "unhandled infinite primary tangent angle finding case";
            }
                        
            # angle reduce
            if (defined $theta) {
                while ($theta >  $pi) {$theta-=2*$pi;}
                while ($theta < -$pi) {$theta+=2*$pi;}
            }

            # turn it around if outgoing wasn't the edge, but the edge's twin
            if (!$edge->vertex0()) {
                if ($theta >  $pi) {$theta -= $pi}
                if ($theta <= $pi) {$theta += $pi}
            }
        }
        elsif ($edge->is_secondary()) {
            # infinite linear secondary

            if ($inf_theta_debug) {
                warn '[',($edge + 0) , "] ",
                ($edge->vertex0() ? "[".$edge->vertex0()->x().", ".$edge->vertex0()->y()."] "
                                  : "[".$edge->vertex1()->x().", ".$edge->vertex1()->y()."] "),
                $edge->vertex1()?'IN ':'OUT', "\n       ",

                ($edge->vertex1() && $edge->prev_vd()->is_infinite() && $edge->prev_vd()->is_secondary()) || 
                ($edge->vertex0() && $edge->next_vd()->is_infinite() && $edge->next_vd()->is_secondary())
                ? '2NDLINK':'       ', " ",

                ($edge->vertex1() && $edge->prev_vd()->is_infinite() && $edge->prev_vd()->is_primary()) || 
                ($edge->vertex0() && $edge->next_vd()->is_infinite() && $edge->next_vd()->is_primary())
                ? 'PRILINK':'       ', " ",

                $edge->rot_prev_vd()->is_secondary() ? 'ROTP2ND':'ROTPPRI'," ",
                $edge->rot_next_vd()->is_secondary() ? 'ROTN2ND':'ROTNPRI'," ",

                $edge->rot_prev_vd()->is_finite() ? 'ROTPFIN':'ROTPINF'," ",
                $edge->rot_next_vd()->is_finite() ? 'ROTNFIN':'ROTNINF'," ",

                "\n"
                ;
            }


            # This is meant to catch the cases in the Voronoi diagram
            # corresponding to the "fixed" for medial axis where
            # an infinite secondary shoots out past the end of an input segment
            # touching it with it's side, but not stopping there.
            if (   $edge->vertex0()
                   && $edge->rot_next_vd()->is_primary()
                   && $edge->rot_prev_vd()->is_primary()
                  ) {
                print "\nCASE S8\n" if $inf_theta_debug;
                if ($edge->rot_next_vd()->is_finite()) {
                    $theta = tangent_angle($edge->rot_next_vd()) - $pi/2;
                }
                else {
                    $theta = tangent_angle($edge->rot_prev_vd()) + $pi/2;
                }
            }
            elsif (   $edge->vertex1()
                   && $edge->twin()->rot_next_vd()->is_primary()
                   && $edge->twin()->rot_prev_vd()->is_primary()
                  ) {
                print "\nCASE S9\n" if $inf_theta_debug;
                if ($edge->twin()->rot_next_vd()->is_finite()) {
                    $theta = tangent_angle($edge->twin()->rot_next_vd()) - $pi/2;
                }
                else {
                    $theta = tangent_angle($edge->twin()->rot_prev_vd()) + $pi/2;
                }
                $theta += $pi;
            }


            # convex corner (or "flat corner" for collinear input segments)
            elsif (${$edge->next_vd()} == ${$edge->prev_vd()}) {
                print "CASE S1" if $inf_theta_debug;
                if ($edge->vertex0()) {

                    # for Voronoi only (givies infinite recursion for medial axis)
                    if (   $edge->twin()->prev_vd()->rot_prev_vd()->is_primary()
                        && $edge->twin()->prev_vd()->rot_next_vd()->is_primary() 
                    ) {
                        print ".01 " if $inf_theta_debug;
                        $theta = tangent_angle($edge->twin()->prev_vd());
                    }
                    # for medial axis
                    else {
                        print ".1 " if $inf_theta_debug;
                        my $p1 = $edge->twin()->vertex1();
                        my $p2 = $edge->twin()->prev_vd()->vertex0();
                        my $seg_angle = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
                        $theta = $seg_angle + $pi/2;
                    }
                }
                elsif ($edge->vertex1()) {

                    # for Voronoi only (givies infinite recursion for medial axis)
                    if (   $edge->twin()->next_vd()->twin()->rot_prev_vd()->is_primary()
                        && $edge->twin()->next_vd()->twin()->rot_next_vd()->is_primary() 
                    ) {
                        print ".02 " if $inf_theta_debug;
                        $theta = tangent_angle($edge->twin()->next_vd());
                    }
                    # for medial axis
                    else {
                        print ".2 " if $inf_theta_debug;
                        my $p1 = $edge->twin()->vertex0();
                        my $p2 = $edge->twin()->next_vd()->vertex1();
                        my $seg_angle = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
                        $theta = $seg_angle + $pi/2;
                    }
                }
            }

            # infinite secondaries on each end of an input segment
            elsif (   $edge->vertex0() 
                   && $edge->next_vd()->is_secondary()
                  ) {
                print "CASE S2" if $inf_theta_debug;

                # for Voronoi only (gives infinite recursion for medial axis)
                if (   $edge->next_vd()->twin()->rot_prev_vd()->is_primary()
                    && $edge->next_vd()->twin()->rot_next_vd()->is_primary() 
                ) {
                    print ".0 " if $inf_theta_debug;
                    $theta = tangent_angle($edge->next_vd()->twin());
                }
                # for medial axis
                else {
                    print ".1 " if $inf_theta_debug;
                    my $p1 = $edge->vertex0();
                    my $p2 = $edge->next_vd()->vertex1();
                    my $seg_angle = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
                    $theta = $seg_angle - $pi/2;
                }
            }
            elsif (   $edge->vertex1() 
                   && $edge->prev_vd()->is_secondary()
                  ) {
                print "CASE S3" if $inf_theta_debug;

                # for Voronoi only (gives infinite recursion for medial axis)
                if (   $edge->prev_vd()->rot_prev_vd()->is_primary()
                    && $edge->prev_vd()->rot_next_vd()->is_primary() 
                ) {
                    print ".0 " if $inf_theta_debug;
                    $theta = tangent_angle($edge->prev_vd()) + $pi;
                }
                # for medial axis
                else {
                    print ".1 " if $inf_theta_debug;
                    my $p1 = $edge->vertex1();
                    my $p2 = $edge->prev_vd()->vertex0();
                    my $seg_angle = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
                    $theta = $seg_angle - $pi/2;
                }
            }
            

            # infinite secondaries where next or prev is infinite primary
            elsif (   $edge->vertex0() 
                   && $edge->next_vd()->is_primary()
                  ) {
                print "CASE S2P" if $inf_theta_debug;
                
                if ($edge->twin()->prev_vd()->is_secondary()) {

                    # for Voronoi only (might give infinite recursion for medial axis)
                    if (   $edge->twin()->prev_vd()->rot_prev_vd()->is_primary()
                        && $edge->twin()->prev_vd()->rot_next_vd()->is_primary() 
                    ) {
                        print ".0 " if $inf_theta_debug;
                        $theta = tangent_angle($edge->twin()->prev());
                    }
                    # for medial axis
                    else {
                        print ".1 " if $inf_theta_debug;
                        my $p1 = $edge->vertex0();
                        my $p2 = $edge->twin()->prev_vd()->vertex0();
                        my $seg_angle = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
                        $theta = $seg_angle + $pi/2;
                    }
                }
                else {die "crap1";}
                    
            }
            elsif (   $edge->vertex1() 
                   && $edge->prev_vd()->is_primary()
                  ) {
                print "CASE S3P" if $inf_theta_debug;

                if ($edge->twin()->next_vd()->is_secondary()) {

                    # for Voronoi only (might give infinite recursion for medial axis)
                    if (   $edge->twin()->next_vd()->twin()->rot_prev_vd()->is_primary()
                        && $edge->twin()->next_vd()->twin()->rot_next_vd()->is_primary() 
                    ) {
                        print ".0 " if $inf_theta_debug;
                        $theta = tangent_angle($edge->twin()->next());
                    }
                    # for medial axis
                    else {
                        print ".1 " if $inf_theta_debug;
                        my $p1 = $edge->vertex1();
                        my $p2 = $edge->twin()->next_vd()->vertex1();
                        my $seg_angle = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
                        $theta = $seg_angle + $pi/2;
                    }
                }
                else {die "crap2";}
            }

            # a corner where on one side it's infinite and on the
            # other it's finite - because around the corner on the
            # finite side there's convexity - but we're interested in
            # what's round the corner on the infinite side

            elsif (   $edge->vertex1() 
                && $edge->next_vd()->is_secondary()
                && $edge->next_vd()->is_finite()
               ) {
                print "CASE S4" if $inf_theta_debug;
                my $p1 = $edge->twin()->vertex0();
                my $p2 = $edge->twin()->next_vd()->vertex1();
                my $seg_angle = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
                $theta = $seg_angle + $pi/2;
            }
            elsif (   $edge->vertex0() 
                && $edge->prev_vd()->is_secondary()
                && $edge->prev_vd()->is_finite()
               ) {
                print "CASE S5" if $inf_theta_debug;
                my $p1 = $edge->twin()->vertex1();
                my $p2 = $edge->twin()->prev_vd()->vertex0();
                my $seg_angle = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
                $theta = $seg_angle + $pi/2;
            }

            # Similar to "flat corners" we've created for the 
            # medial axis upstream, we've also created corner-like 
            # situations where previously an infinite secondary shot 
            # out past the end of an input segment that was the end 
            # of an open path of segments, rather than a polygon.
            # As long as we're figuring tangent angles for medial axis
            # edges in the C++ code, we shouldn't need this. But it's
            # here in case that upstream figuring ever isn't available.
            elsif (   $edge->vertex1()
                   && $edge->next_vd()->is_secondary()
                   && $edge->next_vd()->is_finite()
                   && $edge->next_vd()->next_vd()->is_primary()
                   #&& $edge->next_vd()->next_vd()->is_infinite()
                  ) {
                print "\nCASE S6\n" if $inf_theta_debug;
                my $p1 = $edge->next_vd()->vertex0();
                my $p2 = $edge->next_vd()->vertex1();
                $theta = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
            }
            elsif (   $edge->vertex0()
                   && $edge->prev_vd()->is_secondary()
                   && $edge->prev_vd()->is_finite()
                   && $edge->prev_vd()->prev_vd()->is_primary()
                   #&& $edge->prev_vd()->prev_vd()->is_infinite()
                  ) {
                print "\nCASE S7\n" if $inf_theta_debug;
                my $p1 = $edge->prev_vd()->vertex0();
                my $p2 = $edge->prev_vd()->vertex1();
                $theta = atan2($p2->y() - $p1->y(), $p2->x() - $p1->x());
            }
            else {
                print "failed to find tangent angle for infinite secondary edge\n";
            }

            # angle reduce
            if (defined $theta) {
                while ($theta >  $pi) {$theta-=2*$pi;}
                while ($theta < -$pi) {$theta+=2*$pi;}
            }
        }
        print " THETA: $theta\n" if $inf_theta_debug;
    }

    return $theta;

}

# Given a finite curved Voronoi edge, calculate the quadratic bezier control
# point between the two end points, because a quadratic Bezier is a section of
# a parabola. 

sub calc_Q {
    my $edge = shift;

    my $edge_theta = tangent_angle($edge);
    my $edge_twin_theta = tangent_angle($edge->twin());

    #warn "thetas: $edge_theta, $edge_twin_theta\n";
    if (! defined($edge_theta) || ! defined($edge_twin_theta)
       || $edge_theta =~ /[\+\-]*nan/ || $edge_twin_theta =~ /[\+\-]*nan/
       ) {
        die "an undef or nan theta: edge:$edge_theta, twin:$edge_twin_theta\n";
        return undef;
    }

    my $q=[];
    
    if (   $edge->is_curved() 
        && $edge->is_finite()
       ) {
        
        my $theta_start = $edge_theta;
        my $theta_end = $edge_twin_theta;
        $theta_end += $pi;
        while ($theta_end >  $pi) {$theta_end = $theta_end - ($pi * 2);}
        my $tan_theta_start = tan($theta_start);
        my $tan_theta_end   = tan($theta_end);
        my $delta_tan = $tan_theta_end - $tan_theta_start;
        
        #if (abs($delta_tan) < 0.0001 ) {
        if (abs($delta_tan) == 0 ) {
            $q->[0]=$edge->vertex0()->x()/2 + $edge->vertex1()->x()/2;
            $q->[1]=$edge->vertex0()->y()/2 + $edge->vertex1()->y()/2;
        }
        else {
            my $y0 = $edge->vertex0()->y();
            my $y2 = $edge->vertex1()->y();
            my $x0 = $edge->vertex0()->x();
            my $x2 = $edge->vertex1()->x();
            my $delta_y0_y2 = $y0 - $y2;
            my $delta_x2_x0 = $x2 - $x0;

            my $x1 = $x0 + ( ( $delta_y0_y2 + $tan_theta_end * $delta_x2_x0 ) / $delta_tan );
            my $y1 = $y0 + $tan_theta_start * ($x1 - $x0);
            
            if ($x1 =~ /[\+\-]*nan/ || $y1 =~ /[\+\-]*nan/) { # get rid of this
                die "oh nan! $x1 , $y1\n";
                return [0,0];
            }
            
            $q->[0]= $x1;
            $q->[1]= $y1;
        }
    }
    else {die "don't arrive here";}      
    
    return $q;

}

}

# == and != operators for vertices, edges, cells
# is identity comparison, useful for traversing
# topology, to see if you've been to been by that
# entity before. These are not any kind of geometric
# coordinate tests.

package Boost::Polygon::Voronoi::Vertex;
{
use overload
    '=='   => sub { ${$_[0]} == ${$_[1]} },
    '!='   => sub { ${$_[0]} != ${$_[1]} },
    'bool' => sub { defined $_[0] }
;
}
package Boost::Polygon::Voronoi::Edge;
{
use overload
    '=='   => sub { ${$_[0]} == ${$_[1]} },
    '!='   => sub { ${$_[0]} != ${$_[1]} },
    'bool' => sub { defined $_[0] }
;
}
package Boost::Polygon::Voronoi::Cell;
{
use overload
    '=='   => sub { ${$_[0]} == ${$_[1]} },
    '!='   => sub { ${$_[0]} != ${$_[1]} },
    'bool' => sub { defined $_[0] }
;
}

package Boost::Polygon::MedialAxis::Vertex;
{
use overload
    '=='   => sub { ${$_[0]} == ${$_[1]} },
    '!='   => sub { ${$_[0]} != ${$_[1]} },
    'bool' => sub { defined $_[0] }
;
}
package Boost::Polygon::MedialAxis::Edge;
{
use overload
    '=='   => sub { ${$_[0]} == ${$_[1]} },
    '!='   => sub { ${$_[0]} != ${$_[1]} },
    'bool' => sub { defined $_[0] }
;
}
package Boost::Polygon::MedialAxis::Cell;
{
use overload
    '=='   => sub { ${$_[0]} == ${$_[1]} },
    '!='   => sub { ${$_[0]} != ${$_[1]} },
    'bool' => sub { defined $_[0] }
;
}


# These will provide glue between dynamically created packages
# that auto-scale coordinates and the original packages.
# The main issue is that when making subroutines that go into
# the symbol table, you can't use $self->SUPER::func();
# So for any coordinates, you need to provide a way to get at
# those $self->SUPER::x() like functions.
# So we can set up $self->super_x() here, and the dynamic subclasses
# can use that instead.
# The re-blessing of cells,edges,vertices can be done based on
# the strings we get from ref($_[0]).

package Boost::Polygon::Voronoi::VertexScaled;
{
push @Boost::Polygon::Voronoi::VertexScaled::ISA, 'Boost::Polygon::Voronoi::Vertex';
sub super_x {$_[0]->SUPER::x();}
sub super_y {$_[0]->SUPER::y();}
sub incident_edge {bless $_[0]->SUPER::incident_edge(), ref($_[0])=~s/::VertexScaled/::EdgeScaled/r;}
sub cell {bless $_[0]->SUPER::cell(), ref($_[0])=~s/::VertexScaled/::CellScaled/r;}
}

package Boost::Polygon::Voronoi::EdgeScaled;
{
push @Boost::Polygon::Voronoi::EdgeScaled::ISA, 'Boost::Polygon::Voronoi::Edge';
sub twin     {bless $_[0]->SUPER::twin()    , ref($_[0]);}
sub next     {bless $_[0]->SUPER::next()    , ref($_[0]);}
sub prev     {bless $_[0]->SUPER::prev()    , ref($_[0]);}
sub rot_next {bless $_[0]->SUPER::rot_next(), ref($_[0]);}
sub rot_prev {bless $_[0]->SUPER::rot_prev(), ref($_[0]);}
sub vertex0  {my $v = $_[0]->SUPER::vertex0(); !ref($v) ? undef : bless $v , ref($_[0])=~s/::EdgeScaled/::VertexScaled/r;}
sub vertex1  {my $v = $_[0]->SUPER::vertex1(); !ref($v) ? undef : bless $v , ref($_[0])=~s/::EdgeScaled/::VertexScaled/r;}
sub cell     {bless $_[0]->SUPER::cell()    , ref($_[0])=~s/::EdgeScaled/::CellScaled/r;}
}

package Boost::Polygon::Voronoi::CellScaled;
{
push @Boost::Polygon::Voronoi::CellScaled::ISA, 'Boost::Polygon::Voronoi::Cell';
sub incident_edge {bless $_[0]->SUPER::incident_edge(), ref($_[0])=~s/::CellScaled/::EdgeScaled/r;}
}

# repeat for MedialAxis
package Boost::Polygon::MedialAxis::VertexScaled;
{
push @Boost::Polygon::MedialAxis::VertexScaled::ISA, 'Boost::Polygon::MedialAxis::Vertex';
sub super_x {$_[0]->SUPER::x()}
sub super_y {$_[0]->SUPER::y()}
sub super_r {$_[0]->SUPER::r()}
sub incident_edge {bless $_[0]->SUPER::incident_edge(), ref($_[0])=~s/::VertexScaled/::EdgeScaled/r;}
sub cell {bless $_[0]->SUPER::cell(), ref($_[0])=~s/::VertexScaled/::CellScaled/r;}
}

package Boost::Polygon::MedialAxis::EdgeScaled;
{
push @Boost::Polygon::MedialAxis::EdgeScaled::ISA, 'Boost::Polygon::MedialAxis::Edge';
sub twin     {bless $_[0]->SUPER::twin()    , ref($_[0]);}
sub next     {bless $_[0]->SUPER::next()    , ref($_[0]);}
sub prev     {bless $_[0]->SUPER::prev()    , ref($_[0]);}
sub rot_next {bless $_[0]->SUPER::rot_next(), ref($_[0]);}
sub rot_prev {bless $_[0]->SUPER::rot_prev(), ref($_[0]);}
sub vertex0  {my $v = $_[0]->SUPER::vertex0(); !ref($v) ? undef : bless($v , ref($_[0])=~s/::EdgeScaled/::VertexScaled/r);}
sub vertex1  {my $v = $_[0]->SUPER::vertex1(); !ref($v) ? undef : bless($v , ref($_[0])=~s/::EdgeScaled/::VertexScaled/r);}
sub cell     {bless $_[0]->SUPER::cell()    , ref($_[0])=~s/::EdgeScaled/::CellScaled/r;}
sub foot     {bless $_[0]->SUPER::foot()    , ref($_[0])=~s/::EdgeScaled/::FootScaled/r;}
sub Q        {bless $_[0]->SUPER::Q()       , ref($_[0])=~s/::EdgeScaled/::FootScaled/r;}
}

package Boost::Polygon::MedialAxis::CellScaled;
{
push @Boost::Polygon::MedialAxis::CellScaled::ISA, 'Boost::Polygon::MedialAxis::Cell';
sub incident_edge {bless $_[0]->SUPER::incident_edge(), ref($_[0])=~s/::CellScaled/::EdgeScaled/r;}
}

package Boost::Polygon::MedialAxis::FootScaled;
{
push @Boost::Polygon::MedialAxis::FootScaled::ISA, 'Boost::Polygon::MedialAxis::Foot';
sub super_x {$_[0]->SUPER::x()}
sub super_y {$_[0]->SUPER::y()}
}

1;

__END__
=pod

=head1 NAME

Boost::Polygon::Voronoi - Bindings for the Boost Polygon Voronoi library

=head1 VERSION

version 0.01

=head1 SYNOPSIS

    use Boost::Polygon::Voronoi; # etc.

=head1 ABSTRACT

Bindings for the C++ Boost Polygon Voronoi library, with an additional Medial Axis implementation.

=head1 METHODS

=head2 addPolygon

addPolygon

=head2 addPolygons

addPolygons

=head2 addPolyline

addPolyline

=head2 addPolylines

addPolylines

=head2 addSegment

addSegment

=head2 addSegments

addSegments

=head2 addPoint

addPoint

=head2 addPoints

addPoints

=head2 validateSegments

validateSegments

=head2 calc_Q

calc_Q

=head2 tangent_angle

tangent_angle

=head2 get_curve_focus

get_curve_focus

=for Pod::Coverage addPolygon addPolygons addPolyline addPolylines addSegment addSegments addPoint addPoints validateSegments calc_Q tangent_angle get_curve_focus

=head1 ACKNOWLEDGEMENTS

Thanks to Andrii Sydorchuk for the original C++ Boost::Polygon::Voronoi.

=head1 AUTHOR

Michael E. Sheldrake

=head1 COPYRIGHT AND LICENSE

This software is copyright (c) 2018 by Michael E. Sheldrake.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

=cut

