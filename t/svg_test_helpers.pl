use warnings;
use Test::TempDir::Tiny;

my $pi = 4*atan2(1.0,1.0);
sub tan {return sin($_[0])/cos($_[0]);}

sub get_voronoi_diagram {
    my $builder = shift;
    my $vd = $builder->voronoi_diagram();
    verify_vd($vd);
    return $vd;
}

sub get_medial_axis {
    my $builder = shift;
    my $ma = $builder->medial_axis();
    verify_vd($ma);
    verify_ma($ma);
    return $ma;
}

# should be safe to run on both voronoi diagram and medial axis
sub verify_vd {
    my $vd = shift;
    my @edges = @{$vd->edges()};
    ok(scalar(@edges) == $vd->num_edges(), "got edges");
    my @vertices = @{$vd->vertices()};
    ok(scalar(@vertices) == $vd->num_vertices(), "got vertices");
    my @cells = @{$vd->cells()};
    ok(scalar(@cells) == $vd->num_cells(), "got cells");
}

sub verify_ma {
    my $ma = shift;
}

sub load_test_svg {
    my $builder = shift;
    # load test geometry from svg
    my $svgfn = $0;
    $svgfn.='.svg';

    ok(-e $svgfn, "svg file for test exists");

    my $svg='';
    if (open my $svgfh, '<', $svgfn) {$svg=join('',<$svgfh>);close $svgfh;}

    ok(length($svg) > 0, "load svg file");

    my @test_paths = $svg=~/<path[^>]* class="test_path"[^>]*\/>/msg;
    my @test_ds = map {/ d="([^"]+)"/;$1} @test_paths;

    my $test_scale = 1;
    if ($svg=~/<([^ >]+)\s+[^>]* id="test_scale"[^>]*>(.*?)<\/\1>/ms) {
        $test_scale=$2;
        while ($test_scale=~/<([^ >]+)[^>]*>(.*?)<\/\1>/ms) {$test_scale=$2;}
        $test_scale =~ s/\s*(.*?)\s*/$1/ms;
        $test_scale = eval($test_scale);
        $test_scale=1 if !$test_scale;
    }
    
    ok(@test_ds > 0, "found test paths");
    
    my @builder_path_segs;

    my $bad_path_index=-1;
    my $bad_path_report='';
    my $segcount = 0;

    foreach (my $pi=0;$pi < @test_ds; $pi++) {
        # skip e in split because sometimes it's in numbers like 1.0e-12
        my @segs = split(/([a-df-zA-DF-Z])/, $test_ds[$pi]);
        shift @segs if $segs[0] =~ /^\s*$/;
        # give any Z commands coordinates from preceeding M
        my $curr_m_coords = $segs[1];
        for (my $i=0;$i<@segs;$i+=2) {
            if ($segs[$i]=~/M/) { $curr_m_coords = $segs[$i+1]; }
            elsif ($segs[$i]=~/Z/) { $segs[$i+1] = $curr_m_coords; }
        }
        # check that path is all absolute linear moves
        if (scalar(@segs) != scalar((grep {($_ % 2 == 0 ? ($segs[$_]=~/[MLHVZ]/) : ($segs[$_]=~/[0-9eE\.\- ,]+/))} (0..$#segs)))) {
            $bad_path_index = $pi;
            $bad_path_report = join('.',map {$_ % 2 == 0 ? $segs[$_].'['.($segs[$_]=~/[MLHVZ]/ ? 'ok':'bad').']' : $segs[$_].'['.($segs[$_]=~/[0-9eE\.\- ,]+/ ? 'ok':'bad').']'} (0..$#segs));
            last;
        }
        # turn path into segments ready for builder->insert_segment()
        my @these_builder_path_segs;
        my $last_x;
        my $last_y;
        for (my $i=0;$i<@segs;$i+=2) {
            my $cmd = $segs[$i];
            my $coords_string=$segs[$i+1];
            $coords_string=~s/^\s*(.*?)\s*$/$1/;
            my @coords = map {eval($_) + 0} split(/[ ,]+/, $coords_string);
            if ($cmd eq 'M') {($last_x,$last_y)=($coords[0],$coords[1]);}
            if ($cmd eq 'L' || $cmd eq 'Z') {
                for (my $j=0;$j<@coords;$j+=2) {
                    push @these_builder_path_segs, [$last_x,$last_y,$coords[$j],$coords[$j+1]];
                    ($last_x,$last_y) = ($coords[$j],$coords[$j+1]);
                }
            }
            if ($cmd eq 'H') {
                for (my $j=0;$j<@coords;$+=1) {
                    push @these_builder_path_segs, [$last_x,$last_y,$coords[$j],$last_y];
                    $last_x = $coords[$j];
                }
            }
            if ($cmd eq 'V') {
                for (my $j=0;$j<@coords;$+=1) {
                    push @these_builder_path_segs, [$last_x,$last_y,$last_x,$coords[$j]];
                    $last_x = $coords[$j];
                }
            }
        }
        $segcount += @these_builder_path_segs;
        for my $s (@these_builder_path_segs) {$_*=$test_scale for @$s;}
        push @builder_path_segs, \@these_builder_path_segs;
    }

    ok($bad_path_index < 0, "paths are all absolute linear segments".($bad_path_index > -1?" (bad path[$bad_path_index]: ".$test_ds[$bad_path_index].")\n$bad_path_report":''));
        # make sure path is only absolute linear moves M,L,H,V,Z (z okay too I guess)
    
    ok(scalar(@builder_path_segs) == scalar(@test_ds), "built segment sequences for all test paths");


    my $okay_segsets=0;
    $okay_segsets += (Boost::Polygon::Voronoi::validateSegments($_) ? 1 : 0) for @builder_path_segs;

    ok($okay_segsets == scalar(@builder_path_segs), "all segment sets valid ($okay_segsets of ".scalar(@builder_path_segs).")");


    $bad_path_index = -1;
    my $bad_seg_report;
    my $insert_last_index = 0;
    for (my $i=0;$i<@builder_path_segs;$i++) {
        for (my $j=0;$j<@{$builder_path_segs[$i]};$j++) {
            eval {
                $insert_last_index = $builder->insert_segment(@{$builder_path_segs[$i]->[$j]});
                #warn "\n[",join(",",@{$builder_path_segs[$i]->[$j]}),"]\n";
            };
            if ($@) { 
                $bad_seg_report = 'trouble inserting this segment: ['.
                                  join(',',@{$builder_path_segs[$i]->[$j]}).
                                  ']'.$@;
                last;
            }
        }
        if ($bad_seg_report) { last; }
    }
    
    ok( ! $bad_seg_report, "no bad segments during insert".($bad_seg_report ? ' '.$bad_seg_report : ''));

    ok(($insert_last_index + 1) == $segcount, "inserted all segments into builder (inserted ".($insert_last_index + 1)." of $segcount)");

    my @test_point_tags = $svg=~/<circle[^>]* class="test_point"[^>]*\/>/msg;
    my @test_points;
    foreach my $point_tag (@test_point_tags) {
        my $cx;
        if ($point_tag =~ /cx="([^"]+)"/) {$cx=$1;}
        my $cy;
        if ($point_tag =~ /cy="([^"]+)"/) {$cy=$1;}
        push @test_points, [$cx,$cy] if (defined $cx && defined $cy);
    }

    ok(scalar(@test_points) == scalar(@test_point_tags), "found all test point coordinates (".scalar(@test_points)."/".scalar(@test_point_tags).")");

    $bad_point_index = -1;
    my $bad_point_report;
    my $points_insert_last_index = $insert_last_index; # didn't need to do that
    for (my $i=0;$i<@test_points;$i++) {
        eval {
            $points_insert_last_index = $builder->insert_point(@{$test_points[$i]});
            #warn "\n[",join(",",@{$test_points[$i]}),"]\n";
        };
        if ($@) { 
            $bad_point_report = 'trouble inserting this point: ['.
                              join(',',@{$test_points[$i]}).
                              ']'.$@;
            last;
        }
    if ($bad_point_report) { last; }
    }

    ok( ! $bad_point_report, "no bad points during insert" . ($bad_point_report ? ' ('.$bad_point_report.')' : '' ));

    ok((($points_insert_last_index - $insert_last_index)) == scalar(@test_points), "inserted all points into builder (inserted ".(($points_insert_last_index - $insert_last_index))." of ".scalar(@test_points).")");

}


sub make_result_svg {
    my ($edges,$cells,$vertices,$stroke_width_factor,$vd_or_ma,$out_dir) = @_;
    
    my $isma = $vd_or_ma =~ /medial/i ? 1 : 0;

    my $svg_results = '';

    # load test geometry from svg
    my $svgfn    = $0.'.svg';
    my $svgoutfn;
    my $scriptdir='';
    my $slash = '/';
    if ($svgfn =~ /^(.*?)([\\\/])([^\\\/]+)$/) {$scriptdir=$1;$slash=$2;$svgoutfn = $3;}
    else {$svgoutfn = $svgfn;}
    $svgoutfn .= '.result.'.($isma?'ma':'vd').'.svg';
    my $tmpdir = tempdir("svg_$vd_or_ma");
    $ENV{PERL_TEST_TEMPDIR_TINY_NOCLEANUP}=1;
    $svgoutfn = ($out_dir ? $out_dir : $tmpdir) . $slash . $svgoutfn;


    ok(-e $svgfn, "svg file for test results exists");

    my $svg='';
    if (open my $svgfh, '<', $svgfn) {$svg=join('',<$svgfh>);close $svgfh;}

    ok(length($svg) > 0, "load svg file for test results");

    my $test_scale = 1;
    if ($svg=~/<([^ >]+)\s+[^>]* id="test_scale"[^>]*>(.*?)<\/\1>/ms) {
        $test_scale=$2;
        while ($test_scale=~/<([^ >]+)[^>]*>(.*?)<\/\1>/ms) {$test_scale=$2;}
        $test_scale =~ s/\s*(.*?)\s*/$1/ms;
        $test_scale = eval($test_scale);
        $test_scale=1 if !$test_scale;
    }


    #put this in a style section modifying the style from the external stylesheet.
    #my $basic_style = "stroke-width:".(0.03 * $stroke_width_factor).";";
    my $sw = 0.3 * ($stroke_width_factor || 1);
    my $voronoi_path_styles = "stroke-opacity:0.5;fill:none;stroke:gray;stroke-width:$sw;";
    $svg=~s/(\.vp\s*\{[^\}]*)\}/$1;$voronoi_path_styles}/;
    my $voronoi_cell_styles = "stroke-opacity:1;fill:yellow;fill-opacity:0.3;stroke:white;stroke-width:".($sw*1.1).";";
    $svg=~s/(\.vc\s*\{[^\}]*)\}/$1;$voronoi_cell_styles}/;
    
    my $edges_all = '';


    my $a_point = scalar(@$edges) > 0
                ? ($edges->[0]->vertex0() || $edges->[0]->vertex1())
                : (scalar(@$vertices) > 0 
                   ? $vertices->[0] 
                   : ( (@$cells && scalar(grep {!$_->is_degenerate()} @$cells))
                       ? (grep {!$_->is_degenerate()} @$cells)[0]->incident_edge()->vertex0() 
                       : undef
                     )
                  );

    if (!$a_point) {die "couldn't find an initial vertex for bounds finding";}
    my ($minx,$maxx) = ($a_point->x(),$a_point->x());
    my ($miny,$maxy) = ($a_point->y(),$a_point->y());

    if ($edges) {
        for (@$edges) {
            for (grep $_, ($_->vertex0(),$_->vertex1())) {
                if ($minx > $_->x()) {$minx = $_->x();}
                if ($miny > $_->y()) {$miny = $_->y();}
                if ($maxx < $_->x()) {$maxx = $_->x();}
                if ($maxy < $_->y()) {$maxy = $_->y();}
            }
        }
    }

    my $bbox = {minx=>$minx,miny=>$miny,maxx=>$maxx,maxy=>$maxy,xspan=>($maxx-$minx),yspan=>($maxy-$miny)};
    
    # medial axis
    if ($isma) {

        # black background, to help with inside/outside display
        $svg =~ s/pagecolor="#[^"]+"/pagecolor="#111111"/;

        #my @zero_cells = grep $_->color() == 0, @$cells;        
        #$edges_all .= cell_set_svg(\@zero_cells,$test_scale,$sw,$bbox);
        #my @one_cells = grep $_->color() == 1, @$cells;        
        #$edges_all .= cell_set_svg(\@one_cells,$test_scale,$sw,$bbox);

        $edges_all .= edge_loops_svg(\@$edges,$test_scale,$sw,$bbox);
        my @edges_internal_primary   = grep $_->is_internal() && $_->is_primary(), @$edges;
        my @edges_external_primary   = grep $_->is_external() && $_->is_primary(), @$edges;
        my @edges_internal_secondary = grep $_->is_internal() && $_->is_secondary(), @$edges;
        my @edges_external_secondary = grep $_->is_external() && $_->is_secondary(), @$edges;

        my @edges_limbo_secondary   = grep !$_->is_external() && !$_->is_internal() && $_->is_secondary, @$edges;
        my @edges_limbo_primary   = grep !$_->is_external() && !$_->is_internal() && $_->is_primary(), @$edges;

        $edges_all .= edge_set_svg(\@edges_external_secondary,'external secondary','stroke:orange;color:orange;',$bbox,$test_scale,$sw);
        $edges_all .= edge_set_svg(\@edges_external_primary,  'external primary'  ,'stroke:red;color:blue;'   ,$bbox,$test_scale,$sw);
        $edges_all .= edge_set_svg(\@edges_internal_secondary,'internal secondary','stroke:green;color:blue;' ,$bbox,$test_scale,$sw);
        $edges_all .= edge_set_svg(\@edges_internal_primary,  'internal primary'  ,'stroke:blue;color:blue;'  ,$bbox,$test_scale,$sw);
        
        $edges_all .= edge_set_svg(\@edges_limbo_secondary,  'limbo 2ndary'  ,'stroke:aqua;'  ,$bbox,$test_scale,$sw);
        $edges_all .= edge_set_svg(\@edges_limbo_primary,  'limbo primary'  ,'stroke:aqua;'  ,$bbox,$test_scale,$sw);

    }
    # voronoi diagram
    else {

        $edges_all .= cell_set_svg($cells,$test_scale,$sw,$bbox);

        my @edges_primary   = grep $_->is_primary(), @$edges;
        my @edges_secondary = grep $_->is_secondary(), @$edges;
        $edges_all .= edge_set_svg(\@edges_secondary,'secondary','stroke:green;color:green;' ,$bbox,$test_scale,$sw);
        $edges_all .= edge_set_svg(\@edges_primary,  'primary'  ,'stroke:blue;color:blue;'  ,$bbox,$test_scale,$sw);

    }



    $svg_results .= $edges_all;
    
    # vd just has edges, verts, cells, primary/secondary, infinite?, curved?
    # for stock vd, will have to calc Q for curved edges
    
    # then vd adds internal/external

    # have a way to represent infinte rays
    

    # then maybe some common errors
    # Q that's outside of convex hull of all input geometry - make that Q just a midpoint, but give it a big red circle and put debug info in that circle's attributes
    # - missing foot


    my $substitute=0;
    
    $svg =~ s/(<g[^>]* id="test result"[^>]*?)(\/>)/$1><\/g>/ms;
    $substitute++ if ($svg =~ s/(<g[^>]* id="test result"[^>]*?>)/$1$svg_results/ms);
    ok($substitute, "inserted results");

    my $wrote=0;
    if (open my $fh, '>',$svgoutfn) {
        print $fh $svg;
        close($fh);
        $wrote=1;
    }
    ok($wrote, "svg result file saved");

}
sub edge_set_svg {
    my ($edges, $label, $style, $bbox, $scale, $strokewidth) = @_;
    my $svg = '';

    my $bbdiag=sqrt($bbox->{xspan}**2 + $bbox->{yspan}**2);
    
    #if (!$style) {$style="stroke-width:0.3;stroke:gray;stroke-opacity:0.5;fill:none;";}    
    $svg.='<g inkscape:groupmode="layer" inkscape:label="'.$label.'" id="'.$label.'">'."\n";

    my $indicator_circle_radius = $strokewidth * 2;
    my $q_cir_rad = $strokewidth * 2;

    foreach my $edge (@$edges) {


$svg.='<g>'."\n";

        eval {$edge->is_external()};
        my $isma = $@ ? 0:1;
        # TODO: handle infinite edges where one of $_->vertex0() && $_->vertex1() not defined
        # figure Q for curved

        my $naughty_q = '';        
        
        # infinite edges
        if (!$edge->vertex0() || !$edge->vertex1()) {
            if (!$edge->vertex0() && !$edge->vertex1()) {
                # never seen it. shouldn't happen.
                warn "edge with no vertices [",($$edge),"]";
                next;
            }

            my $theta;
            
            if ($isma && $edge->theta() < 4) { # theta == 4 indicates undefined, for now
                #warn "MA theta defined for infinite edge. Using that. [",$edge->theta(),"]\n";
                $theta = $edge->theta();
            }
            else {
                $Boost::Polygon::Voronoi::last_tangent_angle_case=-1;
                $theta = tangent_angle($edge);  # work it out here then add to C stuff for ma
                $Boost::Polygon::Voronoi::last_tangent_angle_case=-1;
            }
            
            if (defined $theta) {
                my $v0 = $edge->vertex0()
                       ? [$edge->vertex0()->x()/$scale, 
                          $edge->vertex0()->y()/$scale]
                       : [$edge->vertex1()->x()/$scale - ($bbdiag*0.2/$scale) * cos($theta), 
                          $edge->vertex1()->y()/$scale - ($bbdiag*0.2/$scale) * sin($theta)];
                my $v1 = $edge->vertex1() 
                       ? [$edge->vertex1()->x()/$scale, 
                          $edge->vertex1()->y()/$scale] 
                       : [$edge->vertex0()->x()/$scale + ($bbdiag*0.2/$scale) * cos($theta), 
                          $edge->vertex0()->y()/$scale + ($bbdiag*0.2/$scale) * sin($theta)];

                $svg .= '<path class="vp inf" id="'.($$edge).'" previd="'.(${$edge->prev()}).'" nextid="'.(${$edge->next()}).'" style="'.$style.'" ';
                $svg .= 'marker-end="url(#inf_edge_arrow)" ';
                $svg .= 'd="' . sprintf("M%.6f,%.6f",$v0->[0],$v0->[1]);
                $svg .=         sprintf("L%.6f,%.6f",$v1->[0],$v1->[1]);
                $svg .= '" />'."\n";

                if ($isma && $edge->foot()) {
                    $svg .= '<circle id="'.($$edge).'_foot" cx="'.($edge->foot()->x()/$scale).'" cy="'.($edge->foot()->y()/$scale).'" ';
                    $svg .= 'r="'.$q_cir_rad.'" ';
                    $svg .= 'style="fill:green" />'."\n";
                }

            } else {warn "theta not defined, can't show infinite ray\n";}


        }

        # finite edges
        else {
            my $q;
            
            my $midx = ($edge->vertex1()->x() + $edge->vertex0()->x())/2;
            my $midy = ($edge->vertex1()->y() + $edge->vertex0()->y())/2;
            my $length = sqrt(($edge->vertex1()->x() - $edge->vertex0()->x())**2 + ($edge->vertex1()->y() - $edge->vertex0()->y())**2);
            
            # curved - get Q: Bezier control point between end points
            # medial axis has this calculated in C++ code.
            # for stock Voronoi, we calculate it in Perl, for now

            if ($edge->is_curved()) {
                if ($isma && $edge->Q()) {
                    $q = [$edge->Q()->x(),$edge->Q()->y()];
                }
                else { $q = calc_Q($edge);}

                if (!$q) {
                    $q=[$midx,$midy];
                    $naughty_q  = '<circle cx="'.($q->[0]/$scale).'" cy="'.($q->[1]/$scale).'" ';
                    $naughty_q .= 'r="'.$q_cir_rad.'" ';
                    $naughty_q .= 'style="fill:red" error="Couldn\'t find or calculate control point for curved edge."/>'."\n";
                }
                else {
                    my $q_dist = sqrt(($q->[0] - $midx)**2 + ($q->[1] - $midy)**2);
                    if ($q_dist > 10*$length) {
                        $naughty_q  = '<circle cx="'.($q->[0]/$scale).'" cy="'.($q->[1]/$scale).'" ';
                        $naughty_q .= 'r="'.$q_cir_rad.'" ';
                        $naughty_q .= 'style="fill:red" q="'.join(',',@$q).'" midpt="'.join(',',$midx,$midy).'" dist-q-m="'.$q_dist.'" error="Control point for curved edge seemed too far away from the endpoints."/>'."\n";
                    }

                    $naughty_q  = '<circle id="'.($$edge).'_Q" cx="'.($q->[0]/$scale).'" cy="'.($q->[1]/$scale).'" ';
                    $naughty_q .= 'r="'.$q_cir_rad.'" ';
                    $naughty_q .= 'style="fill:gray" q="'.join(',',@$q).'" midpt="'.join(',',$midx,$midy).'" dist-q-m="'.$q_dist.'" note="Control point for curved edge."/>'."\n";

                    if (0 && !$isma) {
                        my $f = get_curve_focus($edge);
                        if ($f) {
                            $naughty_q .= '<circle id="'.($$edge).'_focus" cx="'.($f->x()/$scale).'" cy="'.($f->y()/$scale).'" ';
                            $naughty_q .= 'r="'.($q_cir_rad*10).'" ';
                            $naughty_q .= 'style="fill:aqua" q="'.join(',',@$q).'" midpt="'.join(',',$midx,$midy).'" dist-q-m="'.$q_dist.'" note="Control point for curved edge."/>'."\n";
                        }
                    }
                }
            }
            

            $svg .= '<path class="vp" id="'.($$edge).'" previd="'.(${$edge->prev()}).'" nextid="'.(${$edge->next()}).'" style="'.$style.'" ';
            $svg .= ' color="silver" marker-end="url(#edge_arrow)" '; # color doesn't work for fill of marker here
            $svg .= 'd="'.sprintf("M%.6f,%.6f",($edge->vertex0()->x()/$scale),($edge->vertex0()->y()/$scale));
            $svg .= ($edge->is_curved()
                     ? sprintf("Q%.6f,%.6f,%.6f,%.6f",($q->[0]/$scale),($q->[1]/$scale),($edge->vertex1()->x()/$scale),($edge->vertex1()->y()/$scale))
                     : sprintf("L%.6f,%.6f",($edge->vertex1()->x()/$scale),($edge->vertex1()->y()/$scale))
                    );
            $svg .= '" />'."\n";


            if ($isma && $edge->foot()) {
                $svg .= '<circle id="'.($$edge).'_foot" cx="'.($edge->foot()->x()/$scale).'" cy="'.($edge->foot()->y()/$scale).'" ';
                $svg .= 'r="'.$q_cir_rad.'" ';
                $svg .= 'style="fill:green" />'."\n";
            }


        }



        # show tangents at edge ends, shooting off at angle we're calling theta.
        # try to show tangents for infinite edges too
        if (    
            0
           ) {

            my $theta;
            my $theta_twin;
            
            if (0 && $isma && $edge->theta() < 4) { # theta == 4 indicates undefined, for now
                #warn "MA theta defined for edge. Using that. [",$edge->theta(),"]\n";
                $theta = $edge->theta();
            }
            else {
                $Boost::Polygon::Voronoi::last_tangent_angle_case=-1;
                $theta      = tangent_angle($edge);
                $Boost::Polygon::Voronoi::last_tangent_angle_case=-1;
            }
            if (0 && $isma && $edge->twin()->theta() < 4) { # theta == 4 indicates undefined, for now
                #warn "MA theta defined for edge twin. Using that. [",$edge->twin()->theta(),"]\n";
                $theta_twin = $edge->twin()->theta();
            }
            else {
                $Boost::Polygon::Voronoi::last_tangent_angle_case=-1;
                $theta_twin = tangent_angle($edge->twin());
                $Boost::Polygon::Voronoi::last_tangent_angle_case=-1;
            }
            
            if (defined $theta) {
                my $v0 = $edge->vertex0() 
                       ? [$edge->vertex0()->x()/$scale, 
                          $edge->vertex0()->y()/$scale] 
                       : [$edge->vertex1()->x()/$scale - 200 * cos($theta), 
                          $edge->vertex1()->y()/$scale - 200 * sin($theta)];
                $svg .= '<path class="vp raypoint" is="'.($edge->is_primary()?'PRI':'2ND').' '.($edge->is_curved()?'CRV':'LIN').'" style="'.$style.';stroke-width:'.($strokewidth*5).';stroke:orange;" uid="'.($$edge).'" ';
                $svg .= 'd="'.sprintf("M%.6f,%.6f",$v0->[0],$v0->[1]);
                $svg .= sprintf("L%.6f,%.6f",$v0->[0] + 50*cos($theta),$v0->[1] + 50*sin($theta));
                $svg .= '" />'."\n";
            } else {warn "tangent theta not defined\n";}
            if (defined $theta_twin) {            
                my $v0 = $edge->twin()->vertex0() 
                       ? [$edge->twin()->vertex0()->x()/$scale, 
                          $edge->twin()->vertex0()->y()/$scale] 
                       : [$edge->twin()->vertex1()->x()/$scale - 200 * cos($theta_twin), 
                          $edge->twin()->vertex1()->y()/$scale - 200 * sin($theta_twin)];
                $svg .= '<path class="vp raypoint" is="'.($edge->twin()->is_primary()?'PRI':'2ND').' '.($edge->twin()->is_curved()?'CRV':'LIN').'" style="'.$style.';stroke-width:'.($strokewidth*5).';stroke:yellow;" uid="'.($$edge).'" ';
                $svg .= 'd="'.sprintf("M%.6f,%.6f",$v0->[0],$v0->[1]);
                $svg .= sprintf("L%.6f,%.6f",$v0->[0] + 50*cos($theta_twin),$v0->[1] + 50*sin($theta_twin));
                $svg .= '" />'."\n";
            } else {warn "tangent twin theta not defined\n";}
        }



        $svg .= $naughty_q if $naughty_q;
        $svg .= $infinite_indicator if $infinite_indicator;
        
$svg.='</g>'."\n";

    }

    $svg.='</g>';

    return $svg;
}

sub cell_set_svg {
    my $cells = shift;
    my $scale = shift;
    my $strokewidth=shift;
    my $bbox=shift;

    eval {$edge->is_external()};
    my $isma = $@ ? 0:1;

    my $svg = '';
    
    for (my $i = 0;$i<@$cells;$i++) {
        my $cell = $cells->[$i];
        next if $cell->is_degenerate();
        my $iedge = $cell->incident_edge();
        my $e=$iedge;

        my $show_incident = 0;
        $svg.='<g>'."\n" if $show_incident; # and closed below, if
        if ( $show_incident ) {
            my $edge=$iedge;
            my $note = '['.($$edge).']'.($edge->is_infinite()?'infinite':'NON-INFINITE').' '.($edge->is_primary()?'primary':'secondary');
            my $x = $edge->vertex0() ? $edge->vertex0()->x() : $edge->vertex1()->x();
            my $y = $edge->vertex0() ? $edge->vertex0()->y() : $edge->vertex1()->y();
            my $incident_edge_v0 = '<circle class="incident" ';
            $incident_edge_v0 .= 'thinksItsInfinite="'.( $edge->is_infinite() ? 'YES' : 'NO' ).'" ';
            $incident_edge_v0 .= 'cx="'.($x/$scale).'" cy="'.($y/$scale).'" ';
            $incident_edge_v0 .= 'r="7" ';
            $incident_edge_v0 .= 'marked-external="'.($edge->is_external() ? 'YES':'NO').'" ' if $isma;
            $incident_edge_v0 .= 'style="stroke:green;stroke-width:2;" ';
            $incident_edge_v0 .= 'note="'.$note.'" ';
            $incident_edge_v0 .= 'id="'.($$edge).'" ';
            $incident_edge_v0 .= 'prev="'.(${$edge->prev()}).'" ';
            $incident_edge_v0 .= 'next="'.(${$edge->next()}).'" ';
            $incident_edge_v0 .= '/>'."\n";
            $svg .= $incident_edge_v0;
        }


        my @edges;
        my $cnt=0;

        do {
            push @edges, $e;
            $e = $e->next();
        } while($$e != $$iedge && $cnt++ < 300);
        
        warn "is that a too-many-sided cell?? $cnt" if $cnt > 299;
        
        $svg.='<path class="vc" style="fill:#'.(join '', map {sprintf("%02x",$_*255)} (rand(),rand(),rand())).';" ';

        #$svg.='M'.join('L',map {$_->[0].','.$_->[1]} @poly);
        if (!$edges[0]->vertex0()) {
            $svg .= 'd="'.sprintf("M%.6f,%.6f",($edges[0]->vertex1()->x()/$scale),($edges[0]->vertex1()->y()/$scale));
        }
        else {
            $svg .= 'd="'.sprintf("M%.6f,%.6f",($edges[0]->vertex0()->x()/$scale),($edges[0]->vertex0()->y()/$scale));
        }
        foreach my $edge (@edges) {
            if ($edge->is_curved()) {
                my $q = $isma ? [$edge->Q()->x(),$edge->Q()->y()] : calc_Q($edge);
                $svg .=  sprintf("Q%.6f,%.6f,%.6f,%.6f",$q->[0]/$scale, ($q->[1]/$scale),($edge->vertex1()->x()/$scale),($edge->vertex1()->y()/$scale));
            }
            else {
                if (!$edge->vertex1()) { # then line to where we already are.
                    $svg .= sprintf("L%.6f,%.6f",($edge->vertex0()->x()/$scale),($edge->vertex0()->y()/$scale));
                    }
                else {
                    $svg .= sprintf("L%.6f,%.6f",($edge->vertex1()->x()/$scale),($edge->vertex1()->y()/$scale));
                }
            }

        }
        $svg .= '" />'."\n";

        $svg.='</g>'."\n" if $show_incident;

    }
    return $svg;
}


sub edge_loops_svg {
    my $edges = shift;
    my $scale = shift;
    my $strokewidth=shift;
    my $bbox=shift;

    eval {$edges->[0]->is_external()};
    my $isma = $@ ? 0:1;

    my $svg = '';
    
    my %seen;
    
    for (my $i = 0; $i < @$edges; $i++) {
        next if defined $seen{''.${$edges->[$i]}};

        my $d='';

        $svg.='<path class="vc" ';
        $svg.='style="';
        $svg.='fill-rule:non-zero;';
        if ($edges->[$i]->is_internal()) {
            my $b=255;my $g=128 + 127*rand();my $r=$g;
            $svg.='fill-opacity:0.7;fill:#'.(join '', map {sprintf("%02x",int($_))} ($r,$g,$b)).';';
        }
        elsif ($edges->[$i]->is_external()) {
            my $r=125;my $g=$r*rand();my $b=$g;
            $svg.='fill-opacity:0.7;fill:#'.(join '', map {sprintf("%02x",int($_))} ($r,$g,$b)).';';
        }
        else {            
            $svg.='fill-opacity:0.7;fill:aqua;';
        }
        $svg.='stroke:yellow;stroke-width:'.($strokewidth).';" ';
        $svg.='d="';


        my $iedge = $edges->[$i];
        my $e=$iedge;

        my @loop_edges;
        my @poly_verts;
        my $cnt=0;

        # any combo of these two options should work
        # primary_only not with feet should be the minimal number of edges and points on reconstructed polygon
        # not primary only (including secondary edges) and with feet to reconstruc original polygon
        # will give most edges and polygon points
        # primary only will likely make the outermost an external edge display
        # not display

        my $with_feet = 1;
        my $primary_only = 0;

        do {
            die "RETREAD ON SEEN EDGE" if (defined $seen{''.$$e} && $seen{''.$$e} == 1);
            $seen{''.$$e} = 1;
            
            push @loop_edges, $e if (!$primary_only || $e->is_primary());
                        
            if ($with_feet) {
                if ($e->foot()) {
                    push @poly_verts, [$e->foot()->x(),$e->foot()->y()];
                }
            }
            else {
                if ($e->is_finite() && $e->vertex1() && $e->vertex1()->r() == 0) {
                    push @poly_verts, [$e->vertex1()->x(),$e->vertex1()->y()];
                }
            }
            
            $e = $e->next();
        } while($$e != $$iedge && $cnt++ < 10000);
        warn "is that a too-many-sided cell?? $cnt" if $cnt > 999;

        next if @loop_edges == 0;
        
        # recreates an original input polygon, with winding opposite
        # ma edge loop winding.
        $d.='M'.join('L',map {($_->[0]/$scale).','.($_->[1]/$scale)} reverse @poly_verts).'Z' if @poly_verts > 2;
        
        if (!$loop_edges[0]->vertex0()) {
            $d .= sprintf("M%.6f,%.6f",($loop_edges[0]->vertex1()->x()/$scale),($loop_edges[0]->vertex1()->y()/$scale));
        }
        else {
            $d .= sprintf("M%.6f,%.6f",($loop_edges[0]->vertex0()->x()/$scale),($loop_edges[0]->vertex0()->y()/$scale));
        }
        foreach my $edge (@loop_edges) {
            if ($edge->is_curved()) {
                my $q = $isma ? [$edge->Q()->x(),$edge->Q()->y()] : calc_Q($edge);
                $d .=  sprintf("Q%.6f,%.6f,%.6f,%.6f",$q->[0]/$scale, ($q->[1]/$scale),($edge->vertex1()->x()/$scale),($edge->vertex1()->y()/$scale));
            }
            else {
                if (!$edge->vertex1()) {
                    # then line to where we already are.
                    # TODO make this a fake infinite so we can see it better?
                    #      but would also have make fake infinite boundary 
                    #      polygon above for this to go inside of
                    $d .= sprintf("L%.6f,%.6f",($edge->vertex0()->x()/$scale),($edge->vertex0()->y()/$scale));
                    }
                else {
                    $d .= sprintf("L%.6f,%.6f",($edge->vertex1()->x()/$scale),($edge->vertex1()->y()/$scale));
                }
            }

        }

        $d .= 'Z';
        $svg .= $d.'" />'."\n";

    }


    return $svg;
}

1;