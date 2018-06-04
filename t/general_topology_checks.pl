
sub general_topology_checks {
    my $vd = shift;
    my $counts = shift;

    my @edges    = @{$vd->edges()};
    my @cells    = @{$vd->cells()};
    my @vertices = @{$vd->vertices()};
    
    my $is_ma = ref($vd)=~/medial/i        # when using builder directly
              ? 1 
              : (ref($cells[0])=~/medial/i # when using wrapped api
                 ? 1 : 0
                );
    my $vd_or_ma = $is_ma ? 'MA':'VD';
    
    my $count;

    $count = @edges;
    ok($count == $counts->{edges}, "$vd_or_ma: $counts->{edges} edges? ($count)");

    $count = scalar(grep $_->is_finite(), @edges);
    ok($count == $counts->{finite}, "$vd_or_ma: $counts->{finite} finite edges? ($count)");

    $count = scalar(grep $_->is_infinite(), @edges);
    ok($count == $counts->{infinite}, "$vd_or_ma: $counts->{infinite} infinite edges? ($count)");

    $count = scalar(grep $_->is_finite() && $_->is_secondary(), @edges);
    ok($count == $counts->{finite_secondary}, "$vd_or_ma: $counts->{finite_secondary} finite secondary edges? ($count)");

    $count = scalar(grep $_->is_infinite() && $_->is_secondary(), @edges);
    ok($count == $counts->{infinite_secondary}, "$vd_or_ma: $counts->{infinite_secondary} infinite secondary edges? ($count)");

    $count = scalar(grep $_->is_finite() && $_->is_primary(), @edges);
    ok($count == $counts->{finite_primary}, "$vd_or_ma: $counts->{finite_primary} finite primary edges? ($count)");

    $count = scalar(grep $_->is_infinite() && $_->is_primary(), @edges);
    ok($count == $counts->{infinite_primary}, "$vd_or_ma: $counts->{infinite_primary} infinite primary edges? ($count)");

    $count = scalar(grep $_->is_curved(), @edges);
    ok($count == $counts->{curved}, "$vd_or_ma: $counts->{curved} curved edges? ($count)");

    $count = @vertices;
    ok($count == $counts->{vertices}, "$vd_or_ma: $counts->{vertices} vertices? ($count)");

    $count = @cells;
    ok($count == $counts->{cells}, "$vd_or_ma: $counts->{cells} cells? ($count)");

    $count = grep $_->is_degenerate(), @cells;
    ok($count == $counts->{degenerate_cells}, "$vd_or_ma: $counts->{degenerate_cells} degenerate cells? ($count)");

    my $cell_fail = '';
    my $cell_edges_count_next = 0;
    my $cell_edges_count_prev = 0;
    for (my $i=0;$i < @cells; $i++) {
        next if $cells[$i]->is_degenerate();
        my $e_first = $cells[$i]->incident_edge();
        my $e = $e_first;
        my $limit = defined $counts->{_cell_loop_limit} ? $counts->{_cell_loop_limit} : 1000;
        do {
            $cell_edges_count_next++;
            $e = $e->next();
        } while ($$e != $$e_first && $limit-- > 0);
        if ( ! ($limit > 0) ) {
            $cell_fail .= "(Reached arbitrary loop limit while following edge->next() around cell $i.) ";        
        }

        $e = $e_first;
        $limit = defined $counts->{_cell_loop_limit} ? $counts->{_cell_loop_limit} : 1000;
        do {
            $cell_edges_count_prev++;
            $e = $e->prev();
        } while ($$e != $$e_first && $limit-- > 0);
        if ( ! ($limit > 0) ) {
            $cell_fail .= "(Reached arbitrary loop limit while following edge->prev() around cell $i.) ";        
        }

        if ($cell_fail) {last;}
    }

    ok( ! $cell_fail, "$vd_or_ma: Edge prev() and next() references create closed loops. $cell_fail");

    if ($is_ma) {
        ok(1,"(skipping cell loop edge count for medial axis)");
        ok(1,"(skipping cell loop edge count for medial axis)");
    }
    else {
        ok( $cell_edges_count_next == scalar(@edges), "$vd_or_ma: edge->next() loops edge totals ok. ($cell_edges_count_next/".scalar(@edges).")");
        ok( $cell_edges_count_prev == scalar(@edges), "$vd_or_ma: edge->prev() loops edge totals ok. ($cell_edges_count_prev/".scalar(@edges).")");
    }

    my $twin_fail = '';
    for (my $i=0;$i < @edges; $i++) {
        if (${$edges[$i]} != ${$edges[$i]->twin()->twin()}) {
            $twin_fail = "(The twin of edge $i did not refer back to it as it's twin.)";
            last;
        }
    }
    ok( ! $twin_fail, "$vd_or_ma: Edge twin() references are all reciprocal. $twin_fail");

    my $prev_next_fail = '';
    for (my $i=0;$i < @edges; $i++) {
        if (${$edges[$i]} != ${$edges[$i]->prev()->next()}) {
            $prev_next_fail = "(The prev() of edge $i did not refer back to it as next().)";
            last;
        }
    }
    ok( ! $prev_next_fail, "$vd_or_ma: Edge prev() and next() references are all reciprocal. $prev_next_fail");

    my $next_prev_fail = '';
    for (my $i=0;$i < @edges; $i++) {
        if (${$edges[$i]} != ${$edges[$i]->next()->prev()}) {
            $next_prev_fail = "(The next() of edge $i did not refer back to it as prev().)";
            last;
        }
    }
    ok( ! $next_prev_fail, "$vd_or_ma: Edge next() and prev() references are all reciprocal. $next_prev_fail");

}

1;
