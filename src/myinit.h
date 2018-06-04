#include <boost/polygon/voronoi_builder.hpp>
#include <boost/polygon/voronoi_diagram.hpp>
#include "medial_axis.hpp"

using namespace boost::polygon;

typedef voronoi_diagram<double>    voronoi_diagram_double;
typedef voronoi_vertex<double>     voronoi_vertex_double;
typedef voronoi_edge<double>       voronoi_edge_double;
typedef voronoi_cell<double>       voronoi_cell_double;

typedef medial_axis<double>        medial_axis_double;
typedef medial_axis_vertex<double> medial_axis_vertex_double;
typedef medial_axis_edge<double>   medial_axis_edge_double;
typedef medial_axis_cell<double>   medial_axis_cell_double;
typedef detail::point_2d<double>   medial_axis_foot;

// Use these mostly, so the same typemap code can get the perl
// package names from these typenames.
// One exception - in the .xsp code where you declare the 
// class for each, use the above types. Otherwise I think
// you get C function name conflicts and it doesn't work.
// For instance, in xspp code do like this:
//
//%name{Boost::Polygon::Voronoi::Edge} class voronoi_edge_double
// { ... }
//
// ... but everywhere else pretty much in the typemap stuff
// you use the names below.
typedef voronoi_diagram_double    Boost_Polygon_Voronoi_Diagram;
typedef voronoi_vertex_double     Boost_Polygon_Voronoi_Vertex;
typedef voronoi_edge_double       Boost_Polygon_Voronoi_Edge;
typedef voronoi_cell_double       Boost_Polygon_Voronoi_Cell;

typedef medial_axis_double        Boost_Polygon_MedialAxis;
typedef medial_axis_vertex_double Boost_Polygon_MedialAxis_Vertex;
typedef medial_axis_edge_double   Boost_Polygon_MedialAxis_Edge;
typedef medial_axis_cell_double   Boost_Polygon_MedialAxis_Cell;
typedef medial_axis_foot          Boost_Polygon_MedialAxis_Foot;

