// medial_axis.hpp header file
// derived from
// Boost.Polygon library voronoi_diagram.hpp header file
// which is
//          Copyright Andrii Sydorchuk 2010-2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org for updates, documentation, and revision history.

// Derivative work by Michael E. Sheldrake, Copyright 2013, distributed
// under the same terms as the license above.

// This is essentially boost/polygon/voronoi_diagram.hpp adapted to further 
// process the Voronoi diagram graph structure to make it represent the 
// medial axis of a polygon (with holes) represented by segment input.

#ifndef BOOST_POLYGON_MEDIAL_AXIS
#define BOOST_POLYGON_MEDIAL_AXIS

#include <vector>
#include <utility>
#include <cstdio>
#include <string>
#include <math.h>

#include "boost/lexical_cast.hpp"
#include "boost/polygon/detail/voronoi_ctypes.hpp"
#include "boost/polygon/detail/voronoi_structures.hpp"

#include "boost/polygon/voronoi_geometry_type.hpp"

#define to_str(n) boost::lexical_cast<std::string>( n )

namespace boost {
namespace polygon {

// Forward declarations.
template <typename T>
class medial_axis_edge;
 
// Represents Voronoi cell.
// Data members:
//   1) index of the source within the initial input set
//   2) pointer to the incident edge
//   3) mutable color member
// Cell may contain point or segment site inside.
template <typename T>
class medial_axis_cell {
 public:
  typedef T coordinate_type;
  typedef std::size_t color_type;
  typedef medial_axis_edge<coordinate_type> medial_axis_edge_type;
//  typedef std::size_t source_index_type;
  typedef size_t source_index_type;
  typedef SourceCategory source_category_type;

  medial_axis_cell(source_index_type source_index,
               source_category_type source_category) :
      source_index_(source_index),
      incident_edge_(NULL),
      color_(source_category) {}

  // Returns true if the cell contains point site, false else.
  bool contains_point() const {
    source_category_type source_category = this->source_category();
    return belongs(source_category, GEOMETRY_CATEGORY_POINT);
  }

  // Returns true if the cell contains segment site, false else.
  bool contains_segment() const {
    source_category_type source_category = this->source_category();
    return belongs(source_category, GEOMETRY_CATEGORY_SEGMENT);
  }

  source_index_type source_index() const {
    return source_index_;
  }

  source_category_type source_category() const {
    return static_cast<source_category_type>(color_ & SOURCE_CATEGORY_BITMASK);
  }

  // Degenerate cells don't have any incident edges.
  bool is_degenerate() const { return incident_edge_ == NULL; }

  medial_axis_edge_type* incident_edge() { return incident_edge_; }
  const medial_axis_edge_type* incident_edge() const { return incident_edge_; }
  void incident_edge(medial_axis_edge_type* e) { incident_edge_ = e; }

  color_type color() const { return color_ >> BITS_SHIFT; }
  void color(color_type color) const {
    color_ &= BITS_MASK;
    color_ |= color << BITS_SHIFT;
  }

 private:
  // 5 color bits are reserved.
  enum Bits {
    BITS_SHIFT = 0x5,
    BITS_MASK = 0x1F
  };

  source_index_type source_index_;
  medial_axis_edge_type* incident_edge_;
  mutable color_type color_;
};

// Represents Voronoi vertex.
// Data members:
//   1) vertex coordinates
//   2) radius of a maximal inscribed circle to the polygon at the vertex
//   3) pointer to the incident edge
//   4) mutable color member
template <typename T>
class medial_axis_vertex {
 public:
  typedef T coordinate_type;
  typedef std::size_t color_type;
  typedef medial_axis_edge<coordinate_type> medial_axis_edge_type;

  medial_axis_vertex(const coordinate_type& x, const coordinate_type& y,
                     const coordinate_type& r=0) :
      x_(x),
      y_(y),
      r_(r),
      incident_edge_(NULL),
      color_(0) {}

  const coordinate_type& x() const { return x_; }
  const coordinate_type& y() const { return y_; }
  const coordinate_type& r() const { return r_; }

  bool is_degenerate() const { return incident_edge_ == NULL; }

  medial_axis_edge_type* incident_edge() { return incident_edge_; }
  const medial_axis_edge_type* incident_edge() const { return incident_edge_; }
  void incident_edge(medial_axis_edge_type* e) { incident_edge_ = e; }

  color_type color() const { return color_ >> BITS_SHIFT; }
  void color(color_type color) const {
    color_ &= BITS_MASK;
    color_ |= color << BITS_SHIFT;
  }

 private:
  // 5 color bits are reserved.
  enum Bits {
    BITS_SHIFT = 0x5,
    BITS_MASK = 0x1F
  };

  coordinate_type x_;
  coordinate_type y_;
  coordinate_type r_;
  medial_axis_edge_type* incident_edge_;
  mutable color_type color_;
};

// Half-edge data structure. Represents Voronoi edge.
// Data members:
//   1) pointer to the corresponding cell
//   2) pointer to the vertex that is the starting
//      point of the half-edge
//   3) pointer to the twin edge
//   4) pointer to the CCW next edge
//   5) pointer to the CCW prev edge
//   6) boolean indicating whether foot coordinates have been set
//   7) mutable color member
template <typename T>
class medial_axis_edge {
 public:
  typedef T coordinate_type;
  typedef medial_axis_cell<coordinate_type> medial_axis_cell_type;
  typedef medial_axis_vertex<coordinate_type> medial_axis_vertex_type;
  typedef medial_axis_edge<coordinate_type> medial_axis_edge_type;
  typedef std::size_t color_type;

  medial_axis_edge(bool is_linear, bool is_primary) :
      cell_(NULL),
      vertex_(NULL),
      twin_(NULL),
      next_(NULL),
      prev_(NULL),
      next_vd_(NULL),
      prev_vd_(NULL),
      footset_(false),
      Qset_(false),
      theta_(4), // > pi indicates undefined
      phi_(4),
      color_(0) {
    if (is_linear)
      color_ |= BIT_IS_LINEAR;
    if (is_primary)
      color_ |= BIT_IS_PRIMARY;
  }

  medial_axis_cell_type* cell() { return cell_; }
  const medial_axis_cell_type* cell() const { return cell_; }
  void cell(medial_axis_cell_type* c) { cell_ = c; }

  medial_axis_vertex_type* vertex0() { return vertex_; }
  const medial_axis_vertex_type* vertex0() const { return vertex_; }
  void vertex0(medial_axis_vertex_type* v) { vertex_ = v; }

  medial_axis_vertex_type* vertex1() { return twin_->vertex0(); }
  const medial_axis_vertex_type* vertex1() const { return twin_->vertex0(); }

  medial_axis_edge_type* twin() { return twin_; }
  const medial_axis_edge_type* twin() const { return twin_; }
  void twin(medial_axis_edge_type* e) { twin_ = e; }

  medial_axis_edge_type* next() { return next_; }
  const medial_axis_edge_type* next() const { return next_; }
  void next(medial_axis_edge_type* e) { next_ = e; }

  medial_axis_edge_type* next_vd() { return next_vd_ == NULL ? next_ : next_vd_; }
  const medial_axis_edge_type* next_vd() const { return next_vd_ == NULL ? next_ : next_vd_; }
  void next_vd(medial_axis_edge_type* e) { next_vd_ = e; }

  medial_axis_edge_type* prev() { return prev_; }
  const medial_axis_edge_type* prev() const { return prev_; }
  void prev(medial_axis_edge_type* e) { prev_ = e; }

  medial_axis_edge_type* prev_vd() { return prev_vd_ == NULL ? prev_ : prev_vd_; }
  const medial_axis_edge_type* prev_vd() const { return prev_vd_ == NULL ? prev_ : prev_vd_; }
  void prev_vd(medial_axis_edge_type* e) { prev_vd_ = e; }

  // Returns a pointer to the rotation next edge
  // over the starting point of the half-edge.
  medial_axis_edge_type* rot_next() { return prev_->twin(); }
  const medial_axis_edge_type* rot_next() const { return prev_->twin(); }

  medial_axis_edge_type* rot_next_vd() { return prev_vd_ == NULL ? prev_->twin() : prev_vd_->twin(); }
  const medial_axis_edge_type* rot_next_vd() const { return prev_vd_ == NULL ? prev_->twin() : prev_vd_->twin(); }

  // Returns a pointer to the rotation prev edge
  // over the starting point of the half-edge.
  medial_axis_edge_type* rot_prev() { return twin_->next(); }
  const medial_axis_edge_type* rot_prev() const { return twin_->next(); }

  medial_axis_edge_type* rot_prev_vd() { return twin_->next_vd(); }
  const medial_axis_edge_type* rot_prev_vd() const { return twin_->next_vd(); }

  // Returns true if the edge is finite (segment, parabolic arc).
  // Returns false if the edge is infinite (ray, line).
  bool is_finite() const { return vertex0() && vertex1(); }

  // Returns true if the edge is infinite (ray, line).
  // Returns false if the edge is finite (segment, parabolic arc).
  bool is_infinite() const { return !vertex0() || !vertex1(); }

  // Returns true if the edge is linear (segment, ray, line).
  // Returns false if the edge is curved (parabolic arc).
  bool is_linear() const {
    return (color_ & BIT_IS_LINEAR) ? true : false;
  }

  // Returns true if the edge is curved (parabolic arc).
  // Returns false if the edge is linear (segment, ray, line).
  bool is_curved() const {
    return (color_ & BIT_IS_LINEAR) ? false : true;
  }

  // Returns false if edge goes through the endpoint of the segment.
  // Returns true else.
  bool is_primary() const {
    return (color_ & BIT_IS_PRIMARY) ? true : false;
  }

  // Returns true if edge goes through the endpoint of the segment.
  // Returns false else.
  bool is_secondary() const {
    return (color_ & BIT_IS_PRIMARY) ? false : true;
  }

  // Returns true if edge is within a closed polygon formed by input segments
  // Returns false else.
  bool is_internal() const {
    return (color_ & BIT_IS_INTERNAL) ? true : false;
  }
  void is_internal(bool internal) const {
    if (internal) {
      color_ |= BIT_IS_INTERNAL;
      color_ &= (BITS_MASK ^ BIT_IS_EXTERNAL);
    }
  }

  // Returns true if edge is outside of all closed polygons
  // formed by input segments, or within a hole polygon
  // NEEDS TESTING - untested on nested polygons in holes
  // Returns false else.
  bool is_external() const {
    return (color_ & BIT_IS_EXTERNAL) ? true : false;
  }
  void is_external(bool external) const {
    if (external) {
      color_ |= BIT_IS_EXTERNAL;
      color_ &= (BITS_MASK ^ BIT_IS_INTERNAL);
    }
  }

  color_type color() const { return color_ >> BITS_SHIFT; }
  void color(color_type color) const {
    color_ &= BITS_MASK;
    color_ |= color << BITS_SHIFT;
  }
  
  // foot: where radius from vertex0 touches source segment at a 90 degree angle
  const detail::point_2d<coordinate_type>* foot() const { 
    if (!footset_) {return NULL;}
    return &foot_;
  }
  detail::point_2d<coordinate_type>* foot() { 
    if (!footset_) {return NULL;}
    return &foot_;
  }
  void foot(coordinate_type x, coordinate_type y) {
    footset_ = true;
    foot_.x(x);
    foot_.y(y);
  }
  // Q: Quadratic Bezier control point for curved edges
  const detail::point_2d<coordinate_type>* Q() const { 
    if (!Qset_) {return NULL;}
    return &Q_;
  }
  detail::point_2d<coordinate_type>* Q() { 
    if (!Qset_) {return NULL;}
    return &Q_;
  }
  void Q(coordinate_type x, coordinate_type y) {
    Qset_ = true;
    Q_.x(x);
    Q_.y(y);
  }

  // theta
  const double theta() const { return theta_; }
  double theta() { return theta_; }
  void theta(double theta) { theta_ = theta; }

  // phi
  const double phi() const { return phi_; }
  double phi() { return phi_; }
  void phi(double phi) { phi_ = phi; }


 private:
  // 5 color bits are reserved.
  enum Bits {
    BIT_IS_LINEAR = 0x1,  // linear is opposite to curved
    BIT_IS_PRIMARY = 0x2,  // primary is opposite to secondary
    BIT_IS_INTERNAL = 0x4,  // internal is opposite to external
    BIT_IS_EXTERNAL = 0x8,  // internal is opposite to external

    BITS_SHIFT = 0x5,
    BITS_MASK = 0x1F
  };

  medial_axis_cell_type* cell_;
  medial_axis_vertex_type* vertex_;
  medial_axis_edge_type* twin_;
  medial_axis_edge_type* next_;
  medial_axis_edge_type* prev_;
  medial_axis_edge_type* next_vd_;
  medial_axis_edge_type* prev_vd_;
  mutable color_type color_;
  mutable detail::point_2d<coordinate_type> foot_;
  bool footset_;
  mutable detail::point_2d<coordinate_type> Q_;
  bool Qset_;
  double theta_;
  double phi_;
  //mutable detail::point_2d<default_voronoi_builder::int_type> p1_;

};

template <typename T>
struct medial_axis_traits {
  typedef T coordinate_type;
  typedef medial_axis_cell<coordinate_type> cell_type;
  typedef medial_axis_vertex<coordinate_type> vertex_type;
  typedef medial_axis_edge<coordinate_type> edge_type;
  typedef class {
   public:
    enum { ULPS = 128 };
    bool operator()(const vertex_type& v1, const vertex_type& v2) const {
      return (ulp_cmp(v1.x(), v2.x(), ULPS) ==
              detail::ulp_comparison<T>::EQUAL) &&
             (ulp_cmp(v1.y(), v2.y(), ULPS) ==
              detail::ulp_comparison<T>::EQUAL);
    }
   private:
    typename detail::ulp_comparison<T> ulp_cmp;
  } vertex_equality_predicate_type;
};

// Voronoi output data structure.
// CCW ordering is used on the faces perimeter and around the vertices.
template <typename T, typename TRAITS = medial_axis_traits<T> >
class medial_axis {
 public:
  typedef typename TRAITS::coordinate_type coordinate_type;
  typedef typename TRAITS::cell_type cell_type;
  typedef typename TRAITS::vertex_type vertex_type;
  typedef typename TRAITS::edge_type edge_type;

  typedef std::vector<cell_type> cell_container_type;
  typedef typename cell_container_type::const_iterator const_cell_iterator;

  typedef std::vector<vertex_type> vertex_container_type;
  typedef typename vertex_container_type::const_iterator const_vertex_iterator;

  typedef std::vector<edge_type> edge_container_type;
  typedef typename edge_container_type::const_iterator const_edge_iterator;

  medial_axis() {}

  void clear() {
    cells_.clear();
    vertices_.clear();
    edges_.clear();
  }

  const cell_container_type& cells() const {
    return cells_;
  }

  const vertex_container_type& vertices() const {
    return vertices_;
  }

  const edge_container_type& edges() const {
    return edges_;
  }

  std::string& event_log() const { // available for debugging
    return event_log_;
  }

  std::size_t num_cells() const {
    return cells_.size();
  }

  std::size_t num_edges() const {
    return edges_.size();
  }

  std::size_t num_vertices() const {
    return vertices_.size();
  }

  void _reserve(int num_sites) {
    cells_.reserve(num_sites);
    vertices_.reserve(num_sites << 1);
    edges_.reserve((num_sites << 2) + (num_sites << 1));
  }

  template <typename CT>
  void _process_single_site(const detail::site_event<CT>& site) {
    cells_.push_back(cell_type(site.initial_index(), site.source_category()));
  }

  // Insert a new half-edge into the output data structure.
  // Takes as input left and right sites that form a new bisector.
  // Returns a pair of pointers to a new half-edges.
  template <typename CT>
  std::pair<void*, void*> _insert_new_edge(
      const detail::site_event<CT>& site1,
      const detail::site_event<CT>& site2) {
    //printf("site event A \n");
    // Get sites' indexes.
    int site_index1 = site1.sorted_index();
    int site_index2 = site2.sorted_index();

    bool is_linear = is_linear_edge(site1, site2);
    bool is_primary = is_primary_edge(site1, site2);

    // Create a new half-edge that belongs to the first site.
    edges_.push_back(edge_type(is_linear, is_primary));
    edge_type& edge1 = edges_.back();

    // Create a new half-edge that belongs to the second site.
    edges_.push_back(edge_type(is_linear, is_primary));
    edge_type& edge2 = edges_.back();

    // Add the initial cell during the first edge insertion.
    if (cells_.empty()) {
      cells_.push_back(cell_type(
          site1.initial_index(), site1.source_category()));
    }

    // The second site represents a new site during site event
    // processing. Add a new cell to the cell records.
    cells_.push_back(cell_type(
        site2.initial_index(), site2.source_category()));

    // Set up pointers to cells.
    edge1.cell(&cells_[site_index1]);
    edge2.cell(&cells_[site_index2]);

    // Set up twin pointers.
    edge1.twin(&edge2);
    edge2.twin(&edge1);
    

    /* keep this in sync with same below in other _insert_new_edge() */
    // aliases so we can use the same case code as much as possible
    edge_type & new_edge1 = edge1;
    edge_type & new_edge2 = edge2;

    bool case_report = false;
    
    if (case_report) printf("A: ");

    if (new_edge2.is_secondary()) {
      if (case_report) printf("CASE 0 ");
      // For secondaries, one site is a point at end of the other site,
      // which is a segment. We want the foot to be the point.
      const detail::site_event<CT>& point_site = site1.is_point() ? site1 : site2;
      new_edge2.foot(point_site.x(),point_site.y()); // do 2 as case 0.1 to correspond to similar circle event case numbering
      if (case_report) printf("F0.1[%llu, %d, %d] ",(long long unsigned) &new_edge2, point_site.x(),point_site.y());
      new_edge1.foot(point_site.x(),point_site.y());
      if (case_report) printf("F0.2[%llu, %d, %d] ",(long long unsigned) &new_edge1, point_site.x(),point_site.y());
    } 
    else { // is primary
      if (new_edge2.is_curved()) {
        if (new_edge2.cell()->contains_point()) {
          if (case_report) printf("CASE 1 ");
          new_edge2.foot(site2.x(),site2.y());
          if (case_report) printf("F1.1[%llu, %d, %d] ",(long long unsigned) &new_edge2, site2.x(),site2.y());
        } 
        else { // contains segment
          if (case_report) printf("CASE 2 ");
          new_edge1.foot(site1.x(),site1.y());
          if (case_report) printf("F1.1[%llu, %d, %d] ",(long long unsigned) &new_edge1, site1.x(),site1.y());
        }
      }
      else { // is linear
        if (new_edge1.cell()->contains_point() && new_edge2.cell()->contains_point()) {
          if (case_report) printf("CASE A ");
          new_edge1.foot(site1.x(),site1.y()); // CASE A.1
          if (case_report) printf("FA.1[%llu, %d, %d] ",(long long unsigned) &new_edge1, site1.x(),site1.y());
          new_edge2.foot(site2.x(),site2.y()); // CASE A.2
          if (case_report) printf("FA.2[%llu, %d, %d] ",(long long unsigned) &new_edge2, site2.x(),site2.y());
        }
        else if (new_edge2.prev() != NULL && new_edge1.next() != NULL 
            && new_edge2.prev()->is_linear() && new_edge1.next()->is_curved()) {
          if (case_report) printf("CASE 3 ");
        }
        else if (new_edge2.prev() != NULL && new_edge1.next() != NULL
            && new_edge2.prev()->is_curved() && new_edge1.next()->is_linear()) {
          if (case_report) printf("CASE 4 ");
        }
        else if (new_edge2.prev() != NULL && new_edge1.next() != NULL
            && new_edge2.prev()->is_curved() && new_edge1.next()->is_curved()) {
          if (case_report) printf("CASE 5 ");
        }
        else if (new_edge2.prev() != NULL && new_edge1.next() != NULL
            && new_edge2.prev()->is_linear() && new_edge1.next()->is_linear()) {
          if (new_edge2.vertex0() && new_edge2.vertex0()->r() == 0) { 
            // coming out of corner
            if (case_report) printf("CASE 6 ");
          }
          else if (new_edge2.prev()->vertex0() && new_edge2.prev()->vertex0()->r() == 0) { 
            // prev came out of corner
            if (case_report) printf("CASE 7 ");
          }
          else if (new_edge2.prev()->cell()->contains_point() && new_edge1.next()->cell()->contains_point()) {
            // this might need special handling. it's prob most like is_curved() && is_curved() case
            // maybe it's the same and you can combine them
            if (case_report) printf("CASE 8 ");
          }
          else { // no corner immediately behind
            if (case_report) printf("CASE 9 ");
          }
        }
        else {
          if (case_report) printf("CASE ?? linear "); // ideally you never get here
        }
      }
    }

    if (case_report) printf("\n          ");
    if (case_report) printf("%s e1:%lld e2:%lld ",edge1.is_primary()?"P":"S",(long long unsigned) &edge1, (long long unsigned) &edge2);
    if (case_report) printf("\n          site1: [%d, %d] ",site1.x(),site1.y());
    if (case_report) {if (site1.is_segment()) {printf("[%d, %d] ",site1.x1(),site1.y1());}}
    if (case_report) printf("site2: [%d, %d] ",site2.x(),site2.y());
    if (case_report) {if (site2.is_segment()) {printf("[%d, %d] ",site2.x1(),site2.y1());}}
    if (case_report) printf("\n");



    // Return a pointer to the new half-edge.
    return std::make_pair(&edge1, &edge2);
  }

  // Insert a new half-edge into the output data structure with the
  // start at the point where two previously added half-edges intersect.
  // Takes as input two sites that create a new bisector, circle event
  // that corresponds to the intersection point of the two old half-edges,
  // pointers to those half-edges. Half-edges' direction goes out of the
  // new Voronoi vertex point. Returns a pair of pointers to a new half-edges.
  template <typename CT1, typename CT2>
  std::pair<void*, void*> _insert_new_edge(
      const detail::site_event<CT1>& site1,
      const detail::site_event<CT1>& site3,
      const detail::circle_event<CT2>& circle,
      void* data12, void* data23) {
    edge_type* edge12 = static_cast<edge_type*>(data12);
    edge_type* edge23 = static_cast<edge_type*>(data23);
    //printf("circle event\n");

    // Add a new Voronoi vertex.
    vertices_.push_back(vertex_type(circle.x(), circle.y(), 
                                    circle.lower_x() - circle.x()));
    vertex_type& new_vertex = vertices_.back();
    
    // Update vertex pointers of the old edges.
    edge12->vertex0(&new_vertex);
    edge23->vertex0(&new_vertex);

    bool is_linear = is_linear_edge(site1, site3);
    bool is_primary = is_primary_edge(site1, site3);

    // Add a new half-edge.
    edges_.push_back(edge_type(is_linear, is_primary));
    edge_type& new_edge1 = edges_.back();
    new_edge1.cell(&cells_[site1.sorted_index()]);

    // Add a new half-edge.
    edges_.push_back(edge_type(is_linear, is_primary));
    edge_type& new_edge2 = edges_.back();
    new_edge2.cell(&cells_[site3.sorted_index()]);

    // Update twin pointers.
    new_edge1.twin(&new_edge2);
    new_edge2.twin(&new_edge1);

    // Update vertex pointer.
    new_edge2.vertex0(&new_vertex);

    // Update Voronoi prev/next pointers.
    edge12->prev(&new_edge1);
    new_edge1.next(edge12);
    edge12->twin()->next(edge23);
    edge23->prev(edge12->twin());
    edge23->twin()->next(&new_edge2);
    new_edge2.prev(edge23->twin());

    //set foot
    // It seems possible that we can do all foot-finding in this event processing
    // (here in the circle event, and in the other site even code above).
    // But we're not sure we've completely understood or diagramed all the cases
    // of different edges and vertices being available during these events.
    // 
    // We've captured most cases though. Missing feet are now uncommon, and usually
    // due to degenerate input polygons. Think there are some legit edge cases
    // though that we still need to work out.
    // 
    
    bool case_report = false;

    if (case_report) printf("B: nv[%f %f]\n   ",new_vertex.x(),new_vertex.y());

    if (new_vertex.r() == 0) {
      if (case_report) printf("CASE 0 ");
      new_edge2.foot(new_vertex.x(),new_vertex.y()); // whether primary or secondary, edge coming out of v.r==0 has foot there
      if (case_report) printf("F0.1[%llu, %f, %f] ",(long long unsigned) &new_edge2,new_vertex.x(),new_vertex.y());
      if (new_edge2.is_secondary()) { // if secondary, edge coming into v.r==0 also has foot there
        new_edge1.foot(new_vertex.x(),new_vertex.y());
        if (case_report) printf("F0.2[%llu, %f, %f] ",(long long unsigned) &new_edge1,new_vertex.x(),new_vertex.y());
      }
      if (edge12->is_secondary()) {
        edge12->twin()->foot(new_vertex.x(),new_vertex.y());
        edge12->foot(new_vertex.x(),new_vertex.y());
        if (case_report) printf("F0.3[%llu, %f, %f] ",(long long unsigned) edge12->twin(), new_vertex.x(),new_vertex.y());
      }
      if (edge23->is_secondary()) {
        edge23->twin()->foot(new_vertex.x(),new_vertex.y());
        edge23->foot(new_vertex.x(),new_vertex.y());
        if (case_report) printf("F0.4[%llu, %f, %f] ",(long long unsigned) edge23->twin(), new_vertex.x(),new_vertex.y());
      }      
    } 
    else if (new_edge2.is_primary()) { // is primary
      if (new_edge2.is_curved()) {
        if (new_edge2.cell()->contains_point()) {
          if (case_report) printf("CASE 1 \n");
          new_edge2.foot(site3.x(),site3.y()); // CASE 1.1 // point site is focus of curved new_edge2
          if (case_report) printf("F1.1[%llu, %d, %d] ",(long long unsigned) &new_edge2, site3.x(),site3.y());
        }
        else { // contains segment
          if (case_report) printf("CASE 2 ");
          new_edge1.foot(site1.x(), site1.y()); // CASE 2.1 // point site is focus of curved new_edge1 
          if (case_report) printf("F2.1[%llu, %d, %d] ",(long long unsigned) &new_edge1, site1.x(), site1.y());
          
          if (edge23->is_finite() && edge23->is_linear()) { //edge23 was infinite in open box (3 sides) case... what whould happen there? 
            // CASE 2.2
            double x = edge23->vertex0()->x();
            double y = edge23->vertex0()->y();

            makefoot(x, y, site3.x0(), site3.y0(), site3.x1(), site3.y1());
            new_edge2.foot(x,y);
            if (case_report) printf("F2.2[%llu, %f, %f] ",(long long unsigned) &new_edge2, x,y);

            reflect(x, y, edge23->vertex0()->x(),edge23->vertex0()->y(),
                          edge23->vertex1()->x(),edge23->vertex1()->y());
            edge23->foot(x,y);
            if (case_report) printf("F2.3[%llu, %f, %f] ",(long long unsigned) edge23, x,y);
          }
          if (edge23->is_finite() && edge23->is_curved() && edge23->next() && edge23->next()->is_linear() && edge23->next()->foot()) {
            // CASE 2.3
            edge23->foot(edge23->next()->foot()->x(),edge23->next()->foot()->y());
            if (case_report) printf("F2.4[%llu, %f, %f] ",(long long unsigned) edge23, edge23->next()->foot()->x(),edge23->next()->foot()->y());
          }
        }
      }
      else { // is linear
        if (new_edge1.cell()->contains_point() && new_edge2.cell()->contains_point()) {
          if (case_report) printf("CASE A ");
          new_edge1.foot(site1.x(),site1.y()); // CASE A.1
          if (case_report) printf("FA.1[%llu, %d, %d] ",(long long unsigned) &new_edge1, site1.x(),site1.y());
          new_edge2.foot(site3.x(),site3.y()); // CASE A.2
          if (case_report) printf("FA.2[%llu, %d, %d] ",(long long unsigned) &new_edge2, site3.x(),site3.y());
          
          if (new_edge2.is_linear()
              // turns out not to be general case safe
              // but is legit for this case, always we hope
              && edge12->is_curved()
              && edge23->is_curved()
             ) {
            // There is a foot straight back behind new_edge2
            // at distance == radius, for edge23.
            // "Straight behind" slope is negative recipricol of
            // slope of line between site1 and site3.
            double x  = site1.x();
            double y  = site1.y();
            // negative recipricol slope to get perpendicular angle
            double a = atan2(-(site3.x() - site1.x()), site3.y() - site1.y());
            rotate_2d(x, y, a, new_vertex.x(), new_vertex.y());
            edge23->foot(x,y);
            if (case_report) printf("FA.3[%llu, %f, %f] ",(long long unsigned) edge23, x, y);
          }
          
          // that the following case is needed for some feet but overwrites other
          // existing feet with the wrong thing suggests we should try to find 
          // the unique case that handles the needed ones,  and doesn't see any
          // of the ones already handled by others
          if ( ! new_edge1.next()->twin()->next()->foot()
              && new_edge1.next()->twin()->next()->is_primary()) {
            new_edge1.next()->twin()->next()->foot(site1.x(),site1.y()); // CASE A.3 // coming out of not-curve curve, similar to case 2.2 and 3.1
            if (case_report) printf("FA.4[%llu, %d, %d] ",(long long unsigned) new_edge1.next()->twin()->next(), site1.x(),site1.y());
          }
        }
        else if (new_edge2.prev() != NULL && new_edge1.next() != NULL 
            && new_edge2.prev()->is_linear() && new_edge1.next()->is_curved()) {
          if (case_report) printf("CASE 3 ");
          if (new_edge2.is_primary()) { // CASE 3.1
                                        // similar to CASE 2.2, setting foot on edge coming out of curve
                                        // new_edge2->prev() was a 2ndary, coming out of curve around the
                                        // start end of the site3 segment. Hopefully it's always the start
                                        // end in this case. otherwise can probably pick end by looking at
                                        // site's is_inverse() setting
            new_edge2.foot(site3.x0(),site3.y0());
            if (case_report) printf("F3.1[%llu, %d, %d] ",(long long unsigned) &new_edge2, site3.x0(),site3.y0());
            }
        }
        else if (new_edge2.prev() != NULL && new_edge1.next() != NULL
            && new_edge2.prev()->is_curved() && new_edge1.next()->is_linear()) {
          if (case_report) printf("CASE 4 ");

          // CASE 4.1
          double x = new_edge2.vertex0()->x();
          double y = new_edge2.vertex0()->y();
          makefoot(x, y, site3.x0(), site3.y0(), site3.x1(), site3.y1());
          new_edge2.foot(x,y);
          if (case_report) printf("F4.1[%llu, %f, %f] ",(long long unsigned) &new_edge2, x,y);

        }
        else if (new_edge2.prev() != NULL && new_edge1.next() != NULL
            && new_edge2.prev()->is_curved() && new_edge1.next()->is_curved()) {
          if (case_report) printf("CASE 5 ");
        }
        else if (new_edge2.prev() != NULL && new_edge1.next() != NULL
            && new_edge2.prev()->is_linear() && new_edge1.next()->is_linear()) {
          if (new_edge2.vertex0() && new_edge2.vertex0()->r() == 0) { 
            // coming out of corner
            if (case_report) printf("CASE 6 ");
            new_edge2.foot(site3.x0(),site3.y0()); // CASE 6.1
                                                   // based on one example so far that came up when
                                                   // I commented out some missing-foot-fill-in code in _build()
            if (case_report) printf("F6.1[%llu, %d, %d] ",(long long unsigned) &new_edge2, site3.x0(),site3.y0());
          }
          else if (new_edge2.prev()->vertex0() && new_edge2.prev()->vertex0()->r() == 0) { 
            // prev came out of corner
            if (case_report) printf("CASE 7 "); // was thinking internal, but saw this on a
          }
          else if (new_edge2.prev()->cell()->contains_point() && new_edge1.next()->cell()->contains_point()) {
            // this might need special handling. it's prob most like is_curved() && is_curved() case
            // maybe it's the same and you can combine them
            if (case_report) printf("CASE 8 ");
          }
          else { // no corner immediately behind
            if (case_report) printf("CASE 9 ");
          }
        }
        else {
          if (case_report) printf("CASE ?? linear "); // ideally you never get here
        }
      }
    }

    if (case_report) printf("\n          ");
    if (case_report) printf("%s ne1:%lld ne2:%lld e12:%lld e23:%lld ",
                            new_edge1.is_primary()?"P":"S",
                            (long long unsigned) &new_edge1,
                            (long long unsigned) &new_edge2,
                            (long long unsigned) edge12,
                            (long long unsigned) edge23);
    if (case_report) printf("\n          site1: [%d, %d] ",site1.x(),site1.y());
    if (case_report) {if (site1.is_segment()) {printf("[%d, %d] ",site1.x1(),site1.y1());}}
    if (case_report) printf("site3: [%d, %d] ",site3.x(),site3.y());
    if (case_report) {if (site3.is_segment()) {printf("[%d, %d] ",site3.x1(),site3.y1());}}
    if (case_report) printf("\n");



    // Below are foot-finding and setting cases initially worked out by trial and error.
    // These should all be obsolete. The cases worked out above should find all feet.
    // But until the more methodical approach to cases above, we're keeping the following
    // code arround to possibly help with debugging. (For instance, if the above cases fail
    // to find a foot, but one of the older ad-hoc cases below finds it, that's a clue.)

    /////////////////////
    // START OLD CASES //
    /////////////////////

    bool oldcases = false; // whether to print debug info for old foot-finding cases

    // for edges going into corners

    // edge12 into corner

    // note the cast from float to int for these ==s : LHS is an int type, RHS is float
    // nope: that didn't fix it, and that's what we want to avoid anyway
    if (edge12->vertex1()
        //&& (  (site1.point1().x() == (coordinate_type) edge12->vertex1()->x() && site1.point1().y() == (coordinate_type) edge12->vertex1()->y())
        //   || (site1.point0().x() == (coordinate_type) edge12->vertex1()->x() && site1.point0().y() == (coordinate_type) edge12->vertex1()->y())
        //   )
// If this really works, do same/similar for next if for edge23
// ... catches more cases than the target case, putting some feet not on polygon segments
//     but see if something similarly topological can do the job here
        && edge12->is_primary() 
        && edge12->next() && !edge12->next()->is_primary()
        && (  edge12->next()->cell()->contains_point()
           || edge12->next()->twin()->cell()->contains_point()
           )
        && belongs(site1.source_category(), GEOMETRY_CATEGORY_SEGMENT)
       ) {
      
      double x0 = site1.point0().x();
      double y0 = site1.point0().y();
      double x1 = site1.point1().x();
      double y1 = site1.point1().y();
      double x = new_vertex.x();
      double y = new_vertex.y();
      makefoot(x, y, x0, y0, x1, y1);
      edge12->foot(x, y);

      if (oldcases) printf("OC  1   [%llu] [%f,%f]\n",(long long unsigned) edge12,x,y);
    }

    // edge23 into corner
    if (edge23->vertex1()
        //&& 
        //(  (edge23->vertex1()->x() == site3.point1().x() && edge23->vertex1()->y() == site3.point1().y())
        //|| (edge23->vertex1()->x() == site3.point0().x() && edge23->vertex1()->y() == site3.point0().y())
        //)
        //&& edge23->is_primary() && edge23->next() && !edge23->next()->is_primary()

        && edge23->is_primary() 
        && edge23->next() && !edge23->next()->is_primary()
        && (
           edge23->next()->cell()->contains_point()
           || 
           edge23->next()->twin()->cell()->contains_point()
           )
       && belongs(site3.source_category(), GEOMETRY_CATEGORY_SEGMENT)) {
      
      double x = new_vertex.x();
      double y = new_vertex.y();
      makefoot(x, y, site3.point0().x(), site3.point0().y(),
                     site3.point1().x(), site3.point1().y());

      if (oldcases) printf("OC  2   [%llu] [%f,%f]\n",(long long unsigned) &new_edge2,x,y);
      if (oldcases && new_edge2.foot()) printf("  (foot already set, ignoring)\n");
      if (!new_edge2.foot()) new_edge2.foot(x, y);

    }

    // maybe
    if (edge12->is_primary()
        && edge12->vertex1()
        && belongs(site1.source_category(), GEOMETRY_CATEGORY_POINT)
        && (edge12->vertex1()->x() == site1.point0().x() && edge12->vertex1()->y() == site1.point0().y())
        ) {
      edge12->twin()->foot(edge12->vertex1()->x(), edge12->vertex1()->y());

      if (oldcases) printf("OC  2.5 [%llu] [%f,%f]\n",(long long unsigned) edge12->twin(),edge12->vertex1()->x(), edge12->vertex1()->y());

    }
    // maybe
    if (edge23->is_primary()
        && edge23->vertex1()
        && belongs(site3.source_category(), GEOMETRY_CATEGORY_POINT)
        && (edge23->vertex1()->x() == site3.point0().x() && edge23->vertex1()->y() == site3.point0().y())
        ) {
      edge23->twin()->foot(edge23->vertex1()->x(), edge23->vertex1()->y());

      if (oldcases) printf("OC  3   [%llu] [%f,%f]\n",(long long unsigned) edge23,edge23->vertex1()->x(), edge23->vertex1()->y());
    }

    // special case derived from visual debug
    if (   belongs(site1.source_category(), GEOMETRY_CATEGORY_POINT)
        && belongs(site3.source_category(), GEOMETRY_CATEGORY_SEGMENT)
       ) {
      if (
         (    (site3.point0().x() == site1.point0().x() && site3.point0().y() == site1.point0().y())
           || (site3.point1().x() == site1.point0().x() && site3.point1().y() == site1.point0().y())
         )
         && edge23->is_primary()
         && (edge23->vertex0()->x() == site1.point0().x() && edge23->vertex0()->y() == site1.point0().y())
         ) {
        
        // visually, looked like this was already there
        // but maybe that was wrongly associated with another edge?
        // so set/reset to see...
        // yeah this wasn't set, even though there appeared to be a foot there
        // so that was probably a wrong foot from something else
        // and this is maybe right
        edge23->foot(site1.point0().x(), site1.point0().y());

        if (oldcases) printf("OC  4.1 [%llu] [%f,%f]\n",(long long unsigned) edge23,(double) site1.point0().x(), (double)site1.point0().y());

        double x0 = site3.point0().x();
        double y0 = site3.point0().y();
        double x1 = site3.point1().x();
        double y1 = site3.point1().y();
        double x = edge23->vertex1()->x();
        double y = edge23->vertex1()->y();
        makefoot(x, y, x0, y0, x1, y1);
        edge23->twin()->foot(x, y);

        if (oldcases) printf("OC  4.2 [%llu] [%f,%f]\n",(long long unsigned) edge23->twin(),x,y);

        reflect(x, y, edge23->vertex0()->x(), edge23->vertex0()->y(),
                      edge23->vertex1()->x(), edge23->vertex1()->y());
        edge23->next()->foot(x, y);

        if (oldcases) printf("OC  4.3 [%llu] [%f,%f]\n",(long long unsigned) edge23->next(),x,y);

      }
    }

    // for straight edges
    if (edge23->vertex1() && !edge23->twin()->foot() 
        && edge23->is_linear() 
        //&& edge23->twin()->is_primary() 
        && belongs(site3.source_category(), GEOMETRY_CATEGORY_SEGMENT)) {
      double x0 = site3.point0().x();
      double y0 = site3.point0().y();
      double x1 = site3.point1().x();
      double y1 = site3.point1().y();
      double x = edge23->vertex1()->x();
      double y = edge23->vertex1()->y();
      makefoot(x, y, x0, y0, x1, y1);
      edge23->twin()->foot(x, y);
      if (oldcases) printf("OC  5   [%llu] [%f,%f]\n",(long long unsigned) edge23->twin(),x,y);

    }

    

    if (edge23->vertex1() && !new_edge2.foot() 
        && new_edge2.is_linear() 
        //&& new_edge2.is_primary() 
        && belongs(site3.source_category(), GEOMETRY_CATEGORY_SEGMENT)) {
      double x0 = site3.point0().x();
      double y0 = site3.point0().y();
      double x1 = site3.point1().x();
      double y1 = site3.point1().y();
      double x = new_vertex.x();
      double y = new_vertex.y();
      makefoot(x, y, x0, y0, x1, y1);
      new_edge2.foot(x, y);

      if (oldcases) printf("OC  6   [%llu] [%f,%f]\n",(long long unsigned) &new_edge2,x,y);
    }

    if (!edge12->foot() 
        && edge12->is_linear() 
        //&& edge12->is_primary() 
        && belongs(site1.source_category(), GEOMETRY_CATEGORY_SEGMENT)
        && edge12->vertex1()
        ) {
      double x0 = site1.point0().x();
      double y0 = site1.point0().y();
      double x1 = site1.point1().x();
      double y1 = site1.point1().y();
      double x = new_vertex.x();
      double y = new_vertex.y();
      makefoot(x, y, x0, y0, x1, y1);
      edge12->foot(x, y);

      if (oldcases) printf("OC  7.1 [%llu] [%f,%f]\n",(long long unsigned) edge12,x,y);

      if (new_vertex.r() == 0 && edge12->vertex1()) {
        edge12->twin()->foot(edge12->vertex1()->x(), edge12->vertex1()->y());
        if (oldcases) printf("OC  7.2 [%llu] [%f,%f]\n",(long long unsigned) edge12->twin(),edge12->vertex1()->x(), edge12->vertex1()->y());
      }

    }

    // didn't see change with this
    // might be redundant - or not right
    // thinking not right, though it is picking up a case that the first
    // foot finding conditional should get
    // ... yeah doesn't seem right
    if (false
        && edge12->vertex1()
        && !edge12->next()->foot() 
        && edge12->is_linear() 
        //&& edge12->is_primary() 
        && belongs(site1.source_category(), GEOMETRY_CATEGORY_SEGMENT)) {
      double x0 = site1.point0().x();
      double y0 = site1.point0().y();
      double x1 = site1.point1().x();
      double y1 = site1.point1().y();
      double x = edge12->vertex1()->x();
      double y = edge12->vertex1()->y();
      makefoot(x, y, x0, y0, x1, y1);
      edge12->next()->foot(x, y);

      if (oldcases) printf("OC  8   [%llu] [%f,%f]\n",(long long unsigned) edge12->next(),x,y);

    }
    
    // trying to fill in corner feet
    // seems to work without creating feet where they don't belong
    if (!edge12->is_primary()
        &&  edge12->vertex1()
        && !edge12->twin()->prev()->is_primary()
        &&  edge12->next()->is_primary()
        && !edge12->next()->foot()
       ) {
      edge12->next()->foot(edge12->vertex1()->x(), edge12->vertex1()->y());

      if (oldcases) printf("OC  9   [%llu] [%f,%f]\n",(long long unsigned) edge12->next(),edge12->vertex1()->x(), edge12->vertex1()->y());

    }
    if (   !edge23->is_primary()
        &&  edge23->vertex1()
        && !edge23->twin()->prev()->is_primary()
        &&  edge23->next()->is_primary()
        && !edge23->next()->foot()
       ) {
      edge23->next()->foot(edge23->vertex1()->x(), edge23->vertex1()->y());

      if (oldcases) printf("OC 10   [%llu] [%f,%f]\n",(long long unsigned) edge23->next(),edge23->vertex1()->x(), edge23->vertex1()->y());
    }

    
    // another based on special case, that hopefully has general utility
    // this worked for that one special case
    if (  !edge23->foot()
        && edge23->vertex1()
        && edge23->is_linear()
        && edge23->twin()->next()->foot()
       ) {
        double x = edge23->twin()->next()->foot()->x();
        double y = edge23->twin()->next()->foot()->y();
        reflect(x, y, edge23->vertex0()->x(), edge23->vertex0()->y(),
                      edge23->vertex1()->x(), edge23->vertex1()->y());
        edge23->foot(x, y);

        if (oldcases) printf("OC 11   [%llu] [%f,%f]\n",(long long unsigned) edge23,x,y);
    }
    // same as above but for edge12 - no test case to demonstrate it yet though
    if (  !edge12->foot()
        && edge12->vertex1()
        && edge12->is_linear()
        && edge12->twin()->next()->foot()
       ) {
        double x = edge12->twin()->next()->foot()->x();
        double y = edge12->twin()->next()->foot()->y();
        reflect(x, y, edge12->vertex0()->x(), edge12->vertex0()->y(),
                      edge12->vertex1()->x(), edge12->vertex1()->y());
        edge12->foot(x, y);

        if (oldcases) printf("OC 12   [%llu] [%f,%f]\n",(long long unsigned) edge12,x,y);
    }

    // another special case
    // got it - that was last of first round of many missing feet
    if (   !edge23->foot()
        && !edge12->is_primary()
        //&& !edge23->is_curved()
        &&  edge23->is_primary()
        &&  edge12->vertex1()
       ) {
      edge23->foot(edge12->vertex1()->x(), edge12->vertex1()->y());

      if (oldcases) printf("OC 13   [%llu] [%f,%f]\n",(long long unsigned) edge23,edge12->vertex1()->x(), edge12->vertex1()->y());

    }
    // might need similar but not same as above for similar edge12 case
    // if that's possible, but should demonstrate or illustrate need for that


    // curved edges
    
    // on t/test this might not have significant effect with or without
    // ... this fixes missing feet on full hex grid test
    if (
        //!edge12->foot() && 
        edge12->is_curved()) {
      if (belongs(site1.source_category(), GEOMETRY_CATEGORY_SEGMENT)) {
        double x0 = site1.point0().x();
        double y0 = site1.point0().y();
        double x1 = site1.point1().x();
        double y1 = site1.point1().y();
        double x = new_vertex.x();
        double y = new_vertex.y();
        makefoot(x, y, x0, y0, x1, y1);
        edge12->foot(x, y);

        if (oldcases) printf("OC 14.1 [%llu] [%f,%f]\n",(long long unsigned) edge12,x,y);

      } else if (belongs(site1.source_category(), GEOMETRY_CATEGORY_POINT)) {

        edge12->foot(site1.point0().x(), site1.point0().y());

        if (oldcases) printf("OC 14.2 [%llu] [%f,%f]\n",(long long unsigned) edge12,(double) site1.point0().x(),(double) site1.point0().y());

      }
    }
    
    if (edge12->twin()->foot() && edge12->twin()->cell()->contains_point()
        && !edge23->foot()) {
        edge23->foot(edge12->twin()->foot()->x(), edge12->twin()->foot()->y());

        if (oldcases) printf("OC 15   [%llu] [%f,%f]\n",(long long unsigned) edge23,edge12->twin()->foot()->x(), edge12->twin()->foot()->y());

    }

    if (edge23->foot() && edge23->cell()->contains_point()
        && edge23->next() && !edge23->next()->foot() 
        && edge23->next()->is_primary()) {
        // rare
        edge23->next()->foot(edge23->foot()->x(), edge23->foot()->y());

        if (oldcases) printf("OC 16   [%llu] [%f,%f]\n",(long long unsigned) edge23->next(),edge23->foot()->x(), edge23->foot()->y());

    }

    if (edge23->twin()->foot() && edge23->twin()->cell()->contains_point()
        && !new_edge2.foot()) {
        new_edge2.foot(edge23->twin()->foot()->x(), edge23->twin()->foot()->y());

        if (oldcases) printf("OC 17   [%llu] [%f,%f]\n",(long long unsigned) &new_edge2,edge23->twin()->foot()->x(), edge23->twin()->foot()->y());

    }
    
    if (edge12->foot() && edge12->cell()->contains_point()
        && !new_edge1.foot()) {
        new_edge1.foot(edge12->foot()->x(), edge12->foot()->y());

        if (oldcases) printf("OC 18   [%llu] [%f,%f]\n",(long long unsigned) &new_edge1,edge12->foot()->x(), edge12->foot()->y());

    }

    // derived from lingering special case on hex fan housing (way above hex mesh)
    // ... but also handles at least three feet in t/test, so this is legit
    // ... yes this is happening in several cases that are also handled by other
    // cases above, and then a few those miss
    if (
       !edge23->foot()
       && edge23->is_curved() && new_edge2.is_curved()
       && edge23->twin()->cell()->contains_point()
       && new_edge2.cell()->contains_point()
       && edge12->foot()
       ) {
      double x = edge12->foot()->x();
      double y = edge12->foot()->y();
      reflect(x, y, edge12->vertex0()->x(), edge12->vertex0()->y(),
                   edge12->vertex1()->x(), edge12->vertex1()->y());
      edge23->foot(x, y);

      if (oldcases) printf("OC 19   [%llu] [%f,%f]\n",(long long unsigned) edge23,x,y);

    }

    // SIMPLE CURVED SEGEMENTS FOR
    // new_edge2
    // THAT SHOULD EASILY BE HANDLED NOW BY NEW CASE CODE
    // but I'm still missing some up there
    
    if (!is_linear) {
      if (   belongs(site1.source_category(), GEOMETRY_CATEGORY_POINT)
          && belongs(site3.source_category(), GEOMETRY_CATEGORY_SEGMENT)
         ) {
        if (new_edge2.vertex0()) {
          double x0 = site3.point0().x();
          double y0 = site3.point0().y();
          double x1 = site3.point1().x();
          double y1 = site3.point1().y();
          double x = new_edge2.vertex0()->x();
          double y = new_edge2.vertex0()->y();

          makefoot(x, y, x0, y0, x1, y1);
         
          if (oldcases) printf("OC 20   [%llu] [%f,%f]\n",(long long unsigned) &new_edge2,x,y);
          if (oldcases && new_edge2.foot()) printf("   (ignore because foot already set)\n");
          if (!new_edge2.foot()) new_edge2.foot(x, y);

        }
      }
      
      if (   belongs(site3.source_category(), GEOMETRY_CATEGORY_POINT)
          && belongs(site1.source_category(), GEOMETRY_CATEGORY_SEGMENT)) {
        
        if (oldcases) printf("OC 21   [%llu] [%f,%f]\n",(long long unsigned) &new_edge2,(double) site3.point0().x(),(double) site3.point0().y());
        if (oldcases && new_edge2.foot()) printf("   (ignore because foot already set)\n");

        if (!new_edge2.foot()) new_edge2.foot(site3.point0().x(), site3.point0().y());

      }
    }

    ///////////////////
    // END OLD CASES //
    ///////////////////    
    
    
    //if (new_edge1.foot() != NULL) {printf("FF 1 ");}
    //if (new_edge2.foot() != NULL) {printf("FF 2 [%d, %d]",new_edge2.foot()->x(),new_edge2.foot()->y());}
    //printf("\n");


    // Return a pointer to the new half-edge.
    return std::make_pair(&new_edge1, &new_edge2);
  }

  void _build() {
    //printf("_build\n");
    // do foot check first thing, so pointer numbers match
    // up with any reports above (before stock stuff below has a chance
    // to shift stuff around in edges_ vector)
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      if (
          //it->is_primary() && 
          it->foot() == NULL
         ) {
        printf("[%llu] : NO FOOT 1: \n",(long long unsigned) &(*it));
        if (it->vertex0()) printf(" v0[%f %f]\n",it->vertex0()->x(),it->vertex0()->y());
        if (it->vertex1()) printf(" v1[%f %f]\n",it->vertex1()->x(),it->vertex1()->y());
        it->foot(0,0); // avoid segfault, give it something bad instead of NULL
      }
    }

    if (0) { // handy data dump to copy-paste between stages while debugging
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      printf("edge %lld: %lld, %lld, %lld, %ld, %ld, %s, %s, %s, %s, %s, %s, v%lld, \n    [%d,%d], \n    [%d,%d]\n",
          (long long unsigned int) &(*it),
          (long long unsigned int) it->twin(),
          (long long unsigned int) it->next(),
          (long long unsigned int) it->prev(),
          it->color(),
          it->cell()->source_index(),
          it->is_curved()?"curved":"      ",
          it->is_finite()?"FIN":"   ",
          it->is_primary()?"PRIM":"    ",
          it->twin() == it->next() ? "next=twin":"twok",
          &(*it) == it->next() ? "next=itself":"nxtok",
          it->cell()->contains_segment() ? "SEG":(it->cell()->contains_point() ? "PNT":"   "),
          (long long unsigned int) it->vertex0(),
          it->vertex0() ? (int) it->vertex0()->x():0,
          it->vertex0() ? (int) it->vertex0()->y():0,
          it->vertex1() ? (int) it->vertex1()->x():0,
          it->vertex1() ? (int) it->vertex1()->y():0
      );
    }
    } // end handy data dump

    // Remove degenerate edges.
    edge_iterator last_edge = edges_.begin();
    for (edge_iterator it = edges_.begin(); it != edges_.end(); it += 2) {
      const vertex_type* v1 = it->vertex0();
      const vertex_type* v2 = it->vertex1();
      if (v1 && v2 && vertex_equality_predicate_(*v1, *v2)) {
        remove_edge(&(*it));
      }
      else {
        if (it != last_edge) {
          edge_type* e1 = &(*last_edge = *it);
          edge_type* e2 = &(*(last_edge + 1) = *(it + 1));

          e1->twin(e2);
          e2->twin(e1);
          if (e1->prev()) {
            e1->prev()->next(e1);
            e2->next()->prev(e2);
          }
          if (e2->prev()) {
            e1->next()->prev(e1);
            e2->prev()->next(e2);
          }
        }
        last_edge += 2;
      }
    }
    edges_.erase(last_edge, edges_.end());

    // Set up incident edge pointers for cells and vertices.
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      it->cell()->incident_edge(&(*it));
      if (it->vertex0()) {
        it->vertex0()->incident_edge(&(*it));
      }
    }

    // Remove degenerate vertices.
    vertex_iterator last_vertex = vertices_.begin();
    for (vertex_iterator it = vertices_.begin(); it != vertices_.end(); ++it) {
      if (it->incident_edge()) {
        if (it != last_vertex) {
          *last_vertex = *it;
          vertex_type* v = &(*last_vertex);
          edge_type* e = v->incident_edge();
          do {
            e->vertex0(v);
            e = e->rot_next();
          } while (e != v->incident_edge());
        }
        ++last_vertex;
      }
    }
    vertices_.erase(last_vertex, vertices_.end());

    // Set up next/prev pointers for infinite edges.
    if (vertices_.empty()) {
      if (!edges_.empty()) {
        // Update prev/next pointers for the line edges.
        edge_iterator edge_it = edges_.begin();
        edge_type* edge1 = &(*edge_it);
        edge1->next(edge1);
        edge1->prev(edge1);
        ++edge_it;
        edge1 = &(*edge_it);
        ++edge_it;

        while (edge_it != edges_.end()) {
          edge_type* edge2 = &(*edge_it);
          ++edge_it;

          edge1->next(edge2);
          edge1->prev(edge2);
          edge2->next(edge1);
          edge2->prev(edge1);

          edge1 = &(*edge_it);
          ++edge_it;
        }

        edge1->next(edge1);
        edge1->prev(edge1);
      }
    } else {
      // Update prev/next pointers for the ray edges.
      for (cell_iterator cell_it = cells_.begin();
         cell_it != cells_.end(); ++cell_it) {
        if (cell_it->is_degenerate())
          continue;
        // Move to the previous edge while
        // it is possible in the CW direction.
        edge_type* left_edge = cell_it->incident_edge();
        while (left_edge->prev() != NULL) {
          left_edge = left_edge->prev();
          // Terminate if this is not a boundary cell.
          if (left_edge == cell_it->incident_edge())
            break;
        }

        if (left_edge->prev() != NULL)
          continue;

        edge_type* right_edge = cell_it->incident_edge();
        while (right_edge->next() != NULL)
          right_edge = right_edge->next();
        left_edge->prev(right_edge);
        right_edge->next(left_edge);
      }
    }

    // The above gets us the complete Voronoi diagram.
    // Now we'll convert that to the medial axis.
    
    if (0) { // handy data dump to copy-paste between stages while debugging
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      printf("edge %lld: %lld, %lld, %lld, %ld, %ld, %s, %s, %s, %s, %s, %s, v%lld, [%d,%d]\n",
          (long long unsigned int) &(*it),
          (long long unsigned int) it->twin(),
          (long long unsigned int) it->next(),
          (long long unsigned int) it->prev(),
          it->color(),
          it->cell()->source_index(),
          it->is_curved()?"curved":"      ",
          it->is_finite()?"FIN":"   ",
          it->is_primary()?"PRIM":"    ",
          it->twin() == it->next() ? "next=twin":"twok",
          &(*it) == it->next() ? "next=itself":"nxtok",
          it->cell()->contains_segment() ? "SEG":(it->cell()->contains_point() ? "PNT":"   "),
          (long long unsigned int) it->vertex0(),
          it->vertex0() ? (int) it->vertex0()->x():0,
          it->vertex0() ? (int) it->vertex0()->y():0
      );
    }
    } // end handy data dump


    ///////////////////////////////////////////////////////////////
    //     REPAIR DEGENERACIES FROM COLLINEAR INPUT SEGMENTS     //
    ///////////////////////////////////////////////////////////////

    // We found that you get infinite secondaries passing right through
    // the place where two collinear input segments meet. This makes their
    // terminal vertex appear on an edge "inside" the polygon. There's
    // no Voronoi vertex  at the polygon segment ends, as there would be for
    // a normal corner of the polygon.
    // This is a weird situation for the medial axis.
    // It implies that there is a maximal inscribed circle centered at that
    // vertex inside the polygon with infinite radius, simply because the angle 
    // between the associated input segments is zero. It makes more sense to 
    // treat these "flat corners" like any other corners of the polygon.
    // We will assume that collinear input segments contain information
    // that should be preserved. They will not be considered redundant.
    // Users may screen out collinear segments themselves if they consider
    // them redundant. 
    
    // iterate backwards through the vertex array, because we add vertices to 
    // the end, and we don't want to recheck them when we get there

    for (vertex_iterator it = vertices_.end()-1; it != vertices_.begin() - 1; --it) {
        edge_type * start_edge = it->incident_edge();
        edge_type * edge = start_edge;

        do {

          // Make a vertex out of the foot - first so we can use it with
          // vertex_equality_predicate_(), and later, if we satisfy y any cases
          // here, to insert into the Voronoi diagram.
          const vertex_type& foot_as_vertex = vertex_type(edge->foot()->x(),
                                                          edge->foot()->y(),0);
          
          //printf("4FOOT: [%f, %f]\n",foot_as_vertex.x(),foot_as_vertex.y());
          /*
          printf("[%d][%d][%d][%d][%d][%d][%d][%d][%d]\n",
              edge->next() != NULL ? 1:0,
              edge->prev() != NULL ? 1:0,
              edge->next() == edge->prev() ? 1:0,
              edge->vertex0() ? 1:0, // the outgoing one
              edge->cell()->contains_point() ? 1:0,
              edge->next()->cell()->contains_point() ? 1:0,

              // primaries on either side of the outer secondaries
              edge->twin()->next()->is_primary() ? 1:0,
              edge->next()->twin()->prev()->is_primary() ? 1:0,
              edge->twin()->next()->twin() != edge->next()->twin()->prev() ? 1:0 // ahhhh.
          );
          */

          if (
               edge->is_secondary()
            && edge->next() != edge->prev() // avoid flat corner case
            && edge->twin()->next() != edge->twin()->prev() // avoid flat corner case
            && edge->foot()
            && ! ( (!edge->vertex0() ||
                    vertex_equality_predicate_(foot_as_vertex, *edge->vertex0())
                   )
                   ||
                   (!edge->vertex1() ||
                    vertex_equality_predicate_(foot_as_vertex, *edge->vertex1())
                   )
                 )
             ) {

            // somehow that only gets finite edges, and infinites are
            // picked up by the next elsif
            
            // Those grazing secondaries happen on finite secondaries too.
            // Since each secondary now has a foot, you can check if
            // the foot matches one of the ends in the finite case.
            // If not, insert a vertex that is that foot.
            // ends up being alot like the other two cases, and those are
            // pretty much the same.
            // do the ifs for all three cases, but just to set up
            // the new vertex and the two new edges.
            // Then, if those are defined, do the identical link up code.

            //printf("GOT IT 3 [%s]\n",edge->is_infinite() ? "INF":"FIN");
            //printf("  foot: [%f, %f]\n",foot_as_vertex.x(),foot_as_vertex.y());

            // Make a new Voronoi vertex at the "flat corner" on the polygon.
            vertices_.push_back(foot_as_vertex);
            vertex_type & new_vertex = vertices_.back();

            edges_.push_back(edge_type(true,false));
            edge_type & new_edge1 = edges_.back();
            edges_.push_back(edge_type(true,false));
            edge_type & new_edge2 = edges_.back();
                        
            new_vertex.incident_edge(&new_edge2);

            edge_type * e_out   = edge;
            edge_type * e_in  = edge->twin();

            if ( it->incident_edge() == e_out ) {
              //printf("0 V[%llu] set inc edge\n",(long long unsigned) &(*it));
              it->incident_edge(&new_edge1);
            }

            // link up new secondary edge going out
            new_edge1.twin(&new_edge2);
            new_edge1.prev(e_out->prev());
            new_edge1.next(e_out);
            new_edge1.vertex0(&(*it));
            new_edge1.cell(e_out->cell());
            new_edge1.foot(new_vertex.x(),new_vertex.y());

            // link up new primary edge coming in
            new_edge2.twin(&new_edge1);
            new_edge2.prev(e_in);
            new_edge2.next(e_in->next());
            new_edge2.vertex0(&new_vertex);
            new_edge2.cell(e_in->cell());
            new_edge2.foot(new_vertex.x(),new_vertex.y());

            // link up old outer infinites to new primaries
            e_out->prev(&new_edge1);
            e_in->next(&new_edge2);

            // link up old surrounding primaries to new edges
            new_edge1.prev()->next(&new_edge1);
            new_edge2.next()->prev(&new_edge2);

            // set new vertex0 for old infinites
            e_out->vertex0(&new_vertex);

            // fix edge reference situation for our do{}while() rot_next()
            start_edge = &new_edge1;
            edge = &new_edge1;

          }

          else if (
                 edge->is_secondary()
              && edge->is_infinite()
              && edge->vertex0()
              && edge->rot_next()->is_primary()
              //&& edge->rot_next()->is_infinite()
              && edge->rot_prev()->is_primary()
              //&& edge->rot_prev()->is_finite()
              //&& edge->rot_prev()->is_linear()
             ) {

            //printf("\nGOT IT 1\n\n");

            if (edge->foot() != NULL) {
              //printf("  foot: [%f, %f]\n",foot_as_vertex.x(),foot_as_vertex.y());

              //printf("Inserting edge pair and vertex to repair 'flat corner'\n");

              // Make a new Voronoi vertex at the "flat corner" on the polygon.
              vertices_.push_back(foot_as_vertex);

              // recheck that this copied code from "got it 2"
              // is valid for this "GOT IT 1" case
              // yeah think it ended up being the exact same code
              // so fixup the if()s to make it really one case

              vertex_type & new_vertex = vertices_.back();

              // Make two new secondary edges that will go from
              // the current vertex out to the new vertex,
              // and move the two old infinite edges at the current vertex
              // out to the new vertex.

              edges_.push_back(edge_type(true,false));
              edge_type & new_edge1 = edges_.back();
              edges_.push_back(edge_type(true,false));
              edge_type & new_edge2 = edges_.back();
                        
              new_vertex.incident_edge(&new_edge2);

              edge_type * e_out = edge;
              edge_type * e_in  = edge->twin();

              if ( it->incident_edge() == e_out ) {
                //printf("1 V[%llu] set inc edge\n",(long long unsigned) &(*it));
                it->incident_edge(&new_edge1);
              }

              // link up new secondary edge going out
              new_edge1.twin(&new_edge2);
              new_edge1.prev(e_out->prev());
              new_edge1.next(e_out);
              new_edge1.vertex0(&(*it));
              new_edge1.cell(e_out->cell());
              new_edge1.foot(new_vertex.x(),new_vertex.y());

              // link up new primary edge coming in
              new_edge2.twin(&new_edge1);
              new_edge2.prev(e_in);
              new_edge2.next(e_in->next());
              new_edge2.vertex0(&new_vertex);
              new_edge2.cell(e_in->cell());
              new_edge2.foot(new_vertex.x(),new_vertex.y());

              // link up old outer infinites to new primaries
              e_out->prev(&new_edge1);
              e_in->next(&new_edge2);

              // link up old surrounding primaries to new edges
              new_edge1.prev()->next(&new_edge1);
              new_edge2.next()->prev(&new_edge2);

              // set new vertex0 for old infinites
              e_out->vertex0(&new_vertex);

              // fix edge reference situation for our do{}while() rot_next()
              start_edge = &new_edge1;
              edge = &new_edge1;



            } else {
              // should not be possible, but we may still miss some feet above
              printf("[%llu] : NO FOOT 2: %s \n\n",
                (long long unsigned) edge,
                edge->is_primary() ? "pri" : "2nd");
            }

          }

          else if ( edge->is_secondary()
              && edge->is_infinite()
              && edge->vertex0()
              && edge->rot_prev()->is_primary()
              //&& edge->rot_prev()->is_infinite()
              && edge->rot_next()->is_primary()
              //&& edge->rot_next()->is_finite()
              //&& edge->rot_next()->is_linear()
             ) {

            //printf("\nGOT IT 2\n\n");
            
            if (edge->foot() != NULL) {
              //printf("  foot: [%f, %f]\n",foot_as_vertex.x(),foot_as_vertex.y());


              // Make a new Voronoi vertex at the "flat corner" on the polygon.
              vertices_.push_back(foot_as_vertex);

              vertex_type & new_vertex = vertices_.back();

              // Make two new secondary edges that will go from
              // the current vertex out to the new vertex,
              // and move the two old infinite edges at the current vertex
              // out to the new vertex.
              
              edges_.push_back(edge_type(true,false));
              edge_type & new_edge1 = edges_.back();
              edges_.push_back(edge_type(true,false));
              edge_type & new_edge2 = edges_.back();
                        
              new_vertex.incident_edge(&new_edge2);

              edge_type * e_out   = edge;
              edge_type * e_in  = edge->twin();

              if ( it->incident_edge() == e_out ) {
                it->incident_edge(&new_edge1);
              }

              // link up new secondary edge going out
              new_edge1.twin(&new_edge2);
              new_edge1.prev(e_out->prev());
              new_edge1.next(e_out);
              new_edge1.vertex0(&(*it));
              new_edge1.cell(e_out->cell());
              new_edge1.foot(new_vertex.x(),new_vertex.y());

              // link up new primary edge coming in
              new_edge2.twin(&new_edge1);
              new_edge2.prev(e_in);
              new_edge2.next(e_in->next());
              new_edge2.vertex0(&new_vertex);
              new_edge2.cell(e_in->cell());
              new_edge2.foot(new_vertex.x(),new_vertex.y());

              // link up old outer infinites to new primaries
              e_out->prev(&new_edge1);
              e_in->next(&new_edge2);

              // link up old surrounding primaries to new edges
              new_edge1.prev()->next(&new_edge1);
              new_edge2.next()->prev(&new_edge2);

              // set new vertex0 for old infinites
              e_out->vertex0(&new_vertex);

              // fix edge reference situation for our do{}while() rot_next()
              start_edge = &new_edge1;
              edge = &new_edge1;


            } else {
              // should not be possible, but we may still miss some feet above
              printf("[%llu] : NO FOOT 3: %s \n\n",
                (long long unsigned) edge,
                edge->is_primary() ? "pri" : "2nd");
            }

          }
          else if (   edge->is_secondary()
              // characterizing the two inner secondaries
              && edge->next() != NULL
              && edge->prev() != NULL
              && edge->next() == edge->prev()
              && edge->vertex0() // the outgoing one
              // you would think... but had case where these were both false
              // and the rest seems to work without these, so don't use
              // but leave here so you remember not to use
              // ... ahh, think it happens when there's a 
              // zero-length segment in the the input. So this tolerates that.
              //&& edge->cell()->contains_point()
              //&& edge->next()->cell()->contains_point()

              // primaries on either side of the outer secondaries
              && edge->twin()->next()->is_primary()
              && edge->next()->twin()->prev()->is_primary()
              && edge->twin()->next()->twin() != edge->next()->twin()->prev() // ahhhh.
             ) {

            //printf("Found degenerate edges for 'flat corner' of polygon.\n");

            if (edge->twin()->next()->foot() != NULL) {

              //printf("Inserting edge pair and vertex to repair 'flat corner'\n");
              //printf("  foot: [%f, %f]\n",foot_as_vertex.x(),foot_as_vertex.y());

              // Make a new Voronoi vertex at the "flat corner" on the polygon.
              //vertices_.push_back(vertex_type(edge->twin()->next()->foot()->x(),
              //                                edge->twin()->next()->foot()->y(),
              //                                0));
              vertices_.push_back(foot_as_vertex);
              
              vertex_type & new_vertex = vertices_.back();

              // Make two new primary edges that will go from
              // the current vertex out to the new vertex,
              // and move the four infinite corner edges at the current vertex
              // out to the new vertex. This makes the new vertex look like a
              // convex corner of the polygon, topologically. Downstream code
              // that finds or operates on convex corners should be able to 
              // deal with these "flat corners" in the same manner.
              
              edges_.push_back(edge_type(true,true));
              edge_type & new_edge1 = edges_.back();
              edges_.push_back(edge_type(true,true));
              edge_type & new_edge2 = edges_.back();
                        
              new_vertex.incident_edge(&new_edge2);

              edge_type * e_out_inside   = edge;
              edge_type * e_in_inside  = edge->next();
              edge_type * e_out_outside  = edge->next()->twin();
              edge_type * e_in_outside = edge->twin();

              if (   it->incident_edge() == e_out_inside
                  || it->incident_edge() == e_out_outside
                 ) {
                it->incident_edge(&new_edge1);
              }

              // link up new primary edge going out
              new_edge1.twin(&new_edge2);
              new_edge1.prev(e_out_outside->prev());
              new_edge1.next(e_out_outside);
              new_edge1.vertex0(&(*it));
              new_edge1.cell(e_out_outside->cell());
              new_edge1.foot(new_vertex.x(),new_vertex.y());

              // link up new primary edge coming in
              new_edge2.twin(&new_edge1);
              new_edge2.prev(e_in_outside);
              new_edge2.next(e_in_outside->next());
              new_edge2.vertex0(&new_vertex);
              new_edge2.cell(e_in_outside->cell());
              new_edge2.foot(new_vertex.x(),new_vertex.y());

              // link up old outer infinites to new primaries
              e_out_outside->prev(&new_edge1);
              e_in_outside->next(&new_edge2);

              // link up old surrounding primaries to new edges
              new_edge1.prev()->next(&new_edge1);
              new_edge2.next()->prev(&new_edge2);

              // set new vertex0 for old infinites
              e_out_outside->vertex0(&new_vertex);
              e_out_inside->vertex0(&new_vertex);


              // fix edge reference situation for our do{}while() rot_next()
              start_edge = &new_edge1;
              edge = &new_edge1;

            } else {
              // should not be possible, but we may still miss some feet above
              printf("[%llu] : NO FOOT 4: %s \n\n", 
                (long long unsigned) edge->twin()->next(),
                edge->twin()->next()->is_primary() ? "pri" : "2nd");
            }

          }
          
          //printf("v[%llu] e[%llu]\n",(long long unsigned int) &(*it), (long long unsigned int) edge);
          
          edge = edge->rot_next();
        
        }
        while (edge != start_edge);
    }

    //////////////////////////////////////////////////////////////////////
    //     END OF REPAIR DEGENERACIES FROM COLLINEAR INPUT SEGMENTS     //
    //////////////////////////////////////////////////////////////////////


    if (0) { // handy data dump to copy-paste between stages while debugging
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      printf("edge %lld: %lld, %lld, %lld, %ld, %ld, %s, %s, %s, %s, %s, %s, v%lld\n",
          (long long unsigned int) &(*it),
          (long long unsigned int) it->twin(),
          (long long unsigned int) it->next(),
          (long long unsigned int) it->prev(),
          it->color(),
          it->cell()->source_index(),
          it->is_curved()?"curved":"      ",
          it->is_finite()?"FIN":"   ",
          it->is_primary()?"PRIM":"    ",
          it->twin() == it->next() ? "next=twin":"twok",
          &(*it) == it->next() ? "next=itself":"nxtok",
          it->cell()->contains_segment() ? "SEG":(it->cell()->contains_point() ? "PNT":"   "),
          (long long unsigned int) it->vertex0()
      );
    }
    }

    ///////////////////////////////////////////////////
    //    DETERMINE INTERIOR / EXTERIOR SIDEDNESS    //
    ///////////////////////////////////////////////////
    // Recursively mark edges exterior to the polygon,
    // starting with all the infinite edges.
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      if (it->is_infinite() && it->vertex0() == NULL) { //incoming from infinity
        mark_exterior(&(*it));
      }
    }
    
    // TODO : See if the two loops below can be done in one
    //        (would need mark_exterior() to set up verts_todo)
    //        or see if you can get rid of seperate mark_exterior(),
    //        and just do all this with one loop, seeded with the inf edges.

    std::queue<vertex_type*> verts_todo;

    for (vertex_iterator it = vertices_.begin(); it != vertices_.end(); ++it) {
      if (it->r() == 0) {
        edge_type * se = it->incident_edge();
        edge_type * e = se;
        bool do_v = false;
        bool side = false;
        do {
          if ( ! do_v && (e->is_internal() || e->is_external())) {
            do_v = true;
            se = e;
            side = ! e->is_internal();
            //printf("GONNA\n");
          }
          if (do_v && !(e->is_internal() || e->is_external())) {
            //printf("[%llu] : edge todo\n[%f %f %f] [%f %f %f]\n",
            //  (long long unsigned) e,
            //  e->vertex0()->x(),e->vertex0()->y(),e->vertex0()->r(),
            //  e->vertex1()->x(),e->vertex1()->y(),e->vertex1()->r());
            //printf("premark  : %s%s \n",e->is_internal()?"I":"0",e->is_external()?"X":"0");
            mark_side(side, e, verts_todo);
            //printf("postmark : %s%s \n",e->is_internal()?"I":"0",e->is_external()?"X":"0");
          }
          e = e->rot_next();
        } while (e != se) ;                
      }
    }
    
    //printf("verts_todo size %lu\n",verts_todo.size());

    //int lim = 3;
    while (
           //lim-- > 0 && 
           !verts_todo.empty()
          ) {
      int todo = verts_todo.size();
      //printf("\n##################### BATCH %d (size: %d )\n\n",lim, todo);
      for (int i = 0; i < todo; i++) {
        vertex_type * v = verts_todo.front();
        edge_type * se = v->incident_edge();
        edge_type * e = se;
        int rotcnt = 0;
        do {
          rotcnt++;
          //printf("[%llu] : %s%s \n",(long long unsigned) e, e->is_internal()?"I":"0", e->is_external()?"X":"0");
          e = e->rot_next();
        } while (e != se) ;
        //printf("VERT %d, rot count: %d \n",i, rotcnt);
        
        bool do_v = false;
        bool side = false;
        e = se;
        do {
          if ( (!do_v) && (e->is_internal() || e->is_external())) {
            do_v = true;
            se = e;
            side = e->is_internal() ? true : false;
            side = e->vertex0()->r() == 0 ? !side:side;
            //printf("GONNA %s%s \n",e->is_internal()?"I":"0",e->is_external()?"X":"0");
          }
          //printf("%s : %s%s \n",do_v?"DO":"no",e->is_internal()?"I":"0",e->is_external()?"X":"0");
          if (do_v && !(e->is_internal() || e->is_external())) {
            //printf("mark side while\n");
            mark_side(side, e, verts_todo);
            //printf("[%llu] : set side: %s|%s \n",(long long unsigned) e, e->is_internal() ? "I":"0", e->is_external() ? "X":"0");
          }
          e = e->rot_next();
        } while (e != se) ;
        verts_todo.pop();
      }
    
    //printf("while verts_todo size %lu\n", verts_todo.size());

    }

    /////////////////////////////////////////////////////
    //    CHANGE VORONOI LINKS TO MEDIAL AXIS LINKS    //
    /////////////////////////////////////////////////////

    /////////////
    // At this point we modify the half edge graph to better represent the 
    // the medial axis.
    // The main thing to do is update next() and prev() pointers to follow
    // along the primary edges that represent the medial axis, instead
    // of having them point just to the next/prev within each Voronoi cell.
    //
    // Just look at every vertex with radius==0 and
    // rot_next around it, connecting every edge to it's twin.
    // That's kind of the main difference between the Voronoi diagram and
    // medial axis half edge graphs - VD passes through r==0 nodes, while 
    // MA bounces off of them.

    for (vertex_iterator it = vertices_.begin(); it != vertices_.end(); ++it) {
      if (it->r() == 0) {
        edge_type * se = it->incident_edge();
        edge_type * e = se;
        edge_type * ep;
        do {
          ep = e->prev();
          // preserve old Voronoi links
          e->prev_vd(ep);
          e->twin()->next_vd(e->twin()->next());
          // make new medial axis links
          e->prev(e->twin());
          e->twin()->next(e);
          // do e->rot_next(), (e->prev()->twin()) using old prev() link
          e = ep->twin(); 
        } while (e != se);
      }
    }

    // 

    // theta is the angle direction of the tangent vector at the beginning of the edge.
    //  For curved edges, twins will have different thetas.
    //  Thetas at edge meetings don't always match.
    // phi is the angle between the tangent at the start of the edge (theta)
    //  and the foot.

    // TODO: guess these should be constants set way above somewhere
    const double pi           = 4.0 * atan2(1.0, 1.0);
    const double pi_over_two  = 2.0 * atan2(1.0, 1.0);
    const double pi_over_four =       atan2(1.0, 1.0);
    const double two_pi       = 8.0 * atan2(1.0, 1.0);

    // theta calc (edge initial tangent angle)
    // based on finding slope of directrix first (for curved edges obviously)
    // and using difference in radius lengths at ends of curved edge (parabolic arc)

    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      if (it->is_finite()) {
        if (it->is_curved()) {

          edge_type * edge_with_point = it->cell()->contains_point()
                                      ? &(*it)
                                      : it->twin();

          //# A curved primary segment should always have secondaries
          //# before and after it, coming from the parabola focus and
          //# going to it. Around vertices at the ends of infinite segments,
          //# one of these secondaries may be missing, so use the other.

          edge_type * secondary_with_point = edge_with_point;

          vertex_type * focus;
            
          if (   edge_with_point->prev_vd()->is_secondary() 
              && edge_with_point->prev_vd()->is_finite()
              && edge_with_point->prev_vd()->prev_vd()->is_secondary() 
              && edge_with_point->prev_vd()->prev_vd()->is_finite()
              && edge_with_point->prev_vd()->vertex0()
                 == edge_with_point->prev_vd()->prev_vd()->vertex1()
             ) {
              focus = edge_with_point->prev_vd()->vertex0();
          }
          else if (edge_with_point->next_vd()->is_secondary()
              && edge_with_point->next_vd()->is_finite()
              && edge_with_point->next_vd()->next_vd()->is_secondary() 
              && edge_with_point->next_vd()->next_vd()->is_finite()
              && edge_with_point->next_vd()->vertex1()
                 == edge_with_point->next_vd()->next_vd()->vertex0()
             ) {
             focus = edge_with_point->next_vd()->vertex1();               
          }
          // sometimes the relavent secondary is not immediately adjacent
          else {
              edge_type * se = secondary_with_point;
              do {
                  secondary_with_point = secondary_with_point->prev_vd();
              } while (
                   secondary_with_point->is_primary()
                && secondary_with_point != edge_with_point
                && secondary_with_point != se);
              if (secondary_with_point == se) {
                printf("[%llu] : failed focus finding while determining theta (1)\n",
                  (long long unsigned) &(*it));
              }
              if (secondary_with_point->vertex0()) {
                focus = secondary_with_point->vertex0();
              }
              else {
                // search going the other way, next instead of prev
                secondary_with_point = edge_with_point;
                focus = secondary_with_point->vertex0();
                se = secondary_with_point;
                do {
                    secondary_with_point = secondary_with_point->next_vd();
                } while (
                     secondary_with_point->is_primary()
                  && secondary_with_point != edge_with_point
                  && secondary_with_point != se);
                if (secondary_with_point == se) {
                  printf("[%llu] : failed focus finding while determining theta (2)\n",
                    (long long unsigned) &(*it));
                }
                if (secondary_with_point->vertex1()) {
                  focus = secondary_with_point->vertex1();
                }
                else {
                  printf("[%llu] : failed focus finding while determining theta (3)\n",
                    (long long unsigned) &(*it));
                }
              
              }
              
          }

          vertex_type * p0 = edge_with_point->vertex0();
          vertex_type * p2 = edge_with_point->vertex1();
            
          //# directrix angle = alpha - gamma
          //# alpha = angle of line p0 to p2
          //# gamma = angle opposite delta_r in a right triangle with b as 
          //#         base, r as height, p0 to p2 line as hypotenuse
          //# base b is derived from known delta_r and hypotenuse
          double delta_x = p2->x() - p0->x();
          double delta_y = p2->y() - p0->y();
          double alpha = atan2(delta_y, delta_x);
            
          double delta_r = (  sqrt(pow(p2->x() - focus->x(),2) + pow(p2->y() - focus->y(),2))
                            - sqrt(pow(p0->x() - focus->x(),2) + pow(p0->y() - focus->y(),2)));
          double h = sqrt(pow(delta_x,2) + pow(delta_y,2));
          double b = sqrt(pow(h,2) - pow(delta_r,2));
          double gamma = atan2(delta_r, b);

          double directrix_angle = alpha - gamma;

          while (directrix_angle >  pi) {directrix_angle -= two_pi;}
          while (directrix_angle < -pi) {directrix_angle += two_pi;}
            
          //# Now ready to figure the tangent angle at requested edge's 
          //# vertex0. Although, if we had to use the twin to figure the 
          //# above, we're really figuring the angle at the twin's vertex1, 
          //# which is 180 degrees opposite of what we want.
          //# So we'll turn it around when we're done, if that's that case.
            
          double angle_focus_p = atan2(it->vertex0()->y() - focus->y(),
                                       it->vertex0()->x() - focus->x());

          double angle_focus_p_un_dir_rot = angle_focus_p - directrix_angle;
          while (angle_focus_p_un_dir_rot >  pi) {angle_focus_p_un_dir_rot -= two_pi;}
          while (angle_focus_p_un_dir_rot < -pi) {angle_focus_p_un_dir_rot += two_pi;}

          //# Fixes 180 degree wrongness in 2nd quadrant.
          //# Consider that this signed angle $angle_focus_p_un_dir_rot has value -90 degrees
          //# when the point is at the bottom of the parabola. We then have a 180 degree sweep
          //# both positive and negative from there, for the parabola's two legs (of course
          //# never reaching those limits). So on the right side that's -90 to 90.
          //# On the left that's -90 to -270. But atan2() isn't going to give us the negative
          //# angles we want in the 2nd quadrant. So convert all positive 2nd quadrant angles to
          //# to negatives by subtracting 2*PI.
          if (angle_focus_p_un_dir_rot > pi_over_two) {
            angle_focus_p_un_dir_rot -= two_pi;
          }

          //# sketched this a couple ways for + and - angles and seems right
          double theta_un_dir_rot = pi_over_four + angle_focus_p_un_dir_rot/2;
            
          double theta = theta_un_dir_rot + directrix_angle;

          //# The 180 degree turn-around, if needed.
          if (!it->cell()->contains_point()) {
            theta += pi;
          }
          while (theta >  pi) {theta -= two_pi;}
          while (theta < -pi) {theta += two_pi;}

          it->theta(theta);

        }
        // linear case
        else {
          it->theta(atan2(it->vertex1()->y() - it->vertex0()->y(),
                          it->vertex1()->x() - it->vertex0()->x()));
        }
      }
    }

    // now for tangent angles of infinite edges

    bool inf_theta_debug = false;

    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      if (it->is_infinite()) {
        double theta = 0;
        if (it->is_primary()) {
          //# infinite primary linear
          //# always has a primary behind that we can get the angle from
          //# 3 main cases

          edge_type * outgoing = it->vertex0() ? &(*it) : it->twin();
          
          //printf("outgoing: [%f, %f] \n",outgoing->vertex0()->x(),outgoing->vertex0()->y());
          
          // four primaries meet at a point
          if (     outgoing->rot_next_vd()->is_primary()
                && outgoing->rot_prev_vd()->is_primary()
                && outgoing->rot_next_vd()->rot_next_vd()->is_primary()
                &&    outgoing->rot_next_vd()->rot_next_vd()
                   != outgoing
                &&    outgoing->rot_next_vd()->rot_next_vd()
                   != outgoing->rot_next_vd()
                &&    outgoing->rot_next_vd()->rot_next_vd()
                   != outgoing->rot_prev_vd()
               ) {
            if (inf_theta_debug) printf("CASE P0");
            if (outgoing->rot_next_vd()->rot_next_vd()->is_finite()) {
              if (inf_theta_debug) printf(".1 ");
              theta = outgoing->rot_next_vd()->rot_next_vd()->theta() + pi;
            }
            else if ( outgoing->rot_next_vd()->is_finite()
                   && outgoing->rot_prev_vd()->is_finite()
                    ) {
              if (inf_theta_debug) printf(".2 ");
              double a1 = outgoing->rot_next_vd()->theta();
              double a2 = outgoing->rot_prev_vd()->theta();
              //# fabs() here because we're figuring the counterclockwise
              //# sweep angle, which we add half of to prev() angle to 
              //# get the a1-a2 bisector in the right direction.
              theta = a2 + (fabs(a1 - a2) / 2);
            }
            else {
              printf("unhandled infinite primary tangent angle finding case\n");
              continue;
            }
          }

          // three primaries meet at a point
          else if (outgoing->rot_next_vd()->is_primary()
                && outgoing->rot_prev_vd()->is_primary()
                &&    outgoing->rot_next_vd()->rot_next_vd()->rot_next_vd()
                   == outgoing
               ) {
            if (inf_theta_debug) printf("CASE P0.5");
            if (   outgoing->rot_next_vd()->is_finite()
                && outgoing->rot_prev_vd()->is_finite()
               ) {
              if (inf_theta_debug) printf(".1 ");
              double a1 = outgoing->rot_prev_vd()->theta();
              double a2 = outgoing->rot_next_vd()->theta();
              theta = a1 + ((two_pi - fabs(a2 - a1)) / 2);
            }
            else {
              printf("unhandled infinite primary tangent angle finding case\n");
              continue;
            }
          }

          else if (outgoing->rot_next_vd()->is_primary()) {
            if (inf_theta_debug) printf("CASE P1");
            if (outgoing->rot_prev_vd()->rot_prev_vd()->is_primary()) {
              if (   outgoing->rot_prev_vd()->rot_prev_vd()
                      == outgoing->rot_next_vd()) {
                if (inf_theta_debug) printf(".1 ");
                theta = outgoing->rot_next_vd()->theta() + pi;
              }
              else {
                if (inf_theta_debug) printf(".2 ");
                double a1 = outgoing->rot_next_vd()->theta();
                double a2 = outgoing->rot_prev_vd()->rot_prev_vd()->theta();
                theta = a2 + (fabs(a1 - a2) / 2);
              }
            }
            else {
              printf("unhandled infinite primary tangent angle finding case\n");
              continue;
            }
          }
          else if (outgoing->rot_prev_vd()->is_primary()) {
            if (inf_theta_debug) printf("CASE P2");
            if (outgoing->rot_next_vd()->rot_next_vd()->is_primary()) {
              if (   outgoing->rot_next_vd()->rot_next_vd()
                  == outgoing->rot_prev_vd()) {
                if (inf_theta_debug) printf(".1 ");
                theta = outgoing->rot_prev_vd()->theta();
              }
              else {
                if (inf_theta_debug) printf(".2 ");
                //theta = (  outgoing->rot_prev_vd()->theta()
                //         + outgoing->rot_next_vd()->rot_next_vd()->theta()
                //        ) / 2;
              double a1 = outgoing->rot_next_vd()->rot_next_vd()->theta();
              double a2 = outgoing->rot_prev_vd()->theta();
              theta = a2 + (fabs(a1 - a2) / 2);
              }
            }
            else {
              printf("unhandled infinite primary tangent angle finding case\n");
              continue;
            }
          }
          else if ( outgoing->rot_next_vd()->rot_next_vd()->is_primary()
                 && outgoing->rot_prev_vd()->rot_prev_vd()->is_primary()) {
            if (inf_theta_debug) printf("CASE P3");
            if (   outgoing->rot_next_vd()->rot_next_vd()
                == outgoing->rot_prev_vd()->rot_prev_vd()) {
              if (inf_theta_debug) printf(".1 ");
              theta = outgoing->rot_prev_vd()->rot_prev_vd()->theta() + pi;
            }
            else {
              if (inf_theta_debug) printf(".2 ");
              double a1 = outgoing->rot_next_vd()->rot_next_vd()->theta();
              double a2 = outgoing->rot_prev_vd()->rot_prev_vd()->theta();
              theta = a2 + (fabs(a1 - a2) / 2);
            }
          }
          else {
            printf("failed to find tangent angle for infinite primary edge\n");
            continue;
          }

          while (theta >  pi) { theta -= two_pi; }
          while (theta < -pi) { theta += two_pi; }

          //# turn it around if outgoing wasn't the edge, but the edge's twin
          if (!it->vertex0()) {
            if (theta >  pi) {theta -= pi;}
            if (theta <= pi) {theta += pi;}

          }



        }
        else { // is secondary
          //# infinite linear secondary

          //# convex corner (or "flat corner", fixed up above,
          //  for collinear input segments)
          if (it->next_vd() == it->prev_vd()) {
            if (inf_theta_debug) printf("CASE S1");
            if (it->vertex0()) {
              if (inf_theta_debug) printf(".1 ");
              vertex_type * p1 = it->twin()->vertex1();
              vertex_type * p2 = it->twin()->prev_vd()->vertex0();
              double seg_angle = atan2(p2->y() - p1->y(), p2->x() - p1->x());
              theta = seg_angle + pi_over_two;
            }
            else if (it->vertex1()) {
              if (inf_theta_debug) printf(".2 ");
              vertex_type * p1 = it->twin()->vertex0();
              vertex_type * p2 = it->twin()->next_vd()->vertex1();
              double seg_angle = atan2(p2->y() - p1->y(), p2->x() - p1->x());
              theta = seg_angle + pi_over_two;
            }
          }

          //# infinite secondaries on each end of an input segment
          else if ( it->vertex0() 
                 && it->next_vd()->is_secondary()
                ) {
            if (inf_theta_debug) printf("CASE S2");
            vertex_type * p1 = it->vertex0();
            vertex_type * p2 = it->next_vd()->vertex1();
            double seg_angle = atan2(p2->y() - p1->y(), p2->x() - p1->x());
            theta = seg_angle - pi_over_two;
          }
          else if ( it->vertex1() 
                 && it->prev_vd()->is_secondary()
                ) {
            if (inf_theta_debug) printf("CASE S3");
            vertex_type * p1 = it->vertex1();
            vertex_type * p2 = it->prev_vd()->vertex0();
            double seg_angle = atan2(p2->y() - p1->y(), p2->x() - p1->x());
            theta = seg_angle - pi_over_two;
          }


          //# infinite secondaries where next or prev is infinite primary
          else if (   it->vertex0()
                 && it->next_vd()->is_primary()
                ) {
            if (inf_theta_debug) printf("CASE S2P");
            if (it->twin()->prev_vd()->is_secondary()) {
              if (inf_theta_debug) printf(".1 ");
              vertex_type * p1 = it->vertex0();
              vertex_type * p2 = it->twin()->prev_vd()->vertex0();
              double seg_angle = atan2(p2->y() - p1->y(), p2->x() - p1->x());
              theta = seg_angle + pi/2;
            }
            else {printf("hit dead end in tangent finding\n");}
          }
          else if (   it->vertex1()
                 && it->prev_vd()->is_primary()
                  ) {
            if (inf_theta_debug) printf("CASE S3P");
            if (it->twin()->next_vd()->is_secondary()) {

              if (inf_theta_debug) printf(".1 ");
              vertex_type * p1 = it->vertex1();
              vertex_type * p2 = it->twin()->next_vd()->vertex1();
              double seg_angle = atan2(p2->y() - p1->y(), p2->x() - p1->x());
              theta = seg_angle + pi/2;
            }
            else {printf("hit dead end in tangent finding\n");}
          }

            
          //# a corner where on one side it's infinite and on the
          //# other it's finite - because around the corner on the
          //# finite side there's convexity - but we're interested in
          //# what's around the corner on the infinte side

          else if (   it->vertex1() 
              && it->next_vd()->is_secondary()
              && it->next_vd()->is_finite()
             ) {
            if (inf_theta_debug) printf("CASE S4");
            vertex_type * p1 = it->twin()->vertex0();
            vertex_type * p2 = it->twin()->next_vd()->vertex1();
            double seg_angle = atan2(p2->y() - p1->y(), p2->x() - p1->x());
            theta = seg_angle + pi_over_two;
          }
          else if (   it->vertex0() 
              && it->prev_vd()->is_secondary()
              && it->prev_vd()->is_finite()
             ) {
            if (inf_theta_debug) printf("CASE S5");
            vertex_type * p1 = it->twin()->vertex1();
            vertex_type * p2 = it->twin()->prev_vd()->vertex0();
            double seg_angle = atan2(p2->y() - p1->y(), p2->x() - p1->x());
            theta = seg_angle + pi_over_two;
          }


          else if (   it->vertex1()
               && it->next_vd()->is_secondary()
               && it->next_vd()->is_finite()
               && it->next_vd()->next_vd()->is_primary()
               //&& it->next_vd()->next_vd()->is_infinite()
              ) {
            if (inf_theta_debug) printf("\nCASE S6\n");
            vertex_type * p1 = it->next_vd()->vertex0();
            vertex_type * p2 = it->next_vd()->vertex1();
            theta = atan2(p2->y() - p1->y(), p2->x() - p1->x());
          }
          else if (   it->vertex0()
               && it->prev_vd()->is_secondary()
               && it->prev_vd()->is_finite()
               && it->prev_vd()->prev_vd()->is_primary()
               //&& it->prev_vd()->prev_vd()->is_infinite()
              ) {
            if (inf_theta_debug) printf("\nCASE S7\n");
            vertex_type * p1 = it->prev_vd()->vertex0();
            vertex_type * p2 = it->prev_vd()->vertex1();
            theta = atan2(p2->y() - p1->y(), p2->x() - p1->x());
          }



          else {
            printf("failed to find tangent angle for infinite secondary edge\n");
            if (inf_theta_debug) {
              printf("v0[%d] v1[%d] N[%s %s], NN[%s %s] P[%s %s] PP[%s %s]\n",
              it->vertex0()?1:0,it->vertex1()?1:0,
              it->next()->is_primary()?"PRI":"2ND",it->next()->is_primary()?"INF":"FIN",
              it->next()->next()->is_primary()?"PRI":"2ND",it->next()->next()->is_primary()?"INF":"FIN",
              it->prev()->is_primary()?"PRI":"2ND",it->prev()->is_primary()?"INF":"FIN",
              it->prev()->prev()->is_primary()?"PRI":"2ND",it->prev()->prev()->is_primary()?"INF":"FIN"
              );
            }
            continue;
          }

        while (theta >  pi) { theta -= two_pi; }
        while (theta < -pi) { theta += two_pi; }

        }

      if (inf_theta_debug) printf("\n");
      
      it->theta(theta);

      }
    }
    
    // calculate phi
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {
      
      double phi = 0.0;
      
      bool dbgphi = (false
               && it->vertex0() && it->vertex1()
               && (it->vertex0()->x() - 1310000) < 10000
               && (it->vertex0()->x() - 1310000) > 3500
               && (it->vertex1()->x() - 1310000) < 10000
               && (it->vertex1()->x() - 1310000) > 3500
              );

      if (it->is_finite()) {
      
        double tdx;
        double tdy;
        
        if (dbgphi) printf("4 [%f %f %f]\n  [%f %f %f]\n  [%f %f]\n  the: %f \n    phi: ",it->vertex0()->x(),it->vertex0()->y(),it->vertex0()->r(),it->vertex1()->x(),it->vertex1()->y(),it->vertex1()->r(),it->foot()->x(),it->foot()->y(),it->theta());
        
        if (it->prev() == it->twin()) { // medial axis corner case

          if (dbgphi) printf("CRN ");

          tdx = (double) it->next()->foot()->x() - it->next()->vertex0()->x();
          tdy = (double) it->next()->foot()->y() - it->next()->vertex0()->y();
        }
        else {
          
          if (dbgphi) printf("NRM ");
          
          tdx = (double) it->foot()->x() - it->vertex0()->x();
          tdy = (double) it->foot()->y() - it->vertex0()->y();
        }

        if (tdx == 0 && tdy == 0) {
          phi = it->theta();
          if (dbgphi) printf("00 ");
        }
        else {
          double fpa = atan2(tdy, tdx);
      while (fpa >=  pi) {fpa-=2*pi;}
      while (fpa < -pi) {fpa+=2*pi;}
          phi = fpa - it->theta();

          if (dbgphi) printf("%.17f - %f = %f",fpa,it->theta(), phi);
        }

      }
      else { // is infinite
        phi = it->theta() + pi_over_two;

        if (dbgphi) printf("INF %f + pi/2 = %f ",it->theta(), phi);      
      }

      while (phi >  pi) {phi-=2*pi;}
      while (phi < -pi) {phi+=2*pi;}

      it->phi(phi);

      if (dbgphi) printf(", %f \n",phi);
    }


    // figure Q for curved segments, now that all thetas are available
    // Q is the middle control point for a quadratic Bezier curve.
    for (edge_iterator it = edges_.begin(); it != edges_.end(); ++it) {

      if (it->is_curved()) {
        // Based on solving for x1 of the three Bezier contol points P0, P1, P2,
        // where P0 and P2 are the edge's known endpoints, and we know the
        // angle or slope at those endpoints, which is enough to solve for P1:(x1,y1).
        // (Q is P1)
        detail::point_2d<default_voronoi_builder::int_type> q;
        double theta_start = it->theta();
        double theta_end = it->twin()->theta();
        theta_end += pi;
        while (theta_end >  pi) {theta_end = theta_end - (pi * 2);}
        double tan_theta_start = tan(theta_start);
        double tan_theta_end = tan(theta_end);
        double delta_tan = tan_theta_end - tan_theta_start;
        if ( delta_tan == 0.0 ) {
          // The curved edge must be pretty short and insignificant.
          // Too short to do our math anyway.
          // Intuition suggests the control point might not generally be 
          // tending toward the center point of the line.
          // But everything's so small here that that should be good enough.
          // We don't have enough precision to do any better, probably.
          q.x(it->vertex0()->x()/2 + it->vertex1()->x()/2);
          q.y(it->vertex0()->y()/2 + it->vertex1()->y()/2);
        } else {
          double y0 = (double) it->vertex0()->y();
          double y2 = (double) it->vertex1()->y();
          double x0 = (double) it->vertex0()->x();
          double x2 = (double) it->vertex1()->x();
          double delta_y0_y2 = y0 - y2;
          double delta_x2_x0 = x2 - x0;
          double x1 = x0 + ( ( delta_y0_y2 + tan_theta_end * delta_x2_x0 ) / delta_tan );
          double y1 = y0 + tan(theta_start) * (x1 - x0);
          q.x( (int) x1 );
          q.y( (int) y1 );
        }

        it->Q( (int) q.x(), (int) q.y() );
        
      }
    }

  }

 private:
  typedef typename cell_container_type::iterator cell_iterator;
  typedef typename vertex_container_type::iterator vertex_iterator;
  typedef typename edge_container_type::iterator edge_iterator;
  typedef typename TRAITS::vertex_equality_predicate_type
    vertex_equality_predicate_type;

  template <typename SEvent>
  bool is_primary_edge(const SEvent& site1, const SEvent& site2) const {
    bool flag1 = site1.is_segment();
    bool flag2 = site2.is_segment();
    if (flag1 && !flag2) {
      return (site1.point0() != site2.point0()) &&
             (site1.point1() != site2.point0());
    }
    if (!flag1 && flag2) {
      return (site2.point0() != site1.point0()) &&
             (site2.point1() != site1.point0());
    }
    return true;
  }

  template <typename SEvent>
  bool is_linear_edge(const SEvent& site1, const SEvent& site2) const {
    if (!is_primary_edge(site1, site2)) {
      return true;
    }
    return !(site1.is_segment() ^ site2.is_segment());
  }

  // Remove degenerate edge.
  void remove_edge(edge_type* edge) {
    
    // Are these two ifs necessary?
    // Put these in for debugging, where the problem was something else,
    // but these do fill in/transfer some missing feet.
    // After revising the foot-finding (trying to do it all in the sweepline
    // event processing), see if these are still needed.
    // ... also consider that this will happily put a primary's foot onto
    // a following secondary - not what you meant to do here - this was probably
    // to get consecutive primaries' foot transfered. Definitely want to nix this.
    // okay then...
    /*
    if (edge->foot() && !edge->next()->foot()) {
      edge->next()->foot(edge->foot()->x(), edge->foot()->y());
    }
    if (edge->twin()->foot() && !edge->twin()->next()->foot()) {
      edge->twin()->next()->foot(edge->twin()->foot()->x(), edge->twin()->foot()->y());
    }
    */

    // Update the endpoints of the incident edges to the second vertex.
    vertex_type* vertex = edge->vertex0();
    edge_type* updated_edge = edge->twin()->rot_next();
    while (updated_edge != edge->twin()) {
      updated_edge->vertex0(vertex);
      updated_edge = updated_edge->rot_next();
    }
    
    edge_type* edge1 = edge;
    edge_type* edge2 = edge->twin();

    edge_type* edge1_rot_prev = edge1->rot_prev();
    edge_type* edge1_rot_next = edge1->rot_next();

    edge_type* edge2_rot_prev = edge2->rot_prev();
    edge_type* edge2_rot_next = edge2->rot_next();

    // Update prev/next pointers for the incident edges.
    edge1_rot_next->twin()->next(edge2_rot_prev);
    edge2_rot_prev->prev(edge1_rot_next->twin());
    edge1_rot_prev->prev(edge2_rot_next->twin());
    edge2_rot_next->twin()->next(edge1_rot_prev);

  }

  void mark_side(bool inside, edge_type* edge, std::queue<vertex_type*> &boundary_vertices) {

    //assert(boundary_vertices.size() < 300);
    bool mark_debug = false;

    // halt 1, if already marked
    if (edge->is_external() || edge->is_internal()) {
      if (mark_debug) printf("[%llu] : HALT %s%s \n",(long long unsigned) edge, edge->is_internal()?"1":"0", edge->is_external()?"1":"0");
      return;
    }

    // set internal/external status of this edge
    if (inside) {
      edge->is_internal(true);
    } else {
      edge->is_external(true);
    }

    if (mark_debug) printf("[%llu] : edge did \n[%f %f %f] [%f %f %f]\n",
      (long long unsigned) edge,
      edge->vertex0()->x(),edge->vertex0()->y(),edge->vertex0()->r(),
      edge->vertex1()->x(),edge->vertex1()->y(),edge->vertex1()->r());

    if (edge->vertex0()->r() == 0 || 
        (edge->vertex1()->r() == 0)
        ) {
        if (inside) {
          edge->twin()->is_internal(true);
        } else {
          edge->twin()->is_external(true);
        }

        if (mark_debug) printf("[%llu] : twin didC\n[%f %f %f] [%f %f %f]\n",
          (long long unsigned) edge,
          edge->twin()->vertex0()->x(),edge->twin()->vertex0()->y(),edge->twin()->vertex0()->r(),
          edge->twin()->vertex1()->x(),edge->twin()->vertex1()->y(),edge->twin()->vertex1()->r());
    }

    if (edge->vertex1()->r() == 0 
        && (   (edge->is_secondary() && edge->next()->is_primary() )
            || (edge->is_primary()   && edge->next()->is_secondary())
           )
        ) {
      boundary_vertices.push(edge->vertex1());
      if (mark_debug)  {
        printf("[%llu] : HALT ON r==0, 2nd to prim or prim to 2nd; boundary verts size %lu \n",
               (long long unsigned) edge, boundary_vertices.size()
        ); 
      }
      return;
    }
    // follow chain of next() edges
    else {

      edge_type * e = edge->next();
        
      while (!e->is_internal() && !e->is_external() && e != edge) {
        if (inside) {
          e->is_internal(true);
        }
        else {
          e->is_external(true);
        }
        if (mark_debug) printf("[%llu] : EDGE NEXT SET %s \n",(long long unsigned) e, inside ? "INT":"EXT");


        boundary_vertices.push(e->vertex0());


        if (e->vertex1()->r() == 0 
            //&& e->is_secondary()
           ) {
          if (inside) {
            e->twin()->is_internal(true);
          } else {
            e->twin()->is_external(true);
          }
    
          boundary_vertices.push(edge->vertex1()); // eh? right?

          if (mark_debug) printf("[%llu] : twin didC2\n[%f %f %f] [%f %f %f]\n",
            (long long unsigned) e,
            edge->twin()->vertex0()->x(),edge->twin()->vertex0()->y(),edge->twin()->vertex0()->r(),
            edge->twin()->vertex1()->x(),edge->twin()->vertex1()->y(),edge->twin()->vertex1()->r());
          break;
        }

        e = e->next();

      }
    }
  }

  void mark_exterior(edge_type* edge) {
    //if (edge->color() == 1) {
    if (edge->is_external() || edge->is_internal()) {
      //printf("[%lu] already marked, stop recursion",(unsigned long) edge);
      return;
    }

    _mark_exterior(edge);
    
    vertex_type* v = edge->vertex1();
    if (!v) {
      v = edge->vertex0();
    }

    if (
        v == NULL 
       ) {
      return;
    }
    if (
        edge->is_secondary() 
       ) {
      return;
    }

    // Detect v->r() == 0. That's a more reliable general case indicator
    // of boundary crossing. Than the rot_next() cycle edge pattern we
    // looked for below, after seting e_spin.
    // Consider the center point of a bowtie as
    // a boundary case the spin approach will get wrong but this r==0 gets right.

    if (v->r() == 0) return;

    // Do non-recursive marking on secondary edges around this vertex,
    // as long as there aren't two rotate-adjacent secondaries, which is the 
    // hallmark of crossing through a boundary where two input segments meet,
    // which is what we don't want to do for our medial axis representation.

    edge_type* e_spin = v->incident_edge();

    /*
    if (e_spin == edge || e_spin == edge->twin()) {
        e_spin = e_spin->rot_next();
        if (e_spin == edge || e_spin == edge->twin()) {
          return;
        }
    }
    */

    do {
      if (   e_spin->is_secondary() 
          && e_spin->rot_next()->is_primary()
          && e_spin->rot_prev()->is_primary()
         ) {
        _mark_exterior(e_spin);
      }
        e_spin = e_spin->rot_next();
    } while(e_spin != v->incident_edge());

    // do more recursive marking just on primary edges
    edge_type* e = v->incident_edge();

    if (e == edge || e == edge->twin()) {
      e=e->rot_next();
      if (e == edge || e == edge->twin()) {
        return;
      }
    }

    do {
      if (e->is_primary()) {
        mark_exterior(e);
      }
      e = e->rot_next();
    } while (e != v->incident_edge());
  }

  void _mark_exterior(edge_type* edge) {
    edge->is_external(true);
    edge->twin()->is_external(true);
    //printf("  marked [%lu] [%lu]\n",(unsigned long)edge,(unsigned long)edge->twin());
  }

  void rotate_2d(double &x, double &y, const double theta, const double xo = 0, const double yo = 0) {
    double xp;
    x -= xo;
    y -= yo;
    xp = (x * cos(theta) - y * sin(theta)) + xo;
    y  = (y * cos(theta) + x * sin(theta)) + yo;
    x  = xp;
  }
  template <typename CT>
  void reflect(CT &x, CT &y, const CT x0, const CT y0, const CT x1, const CT y1) {
    double dy = (double) (y1 - y0);
    double dx = (double) (x1 - x0);
    if (dy == 0 && dx == 0) {return;}
    double theta = atan2(dy, dx);
    rotate_2d(x, y, -theta, x0, y0);
    y -= y0;
    y *= -1.0;
    y += y0;
    rotate_2d(x, y, theta, x0, y0);
  }

void makefoot(double & x, double & y, const double x0, const double y0,
                                      const double x1, const double y1) {
    // infinite slope case first
    if (x1 - x0 == 0) {
        x = x0;
    } else {
      double m  = (y1 - y0)/(x1 - x0);
      if (m == 0) {
          y = y0;
      }
      else {
        double intersect_x = ((m * x0) - y0 + ((1 / m) * x) + (y)) / (m + (1 / m));
          double intersect_y = -(x0 - intersect_x) * m + y0;
            x = intersect_x;
            y = intersect_y;
      }
    }
  }

  cell_container_type cells_;
  vertex_container_type vertices_;
  edge_container_type edges_;
  mutable std::string event_log_;
  vertex_equality_predicate_type vertex_equality_predicate_;

  // Disallow copy constructor and operator=
  medial_axis(const medial_axis&);
  void operator=(const medial_axis&);
};
}  // polygon
}  // boost

#endif  // BOOST_POLYGON_MEDIAL_AXIS
