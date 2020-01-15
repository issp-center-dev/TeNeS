#ifndef EDGE_HPP
#define EDGE_HPP

#include <string>
#include <vector>

namespace tenes {

/*
 * axis:
 *
 *  y
 *  ^
 *  |
 *  .->x
 *
 * edge index:
 *
 *   1
 *  0.2
 *   3
 */

struct Edge {
  enum direction { horizontal, vertical };
  int source_site;
  int source_leg;
  int target_site;
  int target_leg;
  direction dir;
  int op_id;
  Edge()
      : source_site(-1), source_leg(-1), target_site(-1), target_leg(-1),
        dir(horizontal), op_id(-1) {}

  /*
   * horizontal:
   *
   *  s-t
   *
   * vertical:
   *
   *  s
   *  |
   *  t
   */
  Edge(int source_site, int target_site, direction dir, int op_id)
      : source_site(source_site), source_leg(dir == horizontal ? 2 : 3),
        target_site(target_site), target_leg(dir == horizontal ? 0 : 1),
        dir(dir), op_id(op_id) {}

  bool is_horizontal() const { return dir == horizontal; }
  bool is_vertical() const { return !is_horizontal(); }
};

using Edges = std::vector<Edge>;

/*
 * Input
 * """
 * 0 1 h 0
 * 0 2 v 0
 * """
 *
 * will be {Edge(0, 1, Edge::horizontal, 0),
 *          Edge(0, 2, Edge::vertival, 0)}
 *
 */
Edges make_edges(std::string const &str);

} // end of namespace tenes

#endif // EDGE_HPP
