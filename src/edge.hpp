#ifndef EDGE_HPP
#define EDGE_HPP

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

#include "util/string.hpp"

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

struct Edge{
  enum direction{ horizontal, vertical };
  int source_site;
  int source_leg;
  int target_site;
  int target_leg;
  direction dir;
  int op_id;
  Edge():source_site(-1), source_leg(-1), target_site(-1), target_leg(-1), dir(horizontal), op_id(-1){}

  /*
   * horizontal:
   *
   *  s-t
   *
   * vertical:
   *
   *  t
   *  |
   *  s
   */
  Edge(int source_site, int target_site, direction dir, int op_id)
      : source_site(source_site), source_leg(dir==horizontal?2:1),
        target_site(target_site), target_leg(dir==horizontal?0:3),
        dir(dir), op_id(op_id){}

  bool is_horizontal() const { return dir==horizontal; }
  bool is_vertical() const { return !is_horizontal(); }
};

typedef std::vector<Edge> Edges;


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
Edges make_edges(std::string const& str){
  constexpr static int NFIELD = 4;
  Edges edges;
  std::stringstream ss(str);
  std::string line;
  while(std::getline(ss, line)){
    line = util::strip(util::drop_comment(line));
    if(line.empty()){continue;}
    auto words = util::split(line);
    if(words.size() != NFIELD){
      std::stringstream msg("");
      msg << words.size() << "fields are found in \"" << line << "\" (required " << NFIELD << " fields.)";
      throw std::invalid_argument(msg.str());
    }
    int source = std::stoi(words[0]);
    int target = std::stoi(words[1]);
    Edge::direction dir = Edge::horizontal;
    switch(words[2][0]){
      case 'h':
        dir = Edge::horizontal;
        break;
      case 'v':
        dir = Edge::vertical;
        break;
      default:
        std::stringstream msg("the third arg in ");
        msg << line << " should be h or v .";
        throw std::invalid_argument(msg.str());
    }
    int opid = std::stoi(words[3]);
    edges.push_back(Edge(source, target, dir, opid));
  }
  return edges;
}

#endif // EDGE_HPP
