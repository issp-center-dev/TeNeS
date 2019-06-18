#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <string>
#include <fstream>
#include <vector>
// Lattice setting

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
 *
 */
class Lattice {
public:
  int LX;
  int LY;
  int N_UNIT;

  std::vector<std::vector<int> > Tensor_list;
  std::vector<std::vector<int> > NN_Tensor;

  Lattice(int X=2, int Y=2) : LX(X), LY(Y), N_UNIT(LX*LY) {
    assert(X > 0);
    assert(Y > 0);

    reset();
  }

  void reset(){
    N_UNIT = LX * LY;
    Tensor_list.assign(LX, std::vector<int>(LY));
    NN_Tensor.assign(N_UNIT, std::vector<int>(4));
    int i=0;
    for(int ix=0; ix<LX; ++ix){
      for(int iy=0; iy<LY; ++iy){
        Tensor_list[ix][iy] = i;
        ++i;
      }
    }
    i=0;
    for(int ix=0; ix<LX; ++ix){
      for(int iy=0; iy<LY; ++iy){
        NN_Tensor[i][0] = Tensor_list[(ix-1+LX)%LX][iy];
        NN_Tensor[i][1] = Tensor_list[ix          ][(iy+1)%LY];
        NN_Tensor[i][2] = Tensor_list[(ix+1)%LX][iy];
        NN_Tensor[i][3] = Tensor_list[ix][(iy-1+LY)%LY];
        ++i;
      }
    }
  }
  void reset(int X, int Y){
    assert(X > 0);
    assert(Y > 0);
    LX = X;
    LY = Y;
    reset();
  }

  void read_parameters(const char *filename) {
    std::ifstream input_file;
    input_file.open(filename, std::ios::in);
    std::string reading_line_buffer;

    while (!input_file.eof()) {
      std::getline(input_file, reading_line_buffer);
      std::stringstream buf(reading_line_buffer);
      std::vector<std::string> result;
      while (buf >> reading_line_buffer) {
        result.push_back(reading_line_buffer);
      }

      if (result.size() > 1) {
        if (result[0].compare("LX") == 0) {
          std::istringstream is(result[1]);
          is >> LX;
        } else if (result[0].compare("LY") == 0) {
          std::istringstream is(result[1]);
          is >> LY;
        } else if (result[0].compare("N_UNIT") == 0) {
          std::istringstream is(result[1]);
          is >> N_UNIT;
        }
        // std::cout<< "## input data: "<<result[0]<<" =
        // "<<result[1]<<std::endl;
      }
    }
    reset();
  }
  void output_parameters(const char *filename, bool append) {
    std::ofstream ofs;
    if (append) {
      ofs.open(filename, std::ios::out | std::ios::app);
    } else {
      ofs.open(filename, std::ios::out);
    }

    ofs << "LX " << LX << std::endl;
    ofs << "LY " << LY << std::endl;
    ofs << "N_UNIT " << N_UNIT << std::endl;
  }

  void output_parameters(const char *filename) {
    output_parameters(filename, false);
  }

  void output_parameters_append(const char *filename) {
    output_parameters(filename, true);
  }

  void Bcast_parameters(MPI_Comm comm) {
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    std::vector<int> params_int(3);

    if (irank == 0) {
      params_int[0] = LX;
      params_int[1] = LY;
      params_int[2] = N_UNIT;

      MPI_Bcast(&params_int.front(), 3, MPI_INT, 0, comm);
    } else {
      MPI_Bcast(&params_int.front(), 3, MPI_INT, 0, comm);

      LX = params_int[0];
      LY = params_int[1];
      N_UNIT = params_int[2];
    }
    reset();
  }
};

/*

class Lattice_skew : public Lattice {
public:
  int LX_ori;
  int LY_ori;
  int LX_diff;
  int LY_diff;

  std::vector<std::vector<int> > Tensor_position;

  Lattice_skew() {
    LX_ori = 2;
    LY_ori = 2;
    LX = LX_ori;
    LY = LY_ori;
    LX_diff = LX_ori;
    LY_diff = LY_ori;
    N_UNIT = LX_ori * LY_ori;

    Tensor_position =
        std::vector<std::vector<int> >(N_UNIT, std::vector<int>(2));
    Tensor_list = std::vector<std::vector<int> >(LX, std::vector<int>(LY));
    NN_Tensor = std::vector<std::vector<int> >(N_UNIT, std::vector<int>(4));
  }

  void read_parameters(const char *filename) {
    std::ifstream input_file;
    input_file.open(filename, std::ios::in);
    std::string reading_line_buffer;

    while (!input_file.eof()) {
      std::getline(input_file, reading_line_buffer);
      std::stringstream buf(reading_line_buffer);
      std::vector<std::string> result;
      while (buf >> reading_line_buffer) {
        result.push_back(reading_line_buffer);
      }

      if (result.size() > 1) {
        if (result[0].compare("LX_ori") == 0) {
          std::istringstream is(result[1]);
          is >> LX_ori;
        } else if (result[0].compare("LY_ori") == 0) {
          std::istringstream is(result[1]);
          is >> LY_ori;
        } else if (result[0].compare("LX") == 0) {
          std::istringstream is(result[1]);
          is >> LX;
        } else if (result[0].compare("LY") == 0) {
          std::istringstream is(result[1]);
          is >> LY;
        } else if (result[0].compare("LX_diff") == 0) {
          std::istringstream is(result[1]);
          is >> LX_diff;
        } else if (result[0].compare("LY_diff") == 0) {
          std::istringstream is(result[1]);
          is >> LY_diff;
        } else if (result[0].compare("N_UNIT") == 0) {
          std::istringstream is(result[1]);
          is >> N_UNIT;
        }
        // std::cout<< "## input data: "<<result[0]<<" =
        // "<<result[1]<<std::endl;
      }
    }
  }
  void output_parameters(const char *filename, bool append) {
    std::ofstream ofs;
    if (append) {
      ofs.open(filename, std::ios::out | std::ios::app);
    } else {
      ofs.open(filename, std::ios::out);
    }

    ofs << "LX_ori " << LX_ori << std::endl;
    ofs << "LY_ori " << LY_ori << std::endl;
    ofs << "LX " << LX << std::endl;
    ofs << "LY " << LY << std::endl;
    ofs << "LX_diff " << LX_diff << std::endl;
    ofs << "LY_diff " << LY_diff << std::endl;
    ofs << "N_UNIT " << N_UNIT << std::endl;
  }

  void output_parameters(const char *filename) {
    output_parameters(filename, false);
  }

  void output_parameters_append(const char *filename) {
    output_parameters(filename, true);
  }

  void Bcast_parameters(MPI_Comm comm) {
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    std::vector<int> params_int(7);

    if (irank == 0) {
      params_int[0] = LX_ori;
      params_int[1] = LY_ori;
      params_int[2] = LX;
      params_int[3] = LY;
      params_int[4] = LX_diff;
      params_int[5] = LY_diff;
      params_int[6] = N_UNIT;

      MPI_Bcast(&params_int.front(), 7, MPI_INT, 0, comm);
    } else {
      MPI_Bcast(&params_int.front(), 7, MPI_INT, 0, comm);

      LX_ori = params_int[0];
      LY_ori = params_int[1];
      LX = params_int[2];
      LY = params_int[3];
      LX_diff = params_int[4];
      LY_diff = params_int[5];
      N_UNIT = params_int[6];
    }
  }
};

*/

#endif // LATTICE_HPP
