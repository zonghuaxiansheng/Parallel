#include "mpi.h"
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

struct Rank2ColorKey {
  int rank_;
  int color_;
  int key_;
};

int main(int args, char** argv) {

  int rank;
  int size;

  // Init MPI
  MPI_Init(&args, &argv);
  // First group
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  char host_name[128];
  gethostname(host_name, 128);

  char* all_name = new char[size*sizeof(char)*128];
  // Collect all host name
  MPI_Gather(host_name, 128, MPI_CHAR, all_name, size*128, MPI_CHAR, 0, MPI_COMM_WORLD);

  Rank2ColorKey* rck_ptr = new Rank2ColorKey[size*sizeof(Rank2ColorKey)];
  Rank2ColorKey rck;
  MPI_Datatype new_type;
  MPI_Datatype type[3] = {MPI_INT, MPI_INT, MPI_INT};
  int len[3] = {1, 1, 1};
  MPI_Aint disp[3];
  disp[0] = &rck_ptr[0].rank_ - &rck_ptr[0];
  disp[0] = &rck_ptr[0].color_ - &rck_ptr[0];
  disp[0] = &rck_ptr[0].key_ - &rck_ptr[0];
  MPI_Type_create_struct(3, len, disp, type, &new_type);

  if (rank == 0) {
    std::map<std::string, std::pair<int, int>> name_map;
    int sub_color = 0;
    for (int i = 0; i < size; i ++) {
      rck_ptr[i].rank_ = i;
      char name_[128];
      strcpy(name_, all_name + size*128);
      std::string name_str(name_);
      std::map<std::string, int>::iterator iter = name_map.find(name_str);
      if (iter != name_map.end()) {
        rck_ptr[i].color_ = iter->second.first;
        rck_ptr[i].key_ = iter->second.second ++;
      } else {
        rck_ptr[i].color_ = sub_color;
        rck_ptr[i].key_ = 0;
        name_map.insert(std::make_pair(name_str, std::make_pair(sub_color ++, 0)));
      }
      // std::cout << "Host name: " << name_ << std::endl;
    }
  }

  MPI_Scatter(rck_ptr, size, new_type, rck, 1, new_type, 0, MPI_COMM_WORLD);

  int cltr_rank = rck.rank_;
  int cltr_color = rck.color_;
  int cltr_key = rck.key_;

  if (cltr_rank != rank) {
    std::cerr << "WTF ?" << std::endl;
    exit(1);
  }

  MPI_Comm L1_World, L2_World, Sub1_World, Sub2_World;
  MPI_Comm_dup(MPI_COMM_WORLD, &L1_World);
  MPI_Comm_dup(MPI_COMM_WORLD, &L2_World);

  MPI_Comm_split(L1_World, cltr_color, cltr_key, &Sub1_World);

  if (cltr_key == 0) {
    MPI_Comm_split(L2_World, cltr_key, cltr_color, &Sub2_World);
  }

  int ibuf;
  if (rank == 0) {
    ibuf = 12345;
  } else {
    ibuf = 0;
  }

  MPI_Bcast(&ibuf, 1, MPI_INT, 0, &Sub2_World);

  MPI_Bcast(&ibuf, 1, MPI_INT, cltr_key, &Sub1_World);

  // Finish
  MPI_Finalize();
  return 0;
}
