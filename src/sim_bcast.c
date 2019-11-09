#include "mpi.h"
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

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
  
  std::map<std::string, int> name_map;
  if (rank == 0) {
    for (int i = 0; i < size; i ++) {
      char name_[128];
      strcpy(name_, all_name + size*128);
      std::string name_str(name_);
      std::map<std::string, int>::iterator iter = name_map.find(name_str);
      if (iter != name_map.end()) {
        iter->second ++;
      } else {
        name_map.insert(std::make_pair(name_str, 1));
      }
      // std::cout << "Host name: " << name_ << std::endl;
    }

    for (auto item : name_map) {
      std::cout << "Host name: " << item.first << " Count: " << item.second << std::endl;
    }
  }

  // delete all_name;

  // Finish
  MPI_Finalize();
  return 0;
}
