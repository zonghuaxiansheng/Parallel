#include <cstring>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "mpi.h"
#include "experiment.h"

namespace ustc_parallel {

	struct Node
	{
		int key_;
		int color_;
	};

	void inline gethostname(int my_rank, char* hostname, int len) {
		char tname[16];
		int group = my_rank % 3;
		if (group == 0) {
			strcpy(tname, "Node0");
		} else if (group == 1) {
			strcpy(tname, "Node1");
		} else {
			strcpy(tname, "Node2");
		}
		strcpy(hostname, tname);
	}

	void CreatePipeLine(int& my_rank, int& psize, MPI_Comm my_comm) {
		MPI_Request request_r, request_s;
		MPI_Status status_r, status_s;
		int pre_rank = my_rank - 1;
		int next_rank = my_rank + 1;
		int tag = 10;
		char next_rbuf[128];
		char cur_rbuf[128];
		char pre_sbuf[128] = "A";
		char cur_sbuf[128] = "A";
		int buf_len = 1;
		int rcnt = 0;

		while (rcnt < 10) {
			if (my_rank == 0) {
				// MPI_Send(sbuf, 128, MPI_CHAR, next_rank, tag, my_comm);
				MPI_Isend(pre_sbuf, 128, MPI_CHAR, next_rank, tag, my_comm, &request_s);
				// MPI_Wait(&request_s, &status_s);
				std::cout << "Rank: " << my_rank
					<< " Sbuf: " << pre_sbuf << " Dst Rank: " << next_rank << std::endl;
			} else if (my_rank == psize - 1) {
				MPI_Irecv(next_rbuf, 128, MPI_CHAR, pre_rank, tag, my_comm, &request_r);
				MPI_Wait(&request_r, &status_r);
				std::cout << "Rank: " << my_rank
					<< " Rbuf: " << next_rbuf << " Src Rank: " << status_r.MPI_SOURCE << std::endl;
			} else {
				MPI_Irecv(next_rbuf, 128, MPI_CHAR, pre_rank, tag, my_comm, &request_r);
				MPI_Isend(pre_sbuf, 128, MPI_CHAR, next_rank, tag, my_comm, &request_s);
				MPI_Wait(&request_r, &status_r);
				// MPI_Wait(&request_s, &status_s);
				std::cout << "Rank: " << my_rank
					<< " Rbuf: " << next_rbuf << " Src Rank: " << status_r.MPI_SOURCE << std::endl;
				strcpy(cur_rbuf, next_rbuf);
				for (int i = 0; i < buf_len; i++) {
					cur_sbuf[i] = cur_rbuf[i] + 1;
				}
				cur_sbuf[buf_len] = '\0';
				strcpy(pre_sbuf, cur_sbuf);
				// MPI_Send(sbuf, 128, MPI_CHAR, next_rank, tag, my_comm);
				std::cout << "Rank: " << my_rank
					<< " Sbuf: " << pre_sbuf << " Dst Rank: " << next_rank << std::endl;
			}
			rcnt ++;
		}
	}

	void CreateSimAlltoAll(int& my_rank, int& psize, MPI_Comm my_comm) {
		MPI_Status status;
		int tag = 10;
		char sbuf[16] = "abc";
		char* rbuf = new char[psize*16];
		memset(rbuf, '\0', psize * 16);

		for (int i = 0; i < psize; i++) {
			if (my_rank != i) {
				MPI_Send(sbuf, 16, MPI_CHAR, i, tag, my_comm);
			}
		}

		for (int i = 0; i < psize; i++) {
			if (my_rank != i) {
				MPI_Recv(rbuf + i * 16, 16, MPI_CHAR, i, tag, my_comm, &status);
			}
		}

		strcpy(rbuf + my_rank * 16, sbuf);

		for (int i = 0; i < psize; i++) {
			char tbuf[16];
			strcpy(tbuf, rbuf + i * 16);
			std::cout << "Rank: " << my_rank
				<< " Buf[" << i << "]: " << tbuf << std::endl;
		}

		delete [] rbuf;
	}

	void CreateSimBcast(int& my_rank, int& psize, MPI_Comm my_comm) {
	
		char hostname[16];
		gethostname(my_rank, hostname, 16);

		char* allname = new char[psize * 16];
		Node* nodes = new Node[psize];
		memset(allname, '\0', psize * 16);

		MPI_Gather(hostname, 16, MPI_CHAR, allname, 16, MPI_CHAR, 0, MPI_COMM_WORLD);

		int sub_color = 0, sub_key = 0;
		int color_tag = 10, key_tag = 20;

		if (my_rank == 0) {
			std::map<std::string, std::pair<int, int>> name_map;
			int sub_color = 0;
			for (int i = 0; i < psize; i++) {
				char tname[16];
				strcpy(tname, allname + i * 16);
				std::string name_str(tname);
				// std::cout << tname << std::endl;
				std::map<std::string, std::pair<int, int>>::iterator iter = name_map.find(name_str);
				if (iter != name_map.end()) {
					nodes[i].color_ = iter->second.first;
					nodes[i].key_ = ++ iter->second.second;
				}
				else {
					nodes[i].color_ = sub_color;
					nodes[i].key_ = 0;
					name_map.insert(std::make_pair(name_str, std::make_pair(sub_color++, 0)));
				}
			}
			for (int i = 0; i < psize; i++) {
				// std::cout << "Node: " << i << " " << nodes[i].color_ << " " << nodes[i].key_ << std::endl;
				if (my_rank != i) {
					MPI_Send(&nodes[i].color_, 1, MPI_INT, i, color_tag, my_comm);
					MPI_Send(&nodes[i].key_, 1, MPI_INT, i, key_tag, my_comm);
				}
			}
			sub_color = nodes[0].color_;
			sub_key = nodes[0].key_;
		} else {
			MPI_Status status;
			MPI_Recv(&sub_color, 1, MPI_INT, 0, color_tag, my_comm, &status);
			MPI_Recv(&sub_key, 1, MPI_INT, 0, key_tag, my_comm, &status);
		}

		std::cout << "Rank: " << my_rank << " Sub Color: " << sub_color << " Sub Key: " << sub_key << std::endl;

		MPI_Barrier(my_comm);

		MPI_Comm leader_comm = MPI_COMM_WORLD, inside_comm = MPI_COMM_WORLD;
		MPI_Comm_split(my_comm, sub_color, sub_key, &inside_comm);
		MPI_Comm_split(my_comm, sub_key, sub_color, &leader_comm);
		
		char sbuf[16];
		if (my_rank == 0) {
			strcpy(sbuf, "sim_bcast");
		} else {
			strcpy(sbuf, "");
		}
		
		std::cout << "Before Split -> Rank: " << my_rank << " Sbuf: " << sbuf << std::endl;
		MPI_Barrier(my_comm);

		MPI_Bcast(sbuf, 16, MPI_CHAR, 0, leader_comm);
		std::cout << "Leader Split -> Rank: " << my_rank << " Sbuf: " << sbuf << std::endl;
		MPI_Barrier(my_comm);

		MPI_Bcast(sbuf, 16, MPI_CHAR, 0, inside_comm);
		std::cout << "Inside Split -> Rank: " << my_rank << " Sbuf: " << sbuf << std::endl;
	}

}